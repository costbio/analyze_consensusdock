#!/usr/bin/env python3
"""
Local AlphaFold Structure Extractor
Extracts AlphaFold models from local AlphaFold database archive for required structures only.

Features:
- Extracts only required structures (from required_structures.csv)
- Extracts structures directly from local tar archive
- No network requests - completely offline
- Parallel extraction using multiple threads
- Progress tracking and detailed logging
- Handles compressed PDB files (.pdb.gz)
- Supports all fragments (F1, F2, F23, etc.)
- Resume capability: Skips already extracted files

Prerequisites:
- required_structures.csv (generated by identify_required_structures.py)
- cavity_mapping.csv (fallback if required_structures.csv not available)
- Local AlphaFold database archive

Usage:
    python extract_alphafold_models.py

Output:
    - alphafold_structures/ folder with extracted AlphaFold structures
    - alphafold_mapping.csv: Updated mapping with structure paths
    - alphafold_extraction.log
"""

import os
import sys
import pandas as pd
import logging
import time
import tarfile
import gzip
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# --- Configuration ---
REQUIRED_STRUCTURES_CSV = "required_structures.csv"  # Input file from identify_required_structures.py
CAVITY_MAPPING_CSV = "cavity_mapping.csv"  # Fallback input file from extract_cavities.py
ALPHAFOLD_ARCHIVE = "/opt/data/alphafold2/UP000005640_9606_HUMAN_v4.tar"  # Local AF database
OUTPUT_AF_DIR = "alphafold_structures"     # Output directory for AlphaFold structures
OUTPUT_MAPPING_CSV = "alphafold_mapping.csv"  # Output mapping with structure paths
LOG_FILE = "alphafold_extraction.log"
NUM_THREADS = 60  # Number of extraction threads
DATABASE_VERSION = "v4"  # AlphaFold database version

# Set up logging
logging.basicConfig(filename=LOG_FILE, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logging.getLogger().addHandler(console_handler)

def extract_and_decompress_pdb(tar_file, member_name, output_path):
    """
    Extract a PDB file from tar archive and decompress if needed.
    
    Parameters
    ----------
    tar_file: tarfile.TarFile
        Open tar file object
    member_name: str
        Name of the member in the tar file
    output_path: str
        Local path to save decompressed file
        
    Returns
    -------
    bool: True if extraction successful, False otherwise
    """
    try:
        # Check if file already exists
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return True
            
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Extract member to memory
        member = tar_file.getmember(member_name)
        extracted_file = tar_file.extractfile(member)
        
        if extracted_file is None:
            logging.warning(f"Could not extract {member_name}")
            return False
        
        # If file is gzipped, decompress it
        if member_name.endswith('.gz'):
            with gzip.open(extracted_file, 'rb') as gz_file:
                with open(output_path, 'wb') as out_file:
                    shutil.copyfileobj(gz_file, out_file)
        else:
            # Direct copy if not compressed
            with open(output_path, 'wb') as out_file:
                shutil.copyfileobj(extracted_file, out_file)
        
        # Verify file was extracted and has content
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return True
        else:
            logging.warning(f"Extracted file is empty: {output_path}")
            return False
            
    except Exception as e:
        logging.warning(f"Extraction failed for {member_name}: {e}")
        return False

def find_structures_in_archive(tar_file, required_structures):
    """
    Find required structures in the tar archive.
    
    Parameters
    ----------
    tar_file: tarfile.TarFile
        Open tar file object
    required_structures: set
        Set of AlphaFold IDs to find (e.g., {'AF-Q96RW7-F1', 'AF-Q96RW7-F23'})
        
    Returns
    -------
    dict: Mapping from AlphaFold ID to tar member name
    """
    logging.info("Scanning archive for required structures...")
    found_structures = {}
    
    # Get all members that are PDB files
    for member in tqdm(tar_file.getmembers(), desc="Scanning archive"):
        if member.isfile() and member.name.endswith('.pdb.gz'):
            # Extract AlphaFold ID from filename
            # Format: AF-A0A024R1R8-F1-model_v4.pdb.gz
            filename = os.path.basename(member.name)
            if filename.startswith('AF-') and '-model_' in filename:
                # Extract everything before '-model_'
                alphafold_id = filename.split('-model_')[0]
                if alphafold_id in required_structures:
                    found_structures[alphafold_id] = member.name
    
    logging.info(f"Found {len(found_structures)}/{len(required_structures)} required structures in archive")
    return found_structures

def process_extraction_task(extraction_task):
    """
    Process a single AlphaFold extraction task.
    
    Parameters
    ----------
    extraction_task: dict
        Dictionary containing extraction information
        
    Returns
    -------
    dict: Result of the extraction
    """
    uniprot_id = extraction_task['uniprot_id']
    fragment_id = extraction_task['fragment_id']
    alphafold_id = extraction_task['alphafold_id']
    archive_path = extraction_task['archive_path']
    member_name = extraction_task['member_name']
    output_path = extraction_task['output_path']
    
    try:
        # Open tar file and extract
        with tarfile.open(archive_path, 'r') as tar:
            extraction_success = extract_and_decompress_pdb(tar, member_name, output_path)
        
        return {
            'success': extraction_success,
            'uniprot_id': uniprot_id,
            'fragment_id': fragment_id,
            'alphafold_id': alphafold_id,
            'output_path': output_path if extraction_success else None,
            'extracted': extraction_success,
            'error': None if extraction_success else 'Extraction failed'
        }
        
    except Exception as e:
        return {
            'success': False,
            'uniprot_id': uniprot_id,
            'fragment_id': fragment_id,
            'alphafold_id': alphafold_id,
            'output_path': None,
            'extracted': False,
            'error': str(e)
        }

def load_required_structures():
    """Load required structures from CSV file with fallback to cavity mapping."""
    # Try to load required structures first (preferred)
    if os.path.exists(REQUIRED_STRUCTURES_CSV):
        try:
            df = pd.read_csv(REQUIRED_STRUCTURES_CSV)
            required_columns = ['UniProt_ID', 'Fragment_ID']
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                logging.error(f"Required columns missing from required structures: {missing_columns}")
                logging.info("Falling back to cavity mapping...")
            else:
                logging.info(f"Loaded {len(df)} required structures from {REQUIRED_STRUCTURES_CSV}")
                return df
                
        except Exception as e:
            logging.error(f"Error loading required structures: {e}")
            logging.info("Falling back to cavity mapping...")
    
    # Fallback to cavity mapping
    if not os.path.exists(CAVITY_MAPPING_CSV):
        logging.error(f"Neither required structures nor cavity mapping file found")
        logging.error("Please run identify_required_structures.py or extract_cavities.py first.")
        sys.exit(1)
    
    try:
        df = pd.read_csv(CAVITY_MAPPING_CSV)
        required_columns = ['UniProt_ID', 'Fragment_ID']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            logging.error(f"Required columns missing from cavity mapping: {missing_columns}")
            sys.exit(1)
        
        logging.info(f"Loaded {len(df)} structures from cavity mapping (fallback)")
        logging.warning("Using all cavity structures - for efficiency, run identify_required_structures.py first")
        return df
        
    except Exception as e:
        logging.error(f"Error loading cavity mapping: {e}")
        sys.exit(1)

def prepare_extraction_tasks(structures_df, archive_path, output_dir):
    """
    Prepare AlphaFold extraction tasks for required structures only.
    
    Returns
    -------
    list: List of extraction tasks
    dict: Statistics about extractions to be performed
    """
    logging.info("Preparing AlphaFold extraction tasks for required structures...")
    
    # Get unique UniProt ID and Fragment ID combinations
    unique_combinations = structures_df[['UniProt_ID', 'Fragment_ID']].drop_duplicates()
    
    # Build set of required AlphaFold IDs
    required_structures = set()
    for _, row in unique_combinations.iterrows():
        alphafold_id = f"AF-{row['UniProt_ID']}-F{row['Fragment_ID']}"
        required_structures.add(alphafold_id)
    
    logging.info(f"Will extract {len(required_structures)} unique AlphaFold structures")
    
    # Scan archive to find available structures
    with tarfile.open(archive_path, 'r') as tar:
        found_structures = find_structures_in_archive(tar, required_structures)
    
    # Prepare extraction tasks and track statistics
    extraction_tasks = []
    skipped_existing = 0
    not_in_archive = 0
    
    for _, row in tqdm(unique_combinations.iterrows(), total=len(unique_combinations), desc="Preparing extraction tasks"):
        uniprot_id = row['UniProt_ID']
        fragment_id = row['Fragment_ID']
        alphafold_id = f"AF-{uniprot_id}-F{fragment_id}"
        
        # Check if structure is available in archive
        if alphafold_id not in found_structures:
            logging.debug(f"Structure not found in archive: {alphafold_id}")
            not_in_archive += 1
            continue
        
        # Generate local paths
        fragment_dir = os.path.join(output_dir, uniprot_id, f"F{fragment_id}")
        output_path = os.path.join(fragment_dir, f"{alphafold_id}-model_{DATABASE_VERSION}.pdb")
        
        # Check if file already exists (skip extraction if present)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            logging.debug(f"AlphaFold structure already exists for {uniprot_id} F{fragment_id}")
            skipped_existing += 1
            continue
        
        extraction_tasks.append({
            'uniprot_id': uniprot_id,
            'fragment_id': fragment_id,
            'alphafold_id': alphafold_id,
            'archive_path': archive_path,
            'member_name': found_structures[alphafold_id],
            'output_path': output_path
        })
    
    stats = {
        'total_combinations': len(unique_combinations),
        'unique_uniprot_ids': len(structures_df['UniProt_ID'].unique()),
        'found_in_archive': len(found_structures),
        'extraction_tasks': len(extraction_tasks),
        'already_extracted': skipped_existing,
        'not_in_archive': not_in_archive
    }
    
    logging.info(f"Found {stats['unique_uniprot_ids']} unique UniProt IDs with {stats['total_combinations']} fragments")
    logging.info(f"Available in archive: {stats['found_in_archive']} structures")
    logging.info(f"Already extracted (skipped): {stats['already_extracted']} structures")
    logging.info(f"Not in archive: {stats['not_in_archive']} structures")
    logging.info(f"Prepared {stats['extraction_tasks']} extraction tasks")
    
    return extraction_tasks, stats

def create_alphafold_mapping(structures_df, output_dir):
    """Create mapping CSV with AlphaFold structure paths."""
    logging.info("Creating AlphaFold mapping...")
    
    updated_rows = []
    for _, row in structures_df.iterrows():
        uniprot_id = row['UniProt_ID']
        fragment_id = row['Fragment_ID']
        pocket_pdb = row.get('Pocket_PDB', None)  # May not exist in required_structures.csv
        
        # Generate AlphaFold structure paths
        alphafold_id = f"AF-{uniprot_id}-F{fragment_id}"
        fragment_dir = os.path.join(output_dir, uniprot_id, f"F{fragment_id}")
        structure_path = os.path.join(fragment_dir, f"{alphafold_id}-model_{DATABASE_VERSION}.pdb")
        
        # Check if structure file exists
        structure_exists = os.path.exists(structure_path) and os.path.getsize(structure_path) > 0
        
        updated_rows.append({
            'UniProt_ID': uniprot_id,
            'Fragment_ID': fragment_id,
            'Cavity_Index': row.get('Cavity_Index', 1),
            'Receptor_PDB': structure_path if structure_exists else None,
            'Pocket_PDB': pocket_pdb,
            'AlphaFold_ID': alphafold_id,
            'AlphaFold_Available': structure_exists
        })
    
    updated_df = pd.DataFrame(updated_rows)
    return updated_df

def main():
    """Main function for AlphaFold model extraction."""
    start_time = time.time()
    
    logging.info("Starting AlphaFold model extraction workflow")
    logging.info(f"Input (preferred): {REQUIRED_STRUCTURES_CSV}")
    logging.info(f"Input (fallback): {CAVITY_MAPPING_CSV}")
    logging.info(f"Archive: {ALPHAFOLD_ARCHIVE}")
    logging.info(f"Output directory: {OUTPUT_AF_DIR}")
    logging.info(f"Output mapping: {OUTPUT_MAPPING_CSV}")
    logging.info(f"Number of threads: {NUM_THREADS}")
    logging.info(f"Database version: {DATABASE_VERSION}")
    
    # Check if archive exists
    if not os.path.exists(ALPHAFOLD_ARCHIVE):
        logging.error(f"AlphaFold archive not found: {ALPHAFOLD_ARCHIVE}")
        logging.error("Please check the path to the AlphaFold database archive.")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_AF_DIR, exist_ok=True)
    
    # Load required structures (with fallback to cavity mapping)
    structures_df = load_required_structures()
    
    # Prepare extraction tasks
    extraction_tasks, stats = prepare_extraction_tasks(structures_df, ALPHAFOLD_ARCHIVE, OUTPUT_AF_DIR)
    
    if not extraction_tasks:
        logging.info("No extraction tasks needed. All required AlphaFold structures are already available.")
        logging.info(f"Found {stats['already_extracted']} already extracted structures")
        logging.info(f"Not in archive: {stats['not_in_archive']} structures")
        if stats['already_extracted'] > 0:
            logging.info("Time saved by reusing existing AlphaFold structures: Significant!")
    else:
        # Perform parallel extractions
        logging.info(f"Starting parallel extraction of {len(extraction_tasks)} AlphaFold structures...")
        logging.info(f"Skipping {stats['already_extracted']} already extracted structures")
        
        successful_extractions = 0
        failed_extractions = 0
        
        with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
            # Submit all tasks
            future_to_task = {executor.submit(process_extraction_task, task): task for task in extraction_tasks}
            
            # Process results with progress bar
            for future in tqdm(as_completed(future_to_task), total=len(extraction_tasks), 
                             desc="Extracting AlphaFold structures"):
                result = future.result()
                
                if result['success']:
                    successful_extractions += 1
                    logging.debug(f"Extracted: {result['uniprot_id']} F{result['fragment_id']}")
                else:
                    failed_extractions += 1
                    logging.warning(f"Failed to extract {result['uniprot_id']} F{result['fragment_id']}: {result['error']}")
        
        logging.info(f"New extractions: {successful_extractions} successful, {failed_extractions} failed")
        logging.info(f"Total structures now available: {stats['already_extracted'] + successful_extractions}")
        logging.info(f"Time saved by skipping existing files: Significant!")
    
    # Create AlphaFold mapping CSV
    alphafold_df = create_alphafold_mapping(structures_df, OUTPUT_AF_DIR)
    alphafold_df.to_csv(OUTPUT_MAPPING_CSV, index=False)
    
    # Summary statistics
    total_structure_files = len([f for f in Path(OUTPUT_AF_DIR).rglob("*.pdb")])
    available_structures = len(alphafold_df[alphafold_df['AlphaFold_Available'] == True])
    
    elapsed_time = time.time() - start_time
    
    logging.info(f"\n--- Extraction Summary ---")
    logging.info(f"Total structure mappings processed: {len(structures_df)}")
    logging.info(f"Unique UniProt IDs: {stats['unique_uniprot_ids']}")
    logging.info(f"Total fragments: {stats['total_combinations']}")
    logging.info(f"Found in archive: {stats['found_in_archive']}")
    logging.info(f"AlphaFold PDB files extracted/found: {total_structure_files}")
    logging.info(f"Available structures: {available_structures}")
    logging.info(f"Updated mapping saved to: {OUTPUT_MAPPING_CSV}")
    logging.info(f"Total processing time: {elapsed_time:.2f} seconds")
    
    # Verify the output
    missing_structures = alphafold_df[alphafold_df['AlphaFold_Available'] == False]
    if len(missing_structures) > 0:
        logging.warning(f"Warning: {len(missing_structures)} entries missing AlphaFold structures")
        logging.warning("These fragments may not be available in the local AlphaFold database")
    else:
        logging.info("All entries have corresponding AlphaFold structures!")

if __name__ == "__main__":
    main()
