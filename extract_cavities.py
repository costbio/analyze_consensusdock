#!/usr/bin/env python3
"""
Cavity Extraction Script
This script extracts cavity/pocket PDB files from AlphaFold cavity tarballs.

This is a separate preprocessing step that can be run independently of the main
docking workflow, allowing for better organization and the ability to cache
extracted cavity files. The receptor PDB information comes from AlphaFold 
structures, not from the cavity tarballs.

Usage:
    python extract_cavities.py

Input:
    - Cavity tarball folder containing .tar.gz files
    
Output:
    - extracted_cavities/ folder with extracted cavity PDB files
    - cavity_mapping.csv: CSV file mapping UniProt IDs to extracted cavities
"""

import os
import sys
import pandas as pd
import tarfile
from tqdm import tqdm
from pathlib import Path
import re
import time
import logging

# --- Configuration ---
# --- IMPORTANT: ADJUST THESE PATHS AS PER YOUR SYSTEM ---
CAVITY_TARBALL_FOLDER = "/opt/data/cavity_space/af_strong_cavities_data_all"
EXTRACT_BASE_DIR = "extracted_cavities"  # Will be created in current working directory
OUTPUT_CSV = "cavity_mapping.csv"
LOG_FILE = "cavity_extraction.log"

# Set up logging
logging.basicConfig(filename=LOG_FILE, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logging.getLogger().addHandler(console_handler)

def extract_cavity_paths_from_tarballs(cavity_tarball_folder, extract_base_dir="extracted_cavities"):
    """
    Iterates through tar.gz files, extracts cavity/pocket PDBs only,
    and returns a dictionary mapping UniProt ID to lists of (pocket_path, fragment_id, cavity_index).
    
    Args:
        cavity_tarball_folder (str): Path to folder containing .tar.gz files
        extract_base_dir (str): Base directory for extraction
        
    Returns:
        dict: Dictionary mapping UniProt IDs to lists of (pocket_pdb, fragment_id, cavity_index) tuples
    """
    logging.info(f"Extracting cavity paths from tarballs in '{cavity_tarball_folder}'...")
    cavity_info = {} # { "Q00577": [(cavity_pdb1, fragment_id1, cavity_index1), (cavity_pdb2, fragment_id2, cavity_index2)], ... }

    os.makedirs(extract_base_dir, exist_ok=True)

    tar_files = [f for f in os.listdir(cavity_tarball_folder) if f.endswith(".tar.gz")]
    
    if not tar_files:
        logging.warning(f"No .tar.gz files found in '{cavity_tarball_folder}'.")
        return cavity_info

    logging.info(f"Found {len(tar_files)} tarball files to process.")
    
    # Statistics tracking
    extracted_count = 0
    skipped_count = 0

    for tar_gz_filename in tqdm(tar_files, desc="Extracting Tarballs", unit="file"):
        full_tar_path = os.path.join(cavity_tarball_folder, tar_gz_filename)
        
        # Extract UniProt ID and fragment info from filename
        # Pattern: AF-Q96RW7-F23-model_v1_cavity_result_1.tar.gz
        match = re.search(r'AF-([A-Z0-9]+)-F(\d+)-model_v1', tar_gz_filename)
        if match:
            uniprot_id = match.group(1)
            fragment_id = int(match.group(2))
        else:
            tqdm.write(f"Warning: Could not extract UniProt ID and fragment from '{tar_gz_filename}'. Skipping.")
            continue

        extract_folder = os.path.join(extract_base_dir, Path(tar_gz_filename).stem.replace(".tar", ""))
        os.makedirs(extract_folder, exist_ok=True)

        # Dictionary to store cavity files: {cavity_num: cavity_pdb}
        cavity_files = {}

        # Check if extraction folder already exists and has files
        existing_cavity_files = []
        if os.path.exists(extract_folder):
            for file_path in Path(extract_folder).glob("*_cavity_*.pdb"):
                cavity_match = re.search(r'_cavity_(\d+)\.pdb$', file_path.name)
                if cavity_match:
                    cavity_num = int(cavity_match.group(1))
                    existing_cavity_files.append((cavity_num, str(file_path)))
        
        if existing_cavity_files:
            # Use existing files instead of extracting
            tqdm.write(f"Found {len(existing_cavity_files)} existing cavity files for '{tar_gz_filename}'. Skipping extraction.")
            for cavity_num, file_path in existing_cavity_files:
                cavity_files[cavity_num] = file_path
            skipped_count += 1
        else:
            # Extract cavity files
            try:
                with tarfile.open(full_tar_path, "r:gz") as tar:
                    # Find and extract only cavity files
                    for member in tar.getmembers():
                        # Match cavity files with any number: _cavity_1.pdb, _cavity_2.pdb, etc.
                        cavity_match = re.search(r'_cavity_(\d+)\.pdb$', member.name)
                        if cavity_match:
                            cavity_num = int(cavity_match.group(1))
                            tar.extract(member, path=extract_folder)
                            pocket_pdb = os.path.join(extract_folder, member.name)
                            cavity_files[cavity_num] = pocket_pdb
                extracted_count += 1
            
            except tarfile.ReadError as e:
                tqdm.write(f"Error reading tar.gz file '{tar_gz_filename}': {e}. Skipping.")
                continue
            except Exception as e:
                tqdm.write(f"An unexpected error occurred with '{tar_gz_filename}': {e}. Skipping.")
                continue
        # Add all cavity files to cavity_info
        cavities_found = 0
        for cavity_num, pocket_pdb in cavity_files.items():
            if uniprot_id not in cavity_info:
                cavity_info[uniprot_id] = []
            cavity_info[uniprot_id].append((pocket_pdb, fragment_id, cavity_num))
            cavities_found += 1
        
        if cavities_found == 0:
            tqdm.write(f"Warning: No cavity files found in '{tar_gz_filename}'. Skipping.")
        else:
            extraction_status = "existing" if existing_cavity_files else "extracted"
            tqdm.write(f"Found {cavities_found} cavity file(s) ({extraction_status}) for fragment F{fragment_id} in '{tar_gz_filename}'.")

    logging.info(f"Finished extracting cavity paths. Found info for {len(cavity_info)} UniProt IDs.")
    logging.info(f"Extraction statistics: {extracted_count} tarballs extracted, {skipped_count} tarballs skipped (already present)")
    return cavity_info

def save_cavity_mapping_to_csv(cavity_info_dict, output_csv):
    """
    Save the cavity mapping information to a CSV file for later use.
    
    Args:
        cavity_info_dict (dict): Dictionary mapping UniProt IDs to cavity files
        output_csv (str): Path to output CSV file
    """
    logging.info(f"Saving cavity mapping to '{output_csv}'...")
    
    rows = []
    for uniprot_id, cavity_list in cavity_info_dict.items():
        for pocket_pdb, fragment_id, cavity_index in cavity_list:
            rows.append({
                'UniProt_ID': uniprot_id,
                'Fragment_ID': fragment_id,
                'Cavity_Index': cavity_index,
                'Pocket_PDB': pocket_pdb
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    
    logging.info(f"Saved {len(rows)} cavity files for {len(cavity_info_dict)} UniProt IDs to '{output_csv}'.")
    return df

def load_cavity_mapping_from_csv(input_csv):
    """
    Load cavity mapping from a previously generated CSV file.
    
    Args:
        input_csv (str): Path to input CSV file
        
    Returns:
        dict: Dictionary mapping UniProt IDs to lists of (pocket_pdb, fragment_id, cavity_index) tuples
    """
    if not os.path.exists(input_csv):
        logging.error(f"Cavity mapping file not found: {input_csv}")
        return {}
    
    try:
        df = pd.read_csv(input_csv)
        required_columns = ['UniProt_ID', 'Fragment_ID', 'Pocket_PDB', 'Cavity_Index']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            logging.error(f"Required columns missing from cavity mapping file: {missing_columns}")
            return {}
        
        cavity_info = {}
        for _, row in df.iterrows():
            uniprot_id = row['UniProt_ID']
            fragment_id = row['Fragment_ID']
            pocket_pdb = row['Pocket_PDB']
            cavity_index = row['Cavity_Index']
            
            if uniprot_id not in cavity_info:
                cavity_info[uniprot_id] = []
            cavity_info[uniprot_id].append((pocket_pdb, fragment_id, cavity_index))
        
        logging.info(f"Loaded cavity mapping for {len(cavity_info)} UniProt IDs from '{input_csv}'.")
        return cavity_info
        
    except Exception as e:
        logging.error(f"Error loading cavity mapping file: {e}")
        return {}

def main():
    """Main function to extract cavities and save mapping."""
    start_time = time.time()
    
    logging.info("Starting cavity extraction workflow.")
    logging.info(f"Cavity Tarball Folder: {CAVITY_TARBALL_FOLDER}")
    logging.info(f"Extract Base Directory: {EXTRACT_BASE_DIR}")
    logging.info(f"Output CSV: {OUTPUT_CSV}")
    
    # Validate input folder
    if not os.path.exists(CAVITY_TARBALL_FOLDER):
        logging.critical(f"Error: Cavity tarball folder not found at '{CAVITY_TARBALL_FOLDER}'. Please correct the path.")
        sys.exit(1)
    
    # Check if extraction was already done
    if os.path.exists(OUTPUT_CSV):
        logging.info(f"Found existing cavity mapping file: {OUTPUT_CSV}")
        response = input("Cavity mapping already exists. Do you want to regenerate it? (y/N): ").strip().lower()
        if response not in ['y', 'yes']:
            logging.info("Using existing cavity mapping. Exiting.")
            # Load and display summary
            cavity_info = load_cavity_mapping_from_csv(OUTPUT_CSV)
            total_cavities = sum(len(cavities) for cavities in cavity_info.values())
            logging.info(f"Summary: {len(cavity_info)} UniProt IDs with {total_cavities} total cavity files.")
            return
        else:
            logging.info("Will regenerate cavity mapping. Existing extracted files will be reused to save time.")
    else:
        logging.info("No existing cavity mapping found. Will extract cavities and create mapping.")
    
    try:
        # Extract cavities
        cavity_data = extract_cavity_paths_from_tarballs(CAVITY_TARBALL_FOLDER, EXTRACT_BASE_DIR)
        
        if not cavity_data:
            logging.warning("No cavity data extracted. Exiting.")
            sys.exit(0)
        
        # Save mapping to CSV
        df = save_cavity_mapping_to_csv(cavity_data, OUTPUT_CSV)
        
        # Summary statistics
        total_cavities = len(df)
        unique_uniprots = len(cavity_data)
        avg_cavities_per_protein = total_cavities / unique_uniprots if unique_uniprots > 0 else 0
        
        elapsed_time = time.time() - start_time
        
        logging.info(f"\n--- Extraction Summary ---")
        logging.info(f"Total UniProt IDs processed: {unique_uniprots}")
        logging.info(f"Total cavity files found: {total_cavities}")
        logging.info(f"Average cavities per protein: {avg_cavities_per_protein:.2f}")
        logging.info(f"Tarballs processed: {len([f for f in os.listdir(CAVITY_TARBALL_FOLDER) if f.endswith('.tar.gz')])}")
        logging.info(f"Time saved by reusing existing files: Significant!")
        logging.info(f"Extraction completed in {elapsed_time:.2f} seconds")
        logging.info(f"Results saved to: {OUTPUT_CSV}")
        logging.info(f"Extracted files location: {EXTRACT_BASE_DIR}")
        
    except Exception as e:
        logging.critical(f"Error during cavity extraction: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
