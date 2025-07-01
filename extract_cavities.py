#!/usr/bin/env python3
"""
Cavity Extraction Script
This script extracts receptor and cavity PDB files from AlphaFold cavity tarballs.

This is a separate preprocessing step that can be run independently of the main
docking workflow, allowing for better organization and the ability to cache
extracted files.

Usage:
    python extract_cavities.py

Input:
    - Cavity tarball folder containing .tar.gz files
    
Output:
    - extracted_cavities/ folder with extracted PDB files
    - cavity_mapping.csv: CSV file mapping UniProt IDs to extracted cavity pairs
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

def extract_pdb_paths_from_tarballs(cavity_tarball_folder, extract_base_dir="extracted_cavities"):
    """
    Iterates through tar.gz files, extracts receptor and pocket PDBs,
    and returns a dictionary mapping UniProt ID to lists of (receptor_path, pocket_path, fragment_id).
    
    Args:
        cavity_tarball_folder (str): Path to folder containing .tar.gz files
        extract_base_dir (str): Base directory for extraction
        
    Returns:
        dict: Dictionary mapping UniProt IDs to lists of (receptor_pdb, pocket_pdb, fragment_id) tuples
    """
    logging.info(f"Extracting PDB paths from tarballs in '{cavity_tarball_folder}'...")
    pdb_info = {} # { "Q00577": [(vacant_pdb1, cavity_pdb1, fragment_id1), (vacant_pdb2, cavity_pdb2, fragment_id2)], ... }

    os.makedirs(extract_base_dir, exist_ok=True)

    tar_files = [f for f in os.listdir(cavity_tarball_folder) if f.endswith(".tar.gz")]
    
    if not tar_files:
        logging.warning(f"No .tar.gz files found in '{cavity_tarball_folder}'.")
        return pdb_info

    logging.info(f"Found {len(tar_files)} tarball files to process.")

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

        # Dictionary to store multiple cavity pairs: {cavity_num: (receptor_pdb, cavity_pdb)}
        cavity_pairs = {}

        try:
            with tarfile.open(full_tar_path, "r:gz") as tar:
                # Find the relevant files and extract them
                for member in tar.getmembers():
                    # Match vacant files with any number: _vacant_1.pdb, _vacant_2.pdb, etc.
                    vacant_match = re.search(r'_vacant_(\d+)\.pdb$', member.name)
                    if vacant_match:
                        cavity_num = int(vacant_match.group(1))
                        tar.extract(member, path=extract_folder)
                        receptor_pdb = os.path.join(extract_folder, member.name)
                        if cavity_num not in cavity_pairs:
                            cavity_pairs[cavity_num] = [None, None]
                        cavity_pairs[cavity_num][0] = receptor_pdb
                    
                    # Match cavity files with any number: _cavity_1.pdb, _cavity_2.pdb, etc.
                    cavity_match = re.search(r'_cavity_(\d+)\.pdb$', member.name)
                    if cavity_match:
                        cavity_num = int(cavity_match.group(1))
                        tar.extract(member, path=extract_folder)
                        pocket_pdb = os.path.join(extract_folder, member.name)
                        if cavity_num not in cavity_pairs:
                            cavity_pairs[cavity_num] = [None, None]
                        cavity_pairs[cavity_num][1] = pocket_pdb
            
            # Add all complete cavity pairs to pdb_info
            valid_pairs_found = 0
            for cavity_num, (receptor_pdb, pocket_pdb) in cavity_pairs.items():
                if receptor_pdb and pocket_pdb:
                    if uniprot_id not in pdb_info:
                        pdb_info[uniprot_id] = []
                    pdb_info[uniprot_id].append((receptor_pdb, pocket_pdb, fragment_id))
                    valid_pairs_found += 1
                else:
                    missing = "receptor" if not receptor_pdb else "pocket"
                    tqdm.write(f"Warning: Missing {missing} PDB for cavity {cavity_num} in '{tar_gz_filename}'.")
            
            if valid_pairs_found == 0:
                tqdm.write(f"Warning: No complete receptor-pocket pairs found in '{tar_gz_filename}'. Skipping.")
            else:
                tqdm.write(f"Found {valid_pairs_found} valid cavity pair(s) for fragment F{fragment_id} in '{tar_gz_filename}'.")

        except tarfile.ReadError as e:
            tqdm.write(f"Error reading tar.gz file '{tar_gz_filename}': {e}. Skipping.")
        except Exception as e:
            tqdm.write(f"An unexpected error occurred with '{tar_gz_filename}': {e}. Skipping.")

    logging.info(f"Finished extracting PDB paths. Found info for {len(pdb_info)} UniProt IDs.")
    return pdb_info

def save_cavity_mapping_to_csv(pdb_info_dict, output_csv):
    """
    Save the cavity mapping information to a CSV file for later use.
    
    Args:
        pdb_info_dict (dict): Dictionary mapping UniProt IDs to cavity pairs
        output_csv (str): Path to output CSV file
    """
    logging.info(f"Saving cavity mapping to '{output_csv}'...")
    
    rows = []
    for uniprot_id, cavity_pairs in pdb_info_dict.items():
        for i, (receptor_pdb, pocket_pdb, fragment_id) in enumerate(cavity_pairs):
            rows.append({
                'UniProt_ID': uniprot_id,
                'Fragment_ID': fragment_id,
                'Cavity_Index': i + 1,
                'Receptor_PDB': receptor_pdb,
                'Pocket_PDB': pocket_pdb
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    
    logging.info(f"Saved {len(rows)} cavity pairs for {len(pdb_info_dict)} UniProt IDs to '{output_csv}'.")
    return df

def load_cavity_mapping_from_csv(input_csv):
    """
    Load cavity mapping from a previously generated CSV file.
    
    Args:
        input_csv (str): Path to input CSV file
        
    Returns:
        dict: Dictionary mapping UniProt IDs to lists of (receptor_pdb, pocket_pdb, fragment_id) tuples
    """
    if not os.path.exists(input_csv):
        logging.error(f"Cavity mapping file not found: {input_csv}")
        return {}
    
    try:
        df = pd.read_csv(input_csv)
        required_columns = ['UniProt_ID', 'Fragment_ID', 'Receptor_PDB', 'Pocket_PDB']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            logging.error(f"Required columns missing from cavity mapping file: {missing_columns}")
            return {}
        
        pdb_info = {}
        for _, row in df.iterrows():
            uniprot_id = row['UniProt_ID']
            fragment_id = row['Fragment_ID']
            receptor_pdb = row['Receptor_PDB']
            pocket_pdb = row['Pocket_PDB']
            
            if uniprot_id not in pdb_info:
                pdb_info[uniprot_id] = []
            pdb_info[uniprot_id].append((receptor_pdb, pocket_pdb, fragment_id))
        
        logging.info(f"Loaded cavity mapping for {len(pdb_info)} UniProt IDs from '{input_csv}'.")
        return pdb_info
        
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
        response = input("Cavity mapping already exists. Do you want to re-extract? (y/N): ").strip().lower()
        if response not in ['y', 'yes']:
            logging.info("Using existing cavity mapping. Exiting.")
            # Load and display summary
            pdb_info = load_cavity_mapping_from_csv(OUTPUT_CSV)
            total_pairs = sum(len(pairs) for pairs in pdb_info.values())
            logging.info(f"Summary: {len(pdb_info)} UniProt IDs with {total_pairs} total cavity pairs.")
            return
    
    try:
        # Extract cavities
        pdb_data = extract_pdb_paths_from_tarballs(CAVITY_TARBALL_FOLDER, EXTRACT_BASE_DIR)
        
        if not pdb_data:
            logging.warning("No cavity data extracted. Exiting.")
            sys.exit(0)
        
        # Save mapping to CSV
        df = save_cavity_mapping_to_csv(pdb_data, OUTPUT_CSV)
        
        # Summary statistics
        total_pairs = len(df)
        unique_uniprots = len(pdb_data)
        avg_pairs_per_protein = total_pairs / unique_uniprots if unique_uniprots > 0 else 0
        
        elapsed_time = time.time() - start_time
        
        logging.info(f"\n--- Extraction Summary ---")
        logging.info(f"Total UniProt IDs processed: {unique_uniprots}")
        logging.info(f"Total cavity pairs extracted: {total_pairs}")
        logging.info(f"Average cavity pairs per protein: {avg_pairs_per_protein:.2f}")
        logging.info(f"Extraction completed in {elapsed_time:.2f} seconds")
        logging.info(f"Results saved to: {OUTPUT_CSV}")
        logging.info(f"Extracted files location: {EXTRACT_BASE_DIR}")
        
    except Exception as e:
        logging.critical(f"Error during cavity extraction: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
