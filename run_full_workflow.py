#!/usr/bin/env python3
"""
Complete Workflow Runner
This script runs the entire consensus docking workflow in the correct order.

Usage:
    python run_full_workflow.py [--skip-uniprot] [--skip-extract] [--skip-download] [--skip-convert]

Example:
    # Run full workflow
    python run_full_workflow.py
    
    # Skip UniProt mapping generation (if already done)
    python run_full_workflow.py --skip-uniprot
    
    # Skip cavity extraction (if already done)  
    python run_full_workflow.py --skip-extract
    
    # Skip AlphaFold download
    python run_full_workflow.py --skip-download
    
    # Skip PDBQT conversion
    python run_full_workflow.py --skip-convert
"""

import os
import sys
import argparse
import subprocess
import logging
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_script(script_name, description):
    """Run a Python script and handle errors."""
    logging.info(f"Starting: {description}")
    start_time = time.time()
    
    try:
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True, check=True)
        elapsed = time.time() - start_time
        logging.info(f"Completed: {description} ({elapsed:.1f}s)")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed: {description}")
        logging.error(f"Error output: {e.stderr}")
        return False
    except FileNotFoundError:
        logging.error(f"Script not found: {script_name}")
        return False

def check_file_exists(filename, description):
    """Check if a required file exists."""
    if os.path.exists(filename):
        logging.info(f"Found: {description} ({filename})")
        return True
    else:
        logging.warning(f"Missing: {description} ({filename})")
        return False

def main():
    """Main workflow runner."""
    parser = argparse.ArgumentParser(description='Run the complete consensus docking workflow')
    parser.add_argument('--skip-uniprot', action='store_true', 
                       help='Skip UniProt mapping generation')
    parser.add_argument('--skip-extract', action='store_true',
                       help='Skip cavity extraction')
    parser.add_argument('--skip-download', action='store_true',
                       help='Skip AlphaFold model download')
    parser.add_argument('--skip-convert', action='store_true',
                       help='Skip PDB to PDBQT conversion')
    parser.add_argument('--test-mode', action='store_true',
                       help='Run batch_dock.py in test mode')
    
    args = parser.parse_args()
    
    logging.info("=== Consensus Docking Workflow Runner ===")
    start_total = time.time()
    
    # Step 1: UniProt Mapping Generation
    if not args.skip_uniprot:
        if not run_script('prepare_uniprot_mapping.py', 'UniProt mapping generation'):
            logging.error("Failed at step 1. Stopping workflow.")
            return 1
    else:
        check_file_exists('uniprot_gene_mapping.csv', 'UniProt gene mapping')
    
    # Step 2: Cavity Extraction  
    if not args.skip_extract:
        if not run_script('extract_cavities.py', 'Cavity extraction'):
            logging.error("Failed at step 2. Stopping workflow.")
            return 1
    else:
        check_file_exists('cavity_mapping.csv', 'Cavity mapping')
    
    # Step 3: AlphaFold Model Download (Recommended)
    if not args.skip_download:
        if not run_script('download_cavityspace_structures.py', 'CavitySpace structure download'):
            logging.warning("AlphaFold download failed. Continuing with extracted cavity data.")
            # Don't stop workflow - this is optional but recommended
    else:
        check_file_exists('alphafold_mapping.csv', 'AlphaFold mapping')
    
    # Step 4: PDB to PDBQT Conversion (Optional but recommended)
    if not args.skip_convert:
        if not run_script('convert_pdb_to_pdbqt.py', 'PDB to PDBQT conversion'):
            logging.warning("PDBQT conversion failed. Docking will convert files on-the-fly.")
            # Don't stop workflow - this is optional
    else:
        check_file_exists('pdbqt_mapping.csv', 'PDBQT mapping')
    
    # Step 5: Batch Docking
    logging.info("Starting: Batch docking execution")
    
    # Check prerequisites
    required_files = [
        ('uniprot_gene_mapping.csv', 'UniProt gene mapping'),
        ('cavity_mapping.csv', 'Cavity mapping')
    ]
    
    missing_files = []
    for filename, description in required_files:
        if not check_file_exists(filename, description):
            missing_files.append(filename)
    
    if missing_files:
        logging.error(f"Missing required files: {missing_files}")
        logging.error("Cannot proceed with batch docking.")
        return 1
    
    # Run batch docking
    batch_script = 'batch_dock.py'
    if args.test_mode:
        logging.info("Running batch docking in test mode")
        # You would need to modify batch_dock.py to accept command line arguments
        # For now, this just runs it normally
    
    if not run_script(batch_script, 'Batch docking execution'):
        logging.error("Failed at final step: batch docking.")
        return 1
    
    # Success!
    total_elapsed = time.time() - start_total
    logging.info("=== Workflow Completed Successfully! ===")
    logging.info(f"Total time: {total_elapsed:.1f}s ({total_elapsed/60:.1f}m)")
    
    # Summary of outputs
    logging.info("\nGenerated files:")
    output_files = [
        'uniprot_gene_mapping.csv',
        'cavity_mapping.csv',
        'alphafold_mapping.csv',
        'pdbqt_mapping.csv',
        'extracted_cavities/',
        'alphafold_models/',
        'converted_pdbqt/',
        'consensus_docking_results/'
    ]
    
    for filename in output_files:
        if os.path.exists(filename):
            logging.info(f"  ✓ {filename}")
        else:
            logging.info(f"  - {filename} (not created)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
