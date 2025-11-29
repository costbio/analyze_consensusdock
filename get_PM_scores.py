#!/usr/bin/env python3
"""
PocketMatch Similarity Calculator for Cavity Files (Parallelized)

This script runs PocketMatch on cavity PDB files extracted from AlphaFold structures
to compute pairwise pocket similarity scores (Pmin and Pmax).

Supports parallel processing by dividing pockets into chunks and using Step3-PM_typeB
for comparisons between chunks and the full set.

Input: Directory containing cavity PDB files in subdirectories
Output: CSV file with pairwise pocket similarity scores
"""

import pandas as pd
import subprocess
import os
import glob
import shutil
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
import tempfile
from functools import partial
from tqdm import tqdm
import concurrent.futures
import json
import logging
import faulthandler
import signal
from biopandas.pdb import PandasPdb

faulthandler.register(signal.SIGUSR1)
# Module-level logger (configured in main entry)
logger = logging.getLogger("pocketmatch")

def _setup_logger(log_path: str):
    """Configure file logger once (idempotent)."""
    logger.setLevel(logging.INFO)
    # Avoid adding multiple handlers on repeated calls
    if not any(isinstance(h, logging.FileHandler) and getattr(h, 'baseFilename', None) == os.path.abspath(log_path) for h in logger.handlers):
        # File handler - captures ALL levels (DEBUG and above)
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.DEBUG)  # File gets everything
        fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        fh.setFormatter(fmt)
        logger.addHandler(fh)
        
        # Console handler (stderr) - only INFO and above by default
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)  # Console only shows INFO+ unless verbose mode
        ch.setFormatter(fmt)
        logger.addHandler(ch)


def collect_cavity_files(cavities_dir):
    """
    Collect all cavity PDB files from subdirectories.
    
    Args:
        cavities_dir: Path to extracted_cavities directory
        
    Returns:
        List of paths to cavity PDB files
    """
    logger.info(f"Collecting cavity PDB files from: {cavities_dir}")
    
    # Find all cavity PDB files (not vacant files)
    cavity_files = list(Path(cavities_dir).glob("*/AF-*_cavity_*.pdb"))
    
    logger.info(f"Found {len(cavity_files):,} cavity PDB files")
    
    return cavity_files


def create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, workspace_id):
    """
    Create an isolated PocketMatch workspace by copying only essential files.
    Since workspaces are created sequentially, we can safely copy files without
    hitting "too many open files" errors.
    
    Args:
        pm_base_dir: Base PocketMatch installation directory
        work_dir: Working directory for this process
        workspace_id: Unique identifier for this workspace
        
    Returns:
        Path to the isolated PocketMatch directory
    """
    isolated_pm_dir = os.path.join(work_dir, f"pm_workspace_{workspace_id}")
    
    # Create workspace directory
    os.makedirs(isolated_pm_dir, exist_ok=True)
    
    # Copy executable file (Step3-PM_typeB is needed for comparisons)
    step3_typeB = os.path.join(pm_base_dir, "Step3-PM_typeB")
    if os.path.exists(step3_typeB):
        dest_step3 = os.path.join(isolated_pm_dir, "Step3-PM_typeB")
        if not os.path.exists(dest_step3):
            shutil.copy2(step3_typeB, dest_step3)
    
    # Create cabbage-file_maker directory structure
    cabbage_dir = os.path.join(isolated_pm_dir, "cabbage-file_maker")
    os.makedirs(cabbage_dir, exist_ok=True)
    
    # Copy essential cabbage-file_maker executables and scripts
    base_cabbage_dir = os.path.join(pm_base_dir, "cabbage-file_maker")
    for item in ["Step0-cabbage.sh", "Step0-cabbage_core", "Step0-END-FILE"]:
        src = os.path.join(base_cabbage_dir, item)
        if os.path.exists(src):
            dest = os.path.join(cabbage_dir, item)
            if not os.path.exists(dest):
                shutil.copy2(src, dest)
    
    # Create empty Sample_pockets directory (will be populated later)
    sample_pockets_dir = os.path.join(cabbage_dir, "Sample_pockets")
    os.makedirs(sample_pockets_dir, exist_ok=True)
    
    # Verify the structure was created successfully
    if not os.path.isdir(cabbage_dir):
        raise RuntimeError(f"Failed to create cabbage-file_maker directory: {cabbage_dir}")
    if not os.path.isdir(sample_pockets_dir):
        raise RuntimeError(f"Failed to create Sample_pockets directory: {sample_pockets_dir}")
    
    return isolated_pm_dir


def extract_complete_residues_from_cavity(args):
    """
    Extract complete residues from AlphaFold structure for a single cavity file.
    Worker function for parallel processing.
    
    Uses BioPandas to properly parse PDB files and ensure complete residue extraction.
    
    Args:
        args: Tuple of (idx, pocket_file, af_base_dir, dest_file)
        
    Returns:
        Tuple of (simple_name, original_name, success_status)
    """
    idx, pocket_file, af_base_dir, dest_file = args
    
    simple_name = f"p{idx+1}.pdb"
    original_name = os.path.basename(pocket_file)
    
    try:
        # Parse cavity file using BioPandas to get residue information
        cavity_ppdb = PandasPdb().read_pdb(str(pocket_file))
        cavity_atoms = cavity_ppdb.df['ATOM']
        
        if len(cavity_atoms) == 0:
            return (simple_name, original_name, False)
        
        # Get unique residues from cavity file (chain + residue number)
        cavity_residues = set(
            zip(cavity_atoms['chain_id'], cavity_atoms['residue_number'])
        )
        
        # Parse filename to extract protein ID and fragment
        # Format: AF-PROTEINID-FX-model_vY_cavity_N.pdb
        pocket_basename = os.path.basename(pocket_file)
        if not (pocket_basename.startswith('AF-') and '_cavity_' in pocket_basename):
            # Not an AlphaFold cavity file - skip it
            return (simple_name, original_name, False)
        
        parts = pocket_basename.split('-')
        if len(parts) < 4:
            return (simple_name, original_name, False)
        
        protein_id = parts[1]  # Extract PROTEINID
        fragment_id = parts[2]  # Extract FX (F1, F2, F3, etc.)
        
        # Find AlphaFold PDB file for the specific fragment
        af_dir = os.path.join(af_base_dir, protein_id, fragment_id)
        af_pdb_pattern = f"AF-{protein_id}-{fragment_id}-model_v*.pdb"
        af_pdb_files = list(Path(af_dir).glob(af_pdb_pattern))
        
        if not af_pdb_files:
            # AlphaFold structure not found - skip this cavity
            return (simple_name, original_name, False)
        
        af_pdb_path = af_pdb_files[0]
        
        # Load AlphaFold structure using BioPandas
        af_ppdb = PandasPdb().read_pdb(str(af_pdb_path))
        af_atoms = af_ppdb.df['ATOM']
        
        if len(af_atoms) == 0:
            logger.debug(f"No ATOM records in AlphaFold structure for {original_name}")
            return (simple_name, original_name, False)
        
        # CRITICAL VALIDATION: Check that ALL cavity residues have CA atoms in AlphaFold structure
        # PocketMatch Step0-cabbage_core allocates memory based on CA count
        # Missing CA atoms cause heap corruption and malloc crashes
        
        # Get residues with CA atoms in AlphaFold structure
        ca_atoms = af_atoms[af_atoms['atom_name'] == 'CA']
        af_residues_with_ca = set(
            zip(ca_atoms['chain_id'], ca_atoms['residue_number'])
        )
        
        # Find which cavity residues have CA atoms in AlphaFold
        valid_residues = cavity_residues & af_residues_with_ca
        
        if len(valid_residues) == 0:
            # No valid residues with CA atoms - this shouldn't happen with AlphaFold
            # but could happen if cavity file is corrupted
            return (simple_name, original_name, False)
        
        if len(valid_residues) < len(cavity_residues):
            # Some cavity residues don't have CA atoms in AlphaFold structure
            # This indicates data corruption - skip this file
            missing = cavity_residues - valid_residues
            logger.warning(f"{original_name} has residues without CA atoms: {missing}")
            return (simple_name, original_name, False)
        
        # Extract all atoms for valid residues from AlphaFold structure
        mask = af_atoms.apply(
            lambda row: (row['chain_id'], row['residue_number']) in valid_residues,
            axis=1
        )
        extracted_atoms = af_atoms[mask].copy()
        
        if len(extracted_atoms) == 0:
            return (simple_name, original_name, False)
        
        # Write extracted atoms to destination file
        # Create a new PandasPdb object and write properly formatted PDB
        output_ppdb = PandasPdb()
        output_ppdb.df['ATOM'] = extracted_atoms
        
        # Write to file with proper PDB formatting (80 character lines)
        output_ppdb.to_pdb(path=dest_file, records=['ATOM'], gz=False)
        
        return (simple_name, original_name, True)
        
    except Exception as e:
        logger.warning(f"Failed to process {original_name}: {e}")
        import traceback
        traceback.print_exc()
        return (simple_name, original_name, False)


def generate_cabbage_file(pm_dir, pocket_files, output_name="outfile.cabbage", af_base_dir=None, show_progress=False, n_threads=None):
    """
    Generate a cabbage file for a set of pocket files.
    
    Args:
        pm_dir: PocketMatch directory
        pocket_files: List of pocket PDB file paths
        output_name: Name for the cabbage output file
        af_base_dir: Base directory containing AlphaFold structures (required for extracting complete residues)
        show_progress: Whether to show progress bar
        n_threads: Number of threads for parallel residue extraction (default: CPU count)
        
    Returns:
        Path to generated cabbage file, dict mapping simplified names to original names
    """
    # Verify workspace structure exists
    cabbage_dir = os.path.join(pm_dir, "cabbage-file_maker")
    if not os.path.exists(cabbage_dir):
        raise RuntimeError(f"cabbage-file_maker directory does not exist: {cabbage_dir}. Workspace may not be properly initialized.")
    
    sample_pockets_dir = os.path.join(cabbage_dir, "Sample_pockets")
    
    # Clean and prepare Sample_pockets directory
    if os.path.exists(sample_pockets_dir):
        for file in Path(sample_pockets_dir).glob("*.pdb"):
            file.unlink()
    else:
        os.makedirs(sample_pockets_dir)
    
    # Extract COMPLETE RESIDUES from AlphaFold structures in parallel
    # Use simple numeric names to avoid issues with long names or special characters
    name_mapping = {}
    
    # Prepare arguments for parallel processing
    if n_threads is None:
        n_threads = cpu_count()
    
    # Use the provided n_threads value directly
    # When called from parallel block processing, this will be 1 (sequential)
    max_workers = n_threads
    
    args_list = [
        (idx, pocket_file, af_base_dir, os.path.join(sample_pockets_dir, f"p{idx+1}.pdb"))
        for idx, pocket_file in enumerate(pocket_files)
    ]
    
    # Process residue extraction in parallel with timeout
    # Use ProcessPoolExecutor for true parallelism (bypasses Python GIL)
    results = []
    logger.debug(f"Extracting complete residues from {len(pocket_files)} cavity files using {max_workers} workers...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks and track futures
        future_to_idx = {executor.submit(extract_complete_residues_from_cavity, args): i for i, args in enumerate(args_list)}
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_idx):
            try:
                result = future.result(timeout=30)  # 30s per file
                results.append(result)
                idx = future_to_idx[future]
                logger.debug(f"Extracted residues from file {idx+1}/{len(pocket_files)}")
            except concurrent.futures.TimeoutError:
                idx = future_to_idx[future]
                logger.error(f"Timeout processing file index {idx}")
                results.append((f"p{idx+1}.pdb", "TIMEOUT", False))
            except Exception as e:
                idx = future_to_idx[future]
                logger.error(f"Error processing file index {idx}: {e}")
                results.append((f"p{idx+1}.pdb", "ERROR", False))
    
    # Build name mapping from results
    for simple_name, original_name, _ in results:
        name_mapping[simple_name] = original_name
    
    # Generate cabbage file using Step0-cabbage.sh
    cabbage_dir = os.path.join(pm_dir, "cabbage-file_maker")
    
    # Verify outfile.cabbage does not exist (should have been cleaned in create_isolated_pocketmatch_workspace)
    outfile_path = os.path.join(cabbage_dir, "outfile.cabbage")
    if os.path.exists(outfile_path):
        raise RuntimeError(f"outfile.cabbage still exists at {outfile_path} - cleanup failed")
    
    # Run Step0-cabbage.sh with Sample_pockets directory as argument
    sample_pockets_dir = os.path.join(cabbage_dir, "Sample_pockets")
    
    if show_progress:
        logger.info("Generating cabbage file...")
    
    logger.debug(f"Running Step0-cabbage.sh with Sample_pockets: {sample_pockets_dir}")
    cmd = ["./Step0-cabbage.sh", sample_pockets_dir]
    # Use cwd parameter instead of os.chdir() to avoid threading issues
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cabbage_dir)
    
    if result.returncode != 0:
        raise RuntimeError(f"Cabbage generation failed: {result.stderr}")
    
    if not os.path.exists(outfile_path):
        raise FileNotFoundError(f"Cabbage file not created: {outfile_path}")
    
    # Rename if needed
    if output_name != "outfile.cabbage":
        new_path = os.path.join(cabbage_dir, output_name)
        shutil.move(outfile_path, new_path)
        cabbage_file = new_path
    else:
        cabbage_file = outfile_path
    
    return cabbage_file, name_mapping


def run_pocketmatch_typeB(pm_dir, cabbage_chunk, cabbage_full):
    """
    Run PocketMatch typeB comparison (chunk vs full set).
    
    Args:
        pm_dir: PocketMatch directory
        cabbage_chunk: Path to cabbage file for chunk
        cabbage_full: Path to cabbage file for full set
        
    Returns:
        Path to output file
    """
    # Copy cabbage files to PM directory if not already there
    chunk_local = os.path.join(pm_dir, os.path.basename(cabbage_chunk))
    full_local = os.path.join(pm_dir, os.path.basename(cabbage_full))
    
    if cabbage_chunk != chunk_local:
        shutil.copy2(cabbage_chunk, chunk_local)
    if cabbage_full != full_local:
        shutil.copy2(cabbage_full, full_local)
    
    # Run PocketMatch typeB
    cmd = ["./Step3-PM_typeB", os.path.basename(chunk_local), os.path.basename(full_local)]
    # Use cwd parameter instead of os.chdir() to avoid threading issues
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=pm_dir)
    
    if result.returncode != 0:
        raise RuntimeError(f"PocketMatch typeB failed with code {result.returncode}. stderr: {result.stderr}, stdout: {result.stdout}")
    
    output_file = os.path.join(pm_dir, "PocketMatch_score.txt")
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"PocketMatch output not created: {output_file}")
    
    return output_file


def process_chunk(chunk_id, chunk_files, all_files, pm_base_dir, work_dir, af_base_dir=None, debug=False):
    """
    Process a single chunk of pocket files in parallel.
    
    Args:
        chunk_id: Unique identifier for this chunk
        chunk_files: List of pocket files in this chunk
        all_files: List of all pocket files
        pm_base_dir: Base PocketMatch directory
        work_dir: Working directory for temporary files
        af_base_dir: Base directory for AlphaFold structures
        debug: Whether to keep temporary files
        
    Returns:
        DataFrame with similarity results for this chunk
    """
    try:
        logger.info(f"Worker {chunk_id}: Processing {len(chunk_files)} pockets...")
        
        # Create isolated workspace
        isolated_pm = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, chunk_id)
        
        # Generate cabbage file for chunk
        cabbage_chunk, chunk_mapping = generate_cabbage_file(isolated_pm, chunk_files, f"chunk_{chunk_id}.cabbage", af_base_dir)
        
        # Generate cabbage file for full set (reuse if exists)
        cabbage_full_path = os.path.join(work_dir, "full_set.cabbage")
        mapping_path = os.path.join(work_dir, "name_mapping.json")
        
        if not os.path.exists(cabbage_full_path):
            # Only generate once (first worker to need it)
            cabbage_full, full_mapping = generate_cabbage_file(isolated_pm, all_files, "full_set.cabbage", af_base_dir)
            shutil.copy2(cabbage_full, cabbage_full_path)
            # Save mapping for reuse
            import json
            with open(mapping_path, 'w') as f:
                json.dump(full_mapping, f)
        else:
            # Copy existing full set cabbage and load mapping
            cabbage_full_dest = os.path.join(isolated_pm, "cabbage-file_maker", "full_set.cabbage")
            shutil.copy2(cabbage_full_path, cabbage_full_dest)
            cabbage_full = cabbage_full_dest
            import json
            with open(mapping_path, 'r') as f:
                full_mapping = json.load(f)
        
        # Run PocketMatch typeB
        output_file = run_pocketmatch_typeB(isolated_pm, cabbage_chunk, cabbage_full)
        
        # Parse results with name mapping
        df_chunk = parse_pocketmatch_output(output_file, name_mapping2=full_mapping)
        
        logger.info(f"Worker {chunk_id}: Completed - {len(df_chunk)} comparisons")
        
        # Cleanup isolated workspace
        if not debug:
            shutil.rmtree(isolated_pm, ignore_errors=True)
        
        return df_chunk
        
    except Exception as e:
        logger.error(f"Worker {chunk_id}: Error - {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()


def build_cabbage_block(block_id, pocket_files, pm_base_dir, work_dir, af_base_dir, max_internal_threads=1, show_progress=False, verbose=False, isolated_pm=None):
    """
    Build a cabbage file for a block (<=100 pockets) in an isolated workspace.
    Returns (cabbage_path, name_mapping, tmp_workspace_path).
    The workspace is cleaned by the caller after consuming outputs if needed.
    """
    logger.info(f"Block {block_id}: Building cabbage with {len(pocket_files)} pockets")
    if verbose:
        logger.debug(f"Block {block_id}: Starting build with {len(pocket_files)} pockets")
    # Use pre-created workspace if provided, otherwise create new one
    if isolated_pm is None:
        isolated_pm = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, f"blk_{block_id}")
    try:
        cabbage_path, name_mapping = generate_cabbage_file(
            isolated_pm,
            pocket_files,
            output_name=f"block_{block_id}.cabbage",
            af_base_dir=af_base_dir,
            show_progress=show_progress,
            n_threads=max_internal_threads,
        )
        # Move outputs to a shared blocks dir for later comparisons
        blocks_dir = os.path.join(work_dir, "cabbage_blocks")
        os.makedirs(blocks_dir, exist_ok=True)
        final_cabbage = os.path.join(blocks_dir, f"block_{block_id}.cabbage")
        shutil.copy2(cabbage_path, final_cabbage)
        # Persist mapping
        mapping_path = os.path.join(blocks_dir, f"block_{block_id}.mapping.json")
        with open(mapping_path, 'w') as f:
            json.dump(name_mapping, f)
        logger.info(f"Block {block_id}: Finished -> {final_cabbage}")
        if verbose:
            logger.debug(f"Block {block_id}: Completed successfully")
        return final_cabbage, name_mapping, isolated_pm
    except Exception as e:
        # Don't delete workspace here - will be cleaned up by caller after all blocks finish
        # Deleting during parallel execution causes race conditions where other blocks
        # try to access files that have been deleted
        logger.error(f"Block {block_id}: Error during build - {e}")
        raise


def run_block_comparison(pair_id, cabbage_a, cabbage_b, mapping_a, mapping_b, pm_base_dir, work_dir, debug=False, verbose=False, isolated_pm=None):
    """
    Run Step3-PM_typeB for two cabbage blocks and parse with separate mappings.
    Returns DataFrame of results for this pair.
    """
    try:
        logger.info(f"Pair {pair_id}: Starting comparison {os.path.basename(cabbage_a)} vs {os.path.basename(cabbage_b)}")
        if verbose:
            logger.debug(f"Pair {pair_id}: Starting comparison")
        # Use pre-created workspace if provided, otherwise create new one
        if isolated_pm is None:
            isolated_pm = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, f"cmp_{pair_id}")
        
        # Copy cabbage files for this comparison into workspace
        local_a = os.path.join(isolated_pm, os.path.basename(cabbage_a))
        local_b = os.path.join(isolated_pm, os.path.basename(cabbage_b))
        shutil.copy2(cabbage_a, local_a)
        shutil.copy2(cabbage_b, local_b)

        # Run comparison
        output_file = run_pocketmatch_typeB(isolated_pm, local_a, local_b)

        # Parse using separate mappings: first set uses mapping_a, second uses mapping_b
        df = parse_pocketmatch_output(output_file, name_mapping1=mapping_a, name_mapping2=mapping_b)
        
        # Optionally save raw results to a separate folder (only if debug mode)
        if debug:
            results_dir = os.path.join(work_dir, "raw_results")
            os.makedirs(results_dir, exist_ok=True)
            result_file = os.path.join(results_dir, f"result_{pair_id}.txt")
            shutil.copy2(output_file, result_file)
        
        # Report completion to log
        completion_msg = f"Pair {pair_id}: Completed with {len(df)} comparisons"
        logger.info(completion_msg)

        # IMMEDIATELY clean up workspace (don't wait for finally block)
        # This prevents disk from getting clogged with thousands of workspaces
        if isolated_pm and not debug:
            try:
                shutil.rmtree(isolated_pm, ignore_errors=False)
                logger.debug(f"Pair {pair_id}: Cleaned up workspace")
            except Exception as cleanup_error:
                logger.warning(f"Pair {pair_id}: Failed to cleanup workspace: {cleanup_error}")
        
        return df
        
    except Exception as e:
        logger.error(f"Pair {pair_id}: Error - {e}")
        import traceback
        import sys
        # Log full traceback to file
        tb_str = ''.join(traceback.format_exception(*sys.exc_info()))
        logger.error(f"Pair {pair_id}: Traceback:\n{tb_str}")
        return pd.DataFrame()
    finally:
        # Final cleanup attempt (in case immediate cleanup above failed)
        if isolated_pm and not debug:
            if os.path.exists(isolated_pm):
                try:
                    shutil.rmtree(isolated_pm, ignore_errors=True)
                except:
                    pass  # Silent fail in finally block


def parse_pocketmatch_output(output_file, name_mapping=None, name_mapping1=None, name_mapping2=None):
    """
    Parse PocketMatch output file into a DataFrame.
    
    Args:
        output_file: Path to PocketMatch_score.txt
        name_mapping: Dict mapping simplified names (p1.pdb) to original names
        
    Returns:
        DataFrame with columns: Pocket1, Pocket2, Pmin, Pmax
    """
    valid_lines = []
    total_lines = 0
    
    with open(output_file, 'r') as f:
        for line in f:
            total_lines += 1
            try:
                fields = line.strip().split()
                if len(fields) >= 4:
                    valid_lines.append(fields[:4])
            except:
                pass
    
    if not valid_lines:
        return pd.DataFrame(columns=['Pocket1', 'Pocket2', 'Pmin', 'Pmax'])
    
    # Create DataFrame
    df = pd.DataFrame(valid_lines, columns=['Pocket1', 'Pocket2', 'Pmin', 'Pmax'])
    
    # Clean up data
    df = df[df['Pmax'] != 'NULL____']
    df = df[df['Pmin'] != 'NULL____']
    
    # Convert to numeric
    df['Pmin'] = pd.to_numeric(df['Pmin'], errors='coerce')
    df['Pmax'] = pd.to_numeric(df['Pmax'], errors='coerce')
    
    # Remove rows with invalid scores
    df = df.dropna(subset=['Pmin', 'Pmax'])
    
    # Remove trailing underscores from pocket names
    df['Pocket1'] = df['Pocket1'].str.rstrip('_')
    df['Pocket2'] = df['Pocket2'].str.rstrip('_')
    
    # Map simplified names back to original names if mapping provided
    # Backward-compat: if name_mapping is provided, use it for both sides
    if name_mapping and not (name_mapping1 or name_mapping2):
        name_mapping1 = name_mapping
        name_mapping2 = name_mapping

    def make_mapper(mapping):
        if not mapping:
            return lambda x: x
        def map_name(x):
            if x in mapping:
                return mapping[x]
            if (x + '.pdb') in mapping:
                return mapping[x + '.pdb']
            if x.endswith('.pdb'):
                base = x[:-4]
                if base in mapping:
                    return mapping[base]
            return x
        return map_name

    df['Pocket1'] = df['Pocket1'].map(make_mapper(name_mapping1))
    df['Pocket2'] = df['Pocket2'].map(make_mapper(name_mapping2))
    
    return df


def get_pm_similarities_parallel(cavities_dir, output_dir, pm_base_dir, af_base_dir=None, n_threads=None, debug=False, verbose=False, test_mode=False):
    """
    Main function to compute PocketMatch similarities with parallelization.
    
    Args:
        cavities_dir: Directory containing extracted cavity files
        output_dir: Directory to save output CSV
        pm_base_dir: Base PocketMatch installation directory
        af_base_dir: Base directory for AlphaFold structures (optional, for complete residue extraction)
        n_threads: Number of parallel threads (default: CPU count)
        debug: Whether to keep temporary files
        verbose: Whether to print detailed progress for each block/comparison
        test_mode: Whether to run in test mode (only first 3 blocks)
        
    Returns:
        DataFrame with pocket similarity scores
    """
    logger.info("="*80)
    logger.info("POCKETMATCH SIMILARITY ANALYSIS (PARALLEL)")
    logger.info("="*80)
    
    # Validate paths
    if not os.path.exists(cavities_dir):
        raise FileNotFoundError(f"Cavities directory not found: {cavities_dir}")
    
    if not os.path.exists(pm_base_dir):
        raise FileNotFoundError(f"PocketMatch directory not found: {pm_base_dir}")
    
    # Create output directory FIRST
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup file logger after output directory exists
    log_path = os.path.join(output_dir, 'pocketmatch_progress.log')
    _setup_logger(log_path)
    if verbose:
        logger.setLevel(logging.DEBUG)
        # Also set console handler to DEBUG level in verbose mode
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                handler.setLevel(logging.DEBUG)
    logger.info("Starting PocketMatch analysis")
    logger.info(f"Logging to: {log_path}")
    
    # Determine number of threads (use conservative 20 instead of all CPUs to avoid resource exhaustion)
    if n_threads is None:
        n_threads = min(20, cpu_count())
    
    logger.info(f"Configuration:")
    logger.info(f"Parallel threads: {n_threads}")
    
    # Collect cavity files
    cavity_files = collect_cavity_files(cavities_dir)
    
    if len(cavity_files) == 0:
        raise ValueError("No cavity PDB files found!")
    
    # Divide files into fixed-size blocks for cabbage generation (max 100 per block)
    block_size = 100
    blocks = [cavity_files[i:i + block_size] for i in range(0, len(cavity_files), block_size)]
    
    # Test mode: limit to first 3 blocks for quick testing
    if test_mode:
        blocks = blocks[:3]
        logger.warning("TEST MODE: Processing only first 3 blocks (up to 300 pockets)")
    
    logger.info(f"Preparing {len(blocks)} cabbage blocks (size <= {block_size})")
    
    # Create temporary working directory
    work_dir = os.path.join(output_dir, "temp_workspaces")
    os.makedirs(work_dir, exist_ok=True)
    
    try:
        # PRE-CREATE all workspaces sequentially to avoid "too many open files" errors
        logger.info(f"Pre-creating {len(blocks)} workspaces sequentially...")
        block_workspaces = {}
        for i in tqdm(range(len(blocks)), desc="   Creating workspaces"):
            ws = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, f"blk_{i}")
            block_workspaces[i] = ws
        
        # Generate all cabbage blocks in parallel using pre-created workspaces
        logger.info("Generating cabbage blocks in parallel...")
        built_blocks = []  # list of (cabbage_path, mapping, tmp_workspace)
        
        # Each block processes files with 1 worker internally (sequential per block)
        # But we run n_threads blocks in parallel
        logger.info(f"Using {n_threads} parallel blocks (each processes files with 1 worker)")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as ex:
            futures = {
                ex.submit(build_cabbage_block, i, blk, pm_base_dir, work_dir, af_base_dir, 1, False, verbose, block_workspaces[i]): i
                for i, blk in enumerate(blocks)
            }
            completed = tqdm(concurrent.futures.as_completed(futures, timeout=7200), total=len(futures), desc="   Blocks")
            for fut in completed:
                try:
                    result = fut.result(timeout=1200)  # 20min per block
                    built_blocks.append(result)
                except concurrent.futures.TimeoutError:
                    block_id = futures[fut]
                    logger.error(f"Block {block_id}: Timeout after 20 minutes")
                except Exception as e:
                    block_id = futures[fut]
                    logger.error(f"Block {block_id}: Error - {e}")
        logger.info(f"Built {len(built_blocks)} cabbage blocks successfully")

        # Clean up temp workspaces used for building (keep only block files)
        for _, _, ws in built_blocks:
            shutil.rmtree(ws, ignore_errors=True)

        # Prepare all pairwise comparisons between blocks (including self)
        total_pairs = (len(built_blocks) * (len(built_blocks) + 1)) // 2
        logger.info(f"Submitting {total_pairs} block comparisons (upper triangle)...")
        pair_jobs = []
        for i, (cab_i, map_i, _) in enumerate(built_blocks):
            for j in range(i, len(built_blocks)):
                cab_j, map_j, _ = built_blocks[j]
                pair_jobs.append((f"{i}_{j}", cab_i, cab_j, map_i, map_j))

        # PRE-CREATE comparison workspaces sequentially
        logger.info(f"Pre-creating {len(pair_jobs)} comparison workspaces sequentially...")
        cmp_workspaces = {}
        for pid, _, _, _, _ in tqdm(pair_jobs, desc="   Creating cmp workspaces"):
            ws = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, f"cmp_{pid}")
            cmp_workspaces[pid] = ws

        # Run comparisons with bounded concurrency and timeouts
        # Memory-efficient processing: save results in batches to avoid RAM overflow
        BATCH_SIZE = 100000  # Save every 100k rows to prevent memory issues
        all_pair_results = []
        batch_counter = 0
        saved_chunks = []
        chunk_dir = os.path.join(output_dir, 'result_chunks')
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as ex:
            futs = {
                ex.submit(run_block_comparison, pid, cab_i, cab_j, map_i, map_j, pm_base_dir, work_dir, debug, verbose, cmp_workspaces[pid]): pid
                for (pid, cab_i, cab_j, map_i, map_j) in pair_jobs
            }
            completed = tqdm(concurrent.futures.as_completed(futs, timeout=14400), total=len(futs), desc="   Comparisons")
            for fut in completed:
                try:
                    result = fut.result(timeout=300)  # 5min per comparison
                    all_pair_results.append(result)
                except concurrent.futures.TimeoutError:
                    pair_id = futs[fut]
                    logger.error(f"Pair {pair_id}: Timeout after 5 minutes")
                    all_pair_results.append(pd.DataFrame())
                except Exception as e:
                    pair_id = futs[fut]
                    logger.error(f"Pair {pair_id}: Error - {e}")
                    all_pair_results.append(pd.DataFrame())
                
                # Check if we should save a batch to disk
                valid_results = [df for df in all_pair_results if df is not None and len(df) > 0]
                if valid_results:
                    total_rows = sum(len(df) for df in valid_results)
                    if total_rows >= BATCH_SIZE:
                        logger.info(f"Saving batch {batch_counter} with ~{total_rows:,} rows to free memory...")
                        os.makedirs(chunk_dir, exist_ok=True)
                        chunk_file = os.path.join(chunk_dir, f'chunk_{batch_counter:04d}.csv')
                        batch_df = pd.concat(valid_results, ignore_index=True)
                        batch_df.to_csv(chunk_file, index=False)
                        saved_chunks.append(chunk_file)
                        batch_counter += 1
                        # Clear memory
                        all_pair_results = []
                        del batch_df
                        logger.info(f"Batch saved, memory cleared. Continuing...")
        
        logger.info(f"Completed {len(futs)} block comparisons")

        # Combine final results (remaining in-memory results + saved chunks)
        logger.info(f"Combining results from {batch_counter} saved batches + remaining in-memory results...")
        final_dfs = []
        
        # Add any remaining in-memory results
        valid_results = [df for df in all_pair_results if df is not None and len(df) > 0]
        if valid_results:
            logger.info(f"Processing {len(valid_results)} remaining in-memory results...")
            final_dfs.extend(valid_results)
        
        # Load saved chunks
        if saved_chunks:
            logger.info(f"Loading {len(saved_chunks)} saved chunk files...")
            for chunk_file in saved_chunks:
                chunk_df = pd.read_csv(chunk_file)
                final_dfs.append(chunk_df)
                logger.debug(f"Loaded {len(chunk_df):,} rows from {os.path.basename(chunk_file)}")
        
        if not final_dfs:
            raise ValueError("No valid results from any comparison. Check error messages above for details.")
        
        # Process and save results in chunks to avoid memory overflow
        logger.info(f"Processing {len(final_dfs)} dataframes in batches...")
        
        # Create final output directory for result files
        results_dir = os.path.join(output_dir, 'PocketMatch_similarities')
        os.makedirs(results_dir, exist_ok=True)
        
        total_rows = 0
        total_duplicates_removed = 0
        output_files = []
        
        # Process in batches to manage memory
        FINAL_BATCH_SIZE = 100000
        current_batch = []
        current_batch_rows = 0
        file_counter = 0
        
        for df in final_dfs:
            current_batch.append(df)
            current_batch_rows += len(df)
            
            if current_batch_rows >= FINAL_BATCH_SIZE:
                # Process and save this batch
                batch_df = pd.concat(current_batch, ignore_index=True)
                initial_len = len(batch_df)
                batch_df = batch_df.drop_duplicates(subset=['Pocket1', 'Pocket2'])
                duplicates_removed = initial_len - len(batch_df)
                total_duplicates_removed += duplicates_removed
                
                # Sort by Pmax
                batch_df = batch_df.sort_values(by='Pmax', ascending=False).reset_index(drop=True)
                
                # Save to file
                output_file = os.path.join(results_dir, f'similarities_part_{file_counter:04d}.csv')
                batch_df.to_csv(output_file, index=False)
                output_files.append(output_file)
                
                logger.info(f"Saved part {file_counter}: {len(batch_df):,} rows to {os.path.basename(output_file)}")
                total_rows += len(batch_df)
                file_counter += 1
                
                # Clear memory
                current_batch = []
                current_batch_rows = 0
                del batch_df
        
        # Process remaining data
        if current_batch:
            batch_df = pd.concat(current_batch, ignore_index=True)
            initial_len = len(batch_df)
            batch_df = batch_df.drop_duplicates(subset=['Pocket1', 'Pocket2'])
            duplicates_removed = initial_len - len(batch_df)
            total_duplicates_removed += duplicates_removed
            
            batch_df = batch_df.sort_values(by='Pmax', ascending=False).reset_index(drop=True)
            
            output_file = os.path.join(results_dir, f'similarities_part_{file_counter:04d}.csv')
            batch_df.to_csv(output_file, index=False)
            output_files.append(output_file)
            
            logger.info(f"Saved part {file_counter}: {len(batch_df):,} rows to {os.path.basename(output_file)}")
            total_rows += len(batch_df)
            del batch_df
        
        # Clean up temporary chunk files
        if saved_chunks and os.path.exists(chunk_dir):
            logger.info(f"Cleaning up temporary chunk files...")
            shutil.rmtree(chunk_dir, ignore_errors=True)
        
        # Log summary statistics
        logger.info(f"="*80)
        logger.info(f"Summary Statistics:")
        logger.info(f"Total output files: {len(output_files)}")
        logger.info(f"Results directory: {results_dir}")
        logger.info(f"Total comparisons: {total_rows:,}")
        if total_duplicates_removed > 0:
            logger.info(f"Removed {total_duplicates_removed:,} duplicate comparisons")
        
        # Calculate overall statistics by reading first file as sample
        if output_files:
            sample_df = pd.read_csv(output_files[0])
            logger.info(f"Sample statistics from first file:")
            logger.info(f"  Pmax range: {sample_df['Pmax'].min():.3f} - {sample_df['Pmax'].max():.3f}")
            logger.info(f"  Pmin range: {sample_df['Pmin'].min():.3f} - {sample_df['Pmin'].max():.3f}")
            logger.info(f"  Mean Pmax: {sample_df['Pmax'].mean():.3f}")
            logger.info(f"  Mean Pmin: {sample_df['Pmin'].mean():.3f}")
            del sample_df
        
        return output_files
        
    finally:
        # Cleanup temporary workspaces
        if not debug:
            logger.info(f"Cleaning up temporary workspaces...")
            shutil.rmtree(work_dir, ignore_errors=True)
            logger.info(f"Cleanup complete")
        else:
            logger.info(f"Debug mode: Keeping temporary workspaces at: {work_dir}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run PocketMatch similarity analysis on cavity files (parallelized).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python get_PM_scores.py \\
    -cavities /media/onur/Elements/cavity_space_consensus_docking/2025_06_29_batch_dock/extracted_cavities \\
    -output /media/onur/Elements/cavity_space_consensus_docking/2025_06_29_batch_dock/pocketmatch_results \\
    -pm_dir /home/onur/software/PocketMatch_v2.1/PocketMatch_v2.1 \\
    -af_structures /media/onur/Elements/cavity_space_consensus_docking/2025_06_29_batch_dock/alphafold_structures \\
    -threads 8

Note: Uses PocketMatch Step3-PM_typeB for parallelization.
      Each thread processes a chunk of pockets against the full set.
      Temporary isolated workspaces are created to avoid conflicts.
        """
    )
    
    parser.add_argument(
        '-cavities',
        type=str,
        required=True,
        help='Path to extracted_cavities directory containing cavity PDB files'
    )
    
    parser.add_argument(
        '-output',
        type=str,
        required=True,
        help='Path to output directory for results'
    )
    
    parser.add_argument(
        '-pm_dir',
        type=str,
        required=True,
        help='Path to PocketMatch base directory'
    )
    
    parser.add_argument(
        '-af_structures',
        type=str,
        required=False,
        default=None,
        help='Path to AlphaFold structures directory (for extracting complete residues)'
    )
    
    parser.add_argument(
        '-threads',
        type=int,
        default=None,
        help=f'Number of parallel threads (default: {cpu_count()} CPUs)'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Debug mode: keep temporary workspaces for inspection'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose mode: print detailed progress for each block and comparison'
    )
    
    parser.add_argument(
        '--test-mode',
        action='store_true',
        help='Test mode: process only first 3 blocks (300 pockets max) for quick testing'
    )
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    cavities_dir = os.path.abspath(args.cavities)
    output_dir = os.path.abspath(args.output)
    pm_base_dir = os.path.abspath(args.pm_dir)
    af_base_dir = os.path.abspath(args.af_structures) if args.af_structures else None
    
    # Run analysis
    try:
        df_sims = get_pm_similarities_parallel(
            cavities_dir,
            output_dir,
            pm_base_dir,
            af_base_dir=af_base_dir,
            n_threads=args.threads,
            debug=args.debug,
            verbose=args.verbose,
            test_mode=args.test_mode
        )
        logger.info("="*80)
        logger.info("ANALYSIS COMPLETE!")
        logger.info("="*80)
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
