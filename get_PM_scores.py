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


def collect_cavity_files(cavities_dir):
    """
    Collect all cavity PDB files from subdirectories.
    
    Args:
        cavities_dir: Path to extracted_cavities directory
        
    Returns:
        List of paths to cavity PDB files
    """
    print(f"\n🔍 Collecting cavity PDB files from: {cavities_dir}")
    
    # Find all cavity PDB files (not vacant files)
    cavity_files = list(Path(cavities_dir).glob("*/AF-*_cavity_*.pdb"))
    
    print(f"✅ Found {len(cavity_files):,} cavity PDB files")
    
    return cavity_files


def create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, workspace_id):
    """
    Create an isolated PocketMatch workspace for parallel processing.
    
    Args:
        pm_base_dir: Base PocketMatch directory to copy from
        work_dir: Working directory for this process
        workspace_id: Unique identifier for this workspace
        
    Returns:
        Path to the isolated PocketMatch directory
    """
    isolated_pm_dir = os.path.join(work_dir, f"pm_workspace_{workspace_id}")
    
    # Copy entire PocketMatch directory to isolated workspace
    shutil.copytree(pm_base_dir, isolated_pm_dir, dirs_exist_ok=True)
    
    return isolated_pm_dir


def generate_cabbage_file(pm_dir, pocket_files, output_name="outfile.cabbage", af_base_dir=None, show_progress=False):
    """
    Generate a cabbage file for a set of pocket files.
    
    Args:
        pm_dir: PocketMatch directory
        pocket_files: List of pocket PDB file paths
        output_name: Name for the cabbage output file
        af_base_dir: Base directory containing AlphaFold structures (required for extracting complete residues)
        show_progress: Whether to show progress bar
        
    Returns:
        Path to generated cabbage file, dict mapping simplified names to original names
    """
    sample_pockets_dir = os.path.join(pm_dir, "cabbage-file_maker", "Sample_pockets")
    
    # Clean and prepare Sample_pockets directory
    if os.path.exists(sample_pockets_dir):
        for file in Path(sample_pockets_dir).glob("*.pdb"):
            file.unlink()
    else:
        os.makedirs(sample_pockets_dir)
    
    # Copy and clean pocket files to Sample_pockets with simplified names
    # Extract COMPLETE RESIDUES from AlphaFold structures (not just cavity atoms)
    # Use simple numeric names to avoid issues with long names or special characters
    name_mapping = {}
    
    iterator = enumerate(pocket_files)
    if show_progress:
        iterator = tqdm(iterator, total=len(pocket_files), desc="   Extracting complete residues", leave=False)
    
    for idx, pocket_file in iterator:
        # Create simple name: p1.pdb, p2.pdb, etc.
        simple_name = f"p{idx+1}.pdb"
        dest_file = os.path.join(sample_pockets_dir, simple_name)
        name_mapping[simple_name] = os.path.basename(pocket_file)
        
        # Parse cavity file to get residues
        residues_to_keep = set()
        with open(pocket_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain = line[21:22]
                    res_num = line[22:26].strip()
                    residues_to_keep.add((chain, res_num))
        
        # Find corresponding AlphaFold structure
        # Extract protein ID from filename: AF-PROTEINID-F1-model_v1_cavity_N.pdb
        pocket_basename = os.path.basename(pocket_file)
        if pocket_basename.startswith('AF-') and '_cavity_' in pocket_basename:
            protein_id = pocket_basename.split('-')[1]  # Extract PROTEINID
            
            # Find AlphaFold PDB file
            af_dir = os.path.join(af_base_dir, protein_id, "F1")
            af_pdb_files = list(Path(af_dir).glob(f"AF-{protein_id}-F1-model_v*.pdb"))
            
            if af_pdb_files:
                af_pdb = af_pdb_files[0]
                # Extract complete residues from AlphaFold structure
                with open(af_pdb, 'r') as infile, open(dest_file, 'w') as outfile:
                    for line in infile:
                        if line.startswith('ATOM'):
                            chain = line[21:22]
                            res_num = line[22:26].strip()
                            if (chain, res_num) in residues_to_keep:
                                # Keep original line, ensure 80 chars
                                line_stripped = line.rstrip('\n\r')
                                line_80 = line_stripped.ljust(80)[:80]
                                outfile.write(line_80 + '\n')
            else:
                # Fallback: use cavity file as-is if AlphaFold structure not found
                with open(pocket_file, 'r') as infile, open(dest_file, 'w') as outfile:
                    for line in infile:
                        if line.startswith('ATOM'):
                            line_stripped = line.rstrip('\n\r')
                            line_80 = line_stripped.ljust(80)[:80]
                            outfile.write(line_80 + '\n')
        else:
            # Not an AlphaFold cavity file, use as-is
            with open(pocket_file, 'r') as infile, open(dest_file, 'w') as outfile:
                for line in infile:
                    if line.startswith('ATOM'):
                        line_stripped = line.rstrip('\n\r')
                        line_80 = line_stripped.ljust(80)[:80]
                        outfile.write(line_80 + '\n')
    
    # Generate cabbage file
    # Process each PDB file individually to avoid memory corruption in Step0-cabbage_core
    original_dir = os.getcwd()
    try:
        cabbage_dir = os.path.join(pm_dir, "cabbage-file_maker")
        os.chdir(cabbage_dir)
        
        # Remove old outfile if exists
        outfile_path = os.path.join(cabbage_dir, "outfile.cabbage")
        if os.path.exists(outfile_path):
            os.remove(outfile_path)
        
        # Process each PDB file individually and append to outfile.cabbage
        sample_pockets_dir = os.path.join(cabbage_dir, "Sample_pockets")
        pdb_files = sorted(Path(sample_pockets_dir).glob("*.pdb"))
        
        iterator = pdb_files
        if show_progress:
            iterator = tqdm(pdb_files, desc="   Generating cabbage file", leave=False)
        
        for pdb_file in iterator:
            # Run Step0-cabbage_core for each file and append
            cmd = ["./Step0-cabbage_core", str(pdb_file)]
            result = subprocess.run(cmd, capture_output=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"Cabbage generation failed for {pdb_file.name}: {result.stderr}")
            
            # Append the output to outfile.cabbage
            with open(outfile_path, 'ab') as outfile:
                outfile.write(result.stdout)
        
        # Add END-FILE marker
        cmd = ["./Step0-END-FILE"]
        result = subprocess.run(cmd, capture_output=True)
        with open(outfile_path, 'ab') as outfile:
            outfile.write(result.stdout)
        
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
        
    finally:
        os.chdir(original_dir)


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
    original_dir = os.getcwd()
    try:
        os.chdir(pm_dir)
        
        # Copy cabbage files to PM directory if not already there
        chunk_local = os.path.join(pm_dir, os.path.basename(cabbage_chunk))
        full_local = os.path.join(pm_dir, os.path.basename(cabbage_full))
        
        if cabbage_chunk != chunk_local:
            shutil.copy2(cabbage_chunk, chunk_local)
        if cabbage_full != full_local:
            shutil.copy2(cabbage_full, full_local)
        
        # Run PocketMatch typeB
        cmd = ["./Step3-PM_typeB", os.path.basename(chunk_local), os.path.basename(full_local)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"PocketMatch typeB failed with code {result.returncode}. stderr: {result.stderr}, stdout: {result.stdout}")
        
        output_file = os.path.join(pm_dir, "PocketMatch_score.txt")
        if not os.path.exists(output_file):
            raise FileNotFoundError(f"PocketMatch output not created: {output_file}")
        
        return output_file
        
    finally:
        os.chdir(original_dir)


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
        print(f"   Worker {chunk_id}: Processing {len(chunk_files)} pockets...")
        
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
        df_chunk = parse_pocketmatch_output(output_file, full_mapping)
        
        print(f"   Worker {chunk_id}: ✅ Completed - {len(df_chunk)} comparisons")
        
        # Cleanup isolated workspace
        if not debug:
            shutil.rmtree(isolated_pm, ignore_errors=True)
        
        return df_chunk
        
    except Exception as e:
        print(f"   Worker {chunk_id}: ❌ Error - {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()


def parse_pocketmatch_output(output_file, name_mapping=None):
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
    if name_mapping:
        # Handle both "p1" and "p1.pdb" formats from PocketMatch output
        def map_name(x):
            # Try exact match first
            if x in name_mapping:
                return name_mapping[x]
            # Try adding .pdb extension
            if x + '.pdb' in name_mapping:
                return name_mapping[x + '.pdb']
            # Try removing .pdb extension
            if x.endswith('.pdb'):
                base = x[:-4]
                if base in name_mapping:
                    return name_mapping[base]
            # Return original if no mapping found
            return x
        
        df['Pocket1'] = df['Pocket1'].map(map_name)
        df['Pocket2'] = df['Pocket2'].map(map_name)
    
    return df


def get_pm_similarities_parallel(cavities_dir, output_dir, pm_base_dir, af_base_dir=None, n_threads=None, debug=False):
    """
    Main function to compute PocketMatch similarities with parallelization.
    
    Args:
        cavities_dir: Directory containing extracted cavity files
        output_dir: Directory to save output CSV
        pm_base_dir: Base PocketMatch installation directory
        af_base_dir: Base directory for AlphaFold structures (optional, for complete residue extraction)
        n_threads: Number of parallel threads (default: CPU count)
        debug: Whether to keep temporary files
        
    Returns:
        DataFrame with pocket similarity scores
    """
    print("="*80)
    print("🔬 POCKETMATCH SIMILARITY ANALYSIS (PARALLEL)")
    print("="*80)
    
    # Validate paths
    if not os.path.exists(cavities_dir):
        raise FileNotFoundError(f"Cavities directory not found: {cavities_dir}")
    
    if not os.path.exists(pm_base_dir):
        raise FileNotFoundError(f"PocketMatch directory not found: {pm_base_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine number of threads
    if n_threads is None:
        n_threads = cpu_count()
    
    print(f"\n⚙️  Configuration:")
    print(f"   Parallel threads: {n_threads}")
    
    # Collect cavity files
    cavity_files = collect_cavity_files(cavities_dir)
    
    if len(cavity_files) == 0:
        raise ValueError("No cavity PDB files found!")
    
    # Divide files into chunks
    chunk_size = max(1, len(cavity_files) // n_threads)
    chunks = [cavity_files[i:i + chunk_size] for i in range(0, len(cavity_files), chunk_size)]
    
    print(f"\n📦 Dividing into {len(chunks)} chunks:")
    print(f"   Chunk size: ~{chunk_size} pockets")
    print(f"   Total comparisons: ~{len(cavity_files) * len(cavity_files):,}")
    
    # Create temporary working directory
    work_dir = os.path.join(output_dir, "temp_workspaces")
    os.makedirs(work_dir, exist_ok=True)
    
    try:
        # Generate full set cabbage file once (shared by all workers)
        print(f"\n🔧 Preparing shared full-set cabbage file...")
        temp_pm = create_isolated_pocketmatch_workspace(pm_base_dir, work_dir, "temp")
        cabbage_full, name_mapping = generate_cabbage_file(temp_pm, cavity_files, "full_set.cabbage", af_base_dir, show_progress=True)
        full_cabbage_path = os.path.join(work_dir, "full_set.cabbage")
        mapping_path = os.path.join(work_dir, "name_mapping.json")
        shutil.copy2(cabbage_full, full_cabbage_path)
        # Save name mapping
        import json
        with open(mapping_path, 'w') as f:
            json.dump(name_mapping, f)
        shutil.rmtree(temp_pm, ignore_errors=True)
        print(f"   ✅ Shared cabbage file ready")
        
        # Process chunks in parallel
        print(f"\n🚀 Processing chunks in parallel...")
        
        process_func = partial(
            process_chunk,
            all_files=cavity_files,
            pm_base_dir=pm_base_dir,
            work_dir=work_dir,
            af_base_dir=af_base_dir,
            debug=debug
        )
        
        with Pool(processes=n_threads) as pool:
            chunk_args = [(i, chunk) for i, chunk in enumerate(chunks)]
            results = pool.starmap(process_func, chunk_args)
        
        # Combine results
        print(f"\n📊 Combining results from {len(results)} chunks...")
        valid_results = [df for df in results if len(df) > 0]
        
        if not valid_results:
            raise ValueError("No valid results from any worker. Check error messages above for details.")
        
        df_sims = pd.concat(valid_results, ignore_index=True)
        
        # Remove duplicates (typeB may create some)
        initial_len = len(df_sims)
        df_sims = df_sims.drop_duplicates(subset=['Pocket1', 'Pocket2'])
        
        if len(df_sims) < initial_len:
            print(f"   Removed {initial_len - len(df_sims):,} duplicate comparisons")
        
        # Sort by Pmax
        df_sims = df_sims.sort_values(by='Pmax', ascending=False).reset_index(drop=True)
        
        # Save results
        output_csv = os.path.join(output_dir, 'PocketMatch_similarities.csv')
        df_sims.to_csv(output_csv, index=False)
        print(f"\n💾 Results saved to: {output_csv}")
        
        # Print summary statistics
        print(f"\n📈 Summary Statistics:")
        print(f"   Total comparisons: {len(df_sims):,}")
        print(f"   Pmax range: {df_sims['Pmax'].min():.3f} - {df_sims['Pmax'].max():.3f}")
        print(f"   Pmin range: {df_sims['Pmin'].min():.3f} - {df_sims['Pmin'].max():.3f}")
        print(f"   Mean Pmax: {df_sims['Pmax'].mean():.3f}")
        print(f"   Mean Pmin: {df_sims['Pmin'].mean():.3f}")
        
        # Show top 10 most similar pockets
        print(f"\n🏆 Top 10 Most Similar Pocket Pairs (by Pmax):")
        print(df_sims.head(10).to_string(index=False))
        
        return df_sims
        
    finally:
        # Cleanup temporary workspaces
        if not debug:
            print(f"\n🧹 Cleaning up temporary workspaces...")
            shutil.rmtree(work_dir, ignore_errors=True)
            print(f"   ✅ Cleanup complete")
        else:
            print(f"\n🔍 Debug mode: Keeping temporary workspaces at: {work_dir}")


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
            debug=args.debug
        )
        print("\n" + "="*80)
        print("✅ ANALYSIS COMPLETE!")
        print("="*80)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
