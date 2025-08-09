#!/usr/bin/env python3
"""
High-Performance Consensus Docking Data Parser and Preparer
===========================================================

This script optimizes the parsing and preparation of consensus docking results
with multiprocessing, progress bars, and performance optimizations.

Features:
- Multiprocessing for file discovery and loading
- Progress bars for all operations
- Memory-efficient data processing
- Fast file I/O with optimizations
- Comprehensive performance monitoring
"""

import os
import sys
import glob
from pathlib import Path
import warnings
import time
import json
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
from functools import partial
import psutil
import pandas as pd
import numpy as np
from tqdm import tqdm

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

# Performance settings
CHUNK_SIZE = 10000  # For chunked processing
MAX_WORKERS_SEARCH = None  # Use all cores for search
MAX_WORKERS_LOAD = 4  # Limit for I/O operations
PROGRESS_UPDATE_INTERVAL = 0.1  # Progress bar update frequency

def get_system_info():
    """Get system information for optimization"""
    return {
        'cpu_count': multiprocessing.cpu_count(),
        'memory_gb': psutil.virtual_memory().total / (1024**3),
        'available_memory_gb': psutil.virtual_memory().available / (1024**3)
    }

def optimized_walk_search(args):
    """
    Ultra-optimized file search with multiple patterns and smart termination.
    """
    search_dir, patterns, max_depth, max_files_per_pattern = args
    results = {pattern: [] for pattern in patterns}
    
    try:
        for root, dirs, files in os.walk(search_dir):
            # Calculate current depth
            current_depth = root[len(search_dir):].count(os.sep)
            if max_depth is not None and current_depth >= max_depth:
                dirs[:] = []  # Don't descend further
                continue
            
            # Skip hidden and system directories for speed
            dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']
            
            # Check files against all patterns simultaneously
            for file in files:
                if file.startswith('.'):  # Skip hidden files
                    continue
                    
                file_lower = file.lower()
                for pattern in patterns:
                    if pattern in file_lower and len(results[pattern]) < max_files_per_pattern:
                        full_path = os.path.join(root, file)
                        results[pattern].append(full_path)
            
            # Early termination if all patterns have enough files
            if all(len(files) >= max_files_per_pattern for files in results.values()):
                break
                
    except (PermissionError, OSError, UnicodeDecodeError):
        pass  # Skip problematic directories
    
    return results

def find_consensus_results_files_optimized(base_directory, max_workers=None):
    """
    Ultra-fast multiprocessing file discovery with intelligent search strategy.
    """
    print(f"\n🔍 ULTRA-OPTIMIZED FILE DISCOVERY")
    print(f"=" * 50)
    print(f"📁 Searching in: {base_directory}")
    
    if not os.path.exists(base_directory):
        print(f"❌ Directory not found: {base_directory}")
        return {}
    
    search_start = time.time()
    
    # Get system info and optimize worker count
    system_info = get_system_info()
    if max_workers is None:
        max_workers = min(system_info['cpu_count'], 16)  # Cap at 16 for I/O bound tasks
    
    print(f"💻 Using {max_workers} cores for parallel search")
    
    # Define search patterns - search all at once for efficiency
    search_patterns = [
        'final_results.csv',
        'results.csv',  # This will catch smina/gold/ledock results.csv
        'smina_results.csv', 
        'gold_results.csv',
        'ledock_results.csv',
        'combined_'  # For combined results
    ]
    
    # Get search directories with intelligent partitioning
    search_dirs = get_search_directories(base_directory, max_workers)
    max_depth = 4  # Slightly deeper search
    max_files_per_pattern = 10000  # Higher limit
    
    print(f"📂 Partitioned search into {len(search_dirs)} segments")
    
    # Parallel file search with all patterns at once
    all_results = {pattern: [] for pattern in search_patterns}
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(search_dirs), desc="🔍 Searching", 
                 bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}") as pbar:
            
            # Submit search jobs
            futures = []
            for search_dir in search_dirs:
                args = (search_dir, search_patterns, max_depth, max_files_per_pattern)
                future = executor.submit(optimized_walk_search, args)
                futures.append(future)
            
            # Collect results
            total_found = 0
            for future in as_completed(futures):
                try:
                    segment_results = future.result(timeout=60)  # Longer timeout
                    for pattern, files in segment_results.items():
                        all_results[pattern].extend(files)
                        total_found += len(files)
                    
                    pbar.set_postfix({'Found': f"{total_found:,}"})
                    pbar.update(1)
                except Exception as e:
                    print(f"⚠️ Search segment failed: {e}")
                    pbar.update(1)
    
    # Categorize results intelligently
    categorized_results = categorize_found_files(all_results)
    
    search_time = time.time() - search_start
    total_files = sum(len(files) for files in categorized_results.values())
    
    print(f"\n✅ File discovery completed in {search_time:.2f} seconds")
    print(f"📊 Files found by category:")
    for category, files in categorized_results.items():
        print(f"  {category}: {len(files):,} files")
    print(f"📋 Total files: {total_files:,}")
    
    if total_files > 0:
        print(f"⚡ Search rate: {total_files/search_time:,.0f} files/sec")
    
    return categorized_results


def get_search_directories(base_directory, max_workers):
    """
    Intelligently partition the search space for optimal parallelization.
    """
    search_dirs = []
    
    try:
        # Get immediate subdirectories
        items = os.listdir(base_directory)
        subdirs = [os.path.join(base_directory, item) for item in items 
                  if os.path.isdir(os.path.join(base_directory, item)) 
                  and not item.startswith('.')]
        
        if len(subdirs) >= max_workers:
            # Use subdirectories directly if we have enough
            search_dirs = subdirs[:max_workers * 2]  # Take a reasonable number
        else:
            # Add the base directory and all subdirectories
            search_dirs = [base_directory] + subdirs
        
        # If still not enough, look one level deeper
        if len(search_dirs) < max_workers and subdirs:
            for subdir in subdirs[:5]:  # Limit to avoid explosion
                try:
                    sub_items = os.listdir(subdir)
                    sub_subdirs = [os.path.join(subdir, item) for item in sub_items 
                                  if os.path.isdir(os.path.join(subdir, item)) 
                                  and not item.startswith('.')]
                    search_dirs.extend(sub_subdirs[:max_workers])
                except (PermissionError, OSError):
                    continue
                    
                if len(search_dirs) >= max_workers * 2:
                    break
    
    except (PermissionError, OSError):
        search_dirs = [base_directory]
    
    return search_dirs if search_dirs else [base_directory]


def categorize_found_files(all_results):
    """
    Intelligently categorize found files into result types.
    """
    categorized = {
        'final_results': [],
        'smina_results': [],
        'gold_results': [],
        'ledock_results': [],
        'combined_results': []
    }
    
    # Process all found files
    for pattern, files in all_results.items():
        for file_path in files:
            file_lower = file_path.lower()
            file_name = os.path.basename(file_lower)
            dir_path = os.path.dirname(file_lower)
            
            # Categorize based on path and filename
            if 'final_results' in file_name:
                categorized['final_results'].append(file_path)
            elif 'combined' in file_name and file_name.endswith('.csv'):
                categorized['combined_results'].append(file_path)
            elif 'smina' in dir_path and 'results.csv' in file_name:
                categorized['smina_results'].append(file_path)
            elif 'gold' in dir_path and 'results.csv' in file_name:
                categorized['gold_results'].append(file_path)
            elif 'ledock' in dir_path and 'results.csv' in file_name:
                categorized['ledock_results'].append(file_path)
            elif 'smina' in file_name:
                categorized['smina_results'].append(file_path)
            elif 'gold' in file_name:
                categorized['gold_results'].append(file_path)
            elif 'ledock' in file_name:
                categorized['ledock_results'].append(file_path)
    
    # Remove duplicates while preserving order
    for category in categorized:
        categorized[category] = list(dict.fromkeys(categorized[category]))
    
    return categorized

def optimized_file_loader(file_info):
    """
    Ultra-fast file loader with smart chunking and memory optimization.
    """
    file_path, expected_size_mb = file_info
    
    try:
        # Measure actual file size
        actual_size_mb = os.path.getsize(file_path) / (1024**2)
        
        # Choose loading strategy based on file size
        read_kwargs = {
            'dtype': {'drugbank_id': 'str', 'uniprot_id': 'str'},
            'na_values': ['', 'NA', 'N/A', 'null', 'NULL', 'None', 'NaN'],
            'keep_default_na': True
        }
        
        if actual_size_mb > 500:  # Very large files
            read_kwargs.update({
                'low_memory': False,
                'engine': 'c',
                'chunksize': None  # Load all at once for very large files
            })
        elif actual_size_mb > 50:  # Medium files
            read_kwargs.update({
                'low_memory': False,
                'engine': 'c'
            })
        else:  # Small files
            read_kwargs.update({
                'engine': 'c'
            })
        
        # Load the CSV
        df = pd.read_csv(file_path, **read_kwargs)
        
        if df.empty:
            return None
        
        # Add metadata efficiently
        metadata = {
            'source_file': os.path.basename(file_path),
            'source_dir': os.path.dirname(file_path),
            'file_size_mb': actual_size_mb
        }
        
        # Determine source type from path analysis
        path_lower = file_path.lower()
        if 'final_results' in path_lower:
            metadata['source_type'] = 'final_results'
        elif '/smina/' in path_lower or 'smina' in os.path.basename(path_lower):
            metadata['source_type'] = 'smina_results'
        elif '/gold/' in path_lower or 'gold' in os.path.basename(path_lower):
            metadata['source_type'] = 'gold_results'
        elif '/ledock/' in path_lower or 'ledock' in os.path.basename(path_lower):
            metadata['source_type'] = 'ledock_results'
        elif 'combined' in path_lower:
            metadata['source_type'] = 'combined_results'
        else:
            metadata['source_type'] = 'unknown'
        
        # Add metadata columns efficiently
        for key, value in metadata.items():
            df[key] = value
        
        # Extract uniprot_id and drugbank_id from directory path
        # Path format: .../DB00390_ATP1A1_Q13286_cavity_3/...
        dir_path = os.path.dirname(file_path)
        dir_name = os.path.basename(dir_path)
        
        # Parse the directory name format: {drugbank_id}_{gene_name}_{uniprot_id}_cavity_{cavity_number}
        import re
        match = re.match(r'(DB\d+)_([A-Z0-9]+)_([A-Z0-9]+)_cavity_(\d+)', dir_name)
        if match:
            drugbank_id, gene_name, uniprot_id, cavity_num = match.groups()
            df['drugbank_id'] = drugbank_id
            df['uniprot_id'] = uniprot_id
            df['gene_name'] = gene_name
            df['cavity_index'] = int(cavity_num)
        else:
            # Fallback: keep existing values or set to Unknown
            if 'drugbank_id' not in df.columns:
                df['drugbank_id'] = 'Unknown'
            if 'uniprot_id' not in df.columns:
                df['uniprot_id'] = 'Unknown'
        
        return {
            'df': df,
            'file_path': file_path,
            'records': len(df),
            'size_mb': actual_size_mb,
            'source_type': metadata['source_type']
        }
        
    except Exception as e:
        return {
            'df': None,
            'file_path': file_path,
            'error': str(e),
            'records': 0,
            'size_mb': 0
        }

def load_consensus_data_optimized(result_files, max_workers=6):
    """
    Ultra-high-performance data loading with intelligent batching and memory management.
    """
    print(f"\n📖 ULTRA-OPTIMIZED DATA LOADING")
    print(f"=" * 50)
    
    # Flatten and prioritize files
    prioritized_files = prioritize_files_for_loading(result_files)
    
    if not prioritized_files:
        print("❌ No files to load")
        return pd.DataFrame()
    
    print(f"📄 Loading {len(prioritized_files):,} files using {max_workers} processes")
    
    load_start = time.time()
    successful_chunks = []
    failed_files = []
    total_size_mb = 0
    total_records = 0
    
    # Estimate file sizes for better scheduling
    file_info_list = []
    for file_path in prioritized_files:
        try:
            size_mb = os.path.getsize(file_path) / (1024**2)
            file_info_list.append((file_path, size_mb))
            total_size_mb += size_mb
        except OSError:
            file_info_list.append((file_path, 0))
    
    print(f"💾 Total data size: {total_size_mb:.1f} MB")
    
    # Load files with progress tracking and error handling
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        with tqdm(total=len(file_info_list), desc="📖 Loading Files", 
                 bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}") as pbar:
            
            # Submit all loading jobs
            future_to_info = {
                executor.submit(optimized_file_loader, file_info): file_info
                for file_info in file_info_list
            }
            
            # Process results as they complete
            for future in as_completed(future_to_info):
                file_info = future_to_info[future]
                file_path = file_info[0]
                
                try:
                    result = future.result(timeout=120)  # 2 minute timeout per file
                    
                    if result and result.get('df') is not None and not result['df'].empty:
                        successful_chunks.append(result['df'])
                        total_records += result['records']
                        
                        # Update progress with detailed info
                        pbar.set_postfix({
                            'Files': len(successful_chunks),
                            'Records': f"{total_records:,}",
                            'Failed': len(failed_files),
                            'MB': f"{sum(r.get('size_mb', 0) for r in [result]):,.1f}"
                        })
                    else:
                        failed_files.append(file_path)
                        if result and 'error' in result:
                            print(f"⚠️ Failed: {os.path.basename(file_path)} - {result['error']}")
                    
                    pbar.update(1)
                    
                except Exception as e:
                    failed_files.append(file_path)
                    print(f"⚠️ Error loading {os.path.basename(file_path)}: {e}")
                    pbar.update(1)
    
    load_time = time.time() - load_start
    
    print(f"\n⚡ Loading phase completed in {load_time:.2f} seconds")
    print(f"✅ Successfully loaded: {len(successful_chunks)} files")
    print(f"❌ Failed to load: {len(failed_files)} files")
    print(f"📊 Total records loaded: {total_records:,}")
    
    if failed_files and len(failed_files) <= 10:
        print(f"⚠️  Failed files: {[os.path.basename(f) for f in failed_files]}")
    elif failed_files:
        print(f"⚠️  Failed files: {len(failed_files)} (showing first 5): {[os.path.basename(f) for f in failed_files[:5]]}")
    
    # Combine DataFrames efficiently
    if successful_chunks:
        print(f"\n🔄 Combining {len(successful_chunks)} DataFrames...")
        combine_start = time.time()
        
        # Sort chunks by size (largest first) for better memory usage
        successful_chunks.sort(key=lambda df: len(df), reverse=True)
        
        # Use concat with optimized settings
        consensus_data = pd.concat(
            successful_chunks, 
            ignore_index=True, 
            sort=False,
            copy=False  # Avoid unnecessary copying
        )
        
        combine_time = time.time() - combine_start
        total_time = load_time + combine_time
        
        # Memory optimization
        memory_mb = consensus_data.memory_usage(deep=True).sum() / (1024**2)
        
        print(f"✅ Data loading completed in {total_time:.2f} seconds")
        print(f"  📖 Loading: {load_time:.2f}s")
        print(f"  🔄 Combining: {combine_time:.2f}s")
        print(f"📊 Final dataset:")
        print(f"  📋 Records: {len(consensus_data):,}")
        print(f"  📊 Columns: {len(consensus_data.columns)}")
        print(f"  💾 Memory: {memory_mb:.1f} MB")
        print(f"⚡ Processing rate: {len(consensus_data)/total_time:,.0f} records/sec")
        
        return consensus_data
    else:
        print("❌ No valid data could be loaded")
        return pd.DataFrame()


def prioritize_files_for_loading(result_files):
    """
    Intelligently prioritize files for loading to optimize processing order.
    """
    prioritized = []
    
    # Priority order: final_results > combined > individual tools
    priority_order = [
        'final_results',
        'combined_results', 
        'smina_results',
        'gold_results',
        'ledock_results'
    ]
    
    for file_type in priority_order:
        if file_type in result_files and result_files[file_type]:
            files = result_files[file_type]
            
            # Sort files by likely importance (size, name, etc.)
            sorted_files = sorted(files, key=lambda f: (
                -os.path.getsize(f) if os.path.exists(f) else 0,  # Larger files first
                -len(os.path.basename(f)),  # Longer names tend to be more descriptive
                f  # Alphabetical as tiebreaker
            ))
            
            prioritized.extend(sorted_files)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_prioritized = []
    for file_path in prioritized:
        if file_path not in seen:
            seen.add(file_path)
            unique_prioritized.append(file_path)
    
    return unique_prioritized

def clean_and_prepare_data_optimized(df):
    """
    Ultra-fast data cleaning and preparation with optimizations.
    """
    if df.empty:
        return df
    
    print(f"\n🧹 OPTIMIZED DATA CLEANING")
    print(f"=" * 50)
    print(f"📊 Input data: {len(df):,} records, {len(df.columns)} columns")
    
    clean_start = time.time()
    
    # Make a copy to avoid modifying the original
    clean_df = df.copy()
    
    with tqdm(total=8, desc="🧹 Cleaning Steps", bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]") as pbar:
        
        # Step 1: Standardize column names
        column_standardization = {
            'drug_id': 'drugbank_id',
            'compound_id': 'drugbank_id',
            'protein_id': 'uniprot_id',
            'target_id': 'uniprot_id',
            'binding_score': 'score',
            'docking_score': 'score',
            'affinity': 'score',
            'binding_affinity': 'score',
        }
        
        for old_name, new_name in column_standardization.items():
            if old_name in clean_df.columns and new_name not in clean_df.columns:
                clean_df[new_name] = clean_df[old_name]
        pbar.update(1)
        
        # Step 2: Convert numeric columns (vectorized)
        numeric_columns = ['score', 'rmsd', 'energy', 'consensus_score']
        for col in numeric_columns:
            if col in clean_df.columns:
                clean_df[col] = pd.to_numeric(clean_df[col], errors='coerce')
        pbar.update(1)
        
        # Step 3: Clean identifiers (vectorized)
        identifier_columns = ['drugbank_id', 'uniprot_id', 'gene_name', 'drug_name', 'protein_name']
        for col in identifier_columns:
            if col in clean_df.columns:
                clean_df[col] = clean_df[col].astype(str).str.strip()
                clean_df[col] = clean_df[col].replace(['nan', 'Unknown', ''], pd.NA)
        
        # Validate extracted IDs
        if 'drugbank_id' in clean_df.columns and 'uniprot_id' in clean_df.columns:
            valid_drugbank = clean_df['drugbank_id'].str.startswith('DB', na=False)
            valid_uniprot = clean_df['uniprot_id'].str.match(r'^[A-Z0-9]{6}$', na=False)
            print(f"🔍 Valid DrugBank IDs: {valid_drugbank.sum():,}/{len(clean_df):,}")
            print(f"🔍 Valid UniProt IDs: {valid_uniprot.sum():,}/{len(clean_df):,}")
        pbar.update(1)
        
        # Step 4: Create categorical variables
        if 'docking_tool' in clean_df.columns:
            clean_df['docking_tool'] = clean_df['docking_tool'].astype('category')
        pbar.update(1)
        
        # Step 5: Add success/failure flags (vectorized)
        if 'score' in clean_df.columns:
            clean_df['docking_successful'] = clean_df['score'].notna()
            clean_df['high_affinity'] = clean_df['score'] < -7.0
        pbar.update(1)
        
        # Step 6: Extract and standardize tool information
        if 'source_type' in clean_df.columns:
            clean_df['primary_tool'] = clean_df['source_type'].str.replace('_results', '').str.title()
            tool_mapping = {
                'Smina': 'Smina',
                'Gold': 'Gold', 
                'Ledock': 'LeDock',
                'Final': 'Consensus'
            }
            clean_df['primary_tool'] = clean_df['primary_tool'].map(tool_mapping).fillna(clean_df['primary_tool'])
        pbar.update(1)
        
        # Step 7: Create compound-target pairs (vectorized)
        if all(col in clean_df.columns for col in ['drugbank_id', 'uniprot_id']):
            clean_df['compound_target_pair'] = clean_df['drugbank_id'].astype(str) + '_' + clean_df['uniprot_id'].astype(str)
        pbar.update(1)
        
        # Step 8: Ensure analysis workflow compatibility
        if 'score' in clean_df.columns and 'Score1' not in clean_df.columns:
            clean_df['Score1'] = clean_df['score']
        
        if 'Score2' not in clean_df.columns:
            if 'consensus_score' in clean_df.columns:
                clean_df['Score2'] = clean_df['consensus_score']
            elif 'energy' in clean_df.columns:
                clean_df['Score2'] = clean_df['energy']
            elif 'Score1' in clean_df.columns:
                clean_df['Score2'] = clean_df['Score1']
        
        if 'rmsd' in clean_df.columns and 'RMSD' not in clean_df.columns:
            clean_df['RMSD'] = clean_df['rmsd']
        elif 'RMSD' not in clean_df.columns:
            clean_df['RMSD'] = float('nan')
        pbar.update(1)
    
    clean_time = time.time() - clean_start
    final_memory_mb = clean_df.memory_usage(deep=True).sum() / (1024**2)
    
    print(f"✅ Cleaning completed in {clean_time:.2f} seconds")
    print(f"📊 Output data: {len(clean_df):,} records, {len(clean_df.columns)} columns")
    print(f"💾 Memory usage: {final_memory_mb:.1f} MB")
    
    return clean_df

def generate_data_summary(df):
    """Generate comprehensive data summary with performance optimization."""
    if df.empty:
        return {}
    
    summary = {}
    
    # Basic statistics (vectorized)
    summary['total_records'] = len(df)
    summary['unique_compounds'] = df['drugbank_id'].nunique() if 'drugbank_id' in df.columns else 0
    summary['unique_targets'] = df['uniprot_id'].nunique() if 'uniprot_id' in df.columns else 0
    summary['unique_pairs'] = df['compound_target_pair'].nunique() if 'compound_target_pair' in df.columns else 0
    
    # ID extraction validation
    if 'drugbank_id' in df.columns and 'uniprot_id' in df.columns:
        valid_drugbank = df['drugbank_id'].str.startswith('DB', na=False).sum()
        valid_uniprot = df['uniprot_id'].str.match(r'^[A-Z0-9]{6}$', na=False).sum()
        
        summary['id_extraction'] = {
            'valid_drugbank_ids': int(valid_drugbank),
            'valid_uniprot_ids': int(valid_uniprot),
            'extraction_success_rate': float(min(valid_drugbank, valid_uniprot) / len(df)) if len(df) > 0 else 0.0
        }
    
    # Analysis workflow compatibility
    required_cols = ['drugbank_id', 'uniprot_id', 'Score1', 'RMSD']
    summary['analysis_compatible'] = all(col in df.columns for col in required_cols)
    summary['missing_required_columns'] = [col for col in required_cols if col not in df.columns]
    
    # Success rates
    if 'docking_successful' in df.columns:
        summary['success_rate'] = df['docking_successful'].mean() * 100
        summary['successful_dockings'] = df['docking_successful'].sum()
    
    # Score statistics
    score_col = 'Score1' if 'Score1' in df.columns else 'score'
    if score_col in df.columns:
        scores = df[score_col].dropna()
        if not scores.empty:
            summary['score_stats'] = {
                'column_used': score_col,
                'mean': float(scores.mean()),
                'median': float(scores.median()),
                'std': float(scores.std()),
                'min': float(scores.min()),
                'max': float(scores.max()),
                'q25': float(scores.quantile(0.25)),
                'q75': float(scores.quantile(0.75)),
                'valid_count': int(len(scores)),
                'missing_count': int(df[score_col].isna().sum())
            }
    
    # RMSD statistics
    if 'RMSD' in df.columns:
        rmsd_values = df['RMSD'].dropna()
        if not rmsd_values.empty:
            summary['rmsd_stats'] = {
                'mean': float(rmsd_values.mean()),
                'median': float(rmsd_values.median()),
                'min': float(rmsd_values.min()),
                'max': float(rmsd_values.max()),
                'valid_count': int(len(rmsd_values)),
                'missing_count': int(df['RMSD'].isna().sum())
            }
    
    # Tool distribution
    if 'primary_tool' in df.columns:
        summary['tool_distribution'] = df['primary_tool'].value_counts().to_dict()
    
    return summary

def export_prepared_data_optimized(df, output_dir="./"):
    """
    High-performance data export with progress tracking.
    """
    if df.empty:
        print("❌ No data to export")
        return {}
    
    print(f"\n💾 OPTIMIZED DATA EXPORT")
    print(f"=" * 50)
    
    # Create output directory if needed
    if output_dir != "./":
        os.makedirs(output_dir, exist_ok=True)
    
    export_start = time.time()
    export_df = df.copy()
    
    # Map columns for analysis compatibility
    column_mapping = {
        'score': 'Score1',
        'binding_affinity': 'Score1',
        'rmsd': 'RMSD',
    }
    
    for old_col, new_col in column_mapping.items():
        if old_col in export_df.columns and new_col not in export_df.columns:
            export_df[new_col] = export_df[old_col]
    
    # Ensure required columns exist
    required_columns = ['drugbank_id', 'uniprot_id', 'Score1', 'RMSD']
    for col in required_columns:
        if col not in export_df.columns:
            if col == 'Score1':
                export_df[col] = export_df.get('score', float('nan'))
            elif col == 'RMSD':
                export_df[col] = export_df.get('rmsd', float('nan'))
            elif col == 'drugbank_id':
                # Try to extract from source_dir if not already extracted
                if 'source_dir' in export_df.columns:
                    import re
                    def extract_drugbank_id(source_dir):
                        if pd.isna(source_dir):
                            return 'Unknown'
                        dir_name = os.path.basename(str(source_dir))
                        match = re.match(r'(DB\d+)_[A-Z0-9]+_[A-Z0-9]+_cavity_\d+', dir_name)
                        return match.group(1) if match else 'Unknown'
                    export_df[col] = export_df['source_dir'].apply(extract_drugbank_id)
                else:
                    export_df[col] = 'Unknown'
            elif col == 'uniprot_id':
                # Try to extract from source_dir if not already extracted
                if 'source_dir' in export_df.columns:
                    import re
                    def extract_uniprot_id(source_dir):
                        if pd.isna(source_dir):
                            return 'Unknown'
                        dir_name = os.path.basename(str(source_dir))
                        match = re.match(r'DB\d+_[A-Z0-9]+_([A-Z0-9]+)_cavity_\d+', dir_name)
                        return match.group(1) if match else 'Unknown'
                    export_df[col] = export_df['source_dir'].apply(extract_uniprot_id)
                else:
                    export_df[col] = 'Unknown'
            else:
                export_df[col] = 'Unknown'
    
    # Add Score2 if missing
    if 'Score2' not in export_df.columns:
        if 'consensus_score' in export_df.columns:
            export_df['Score2'] = export_df['consensus_score']
        elif 'energy' in export_df.columns:
            export_df['Score2'] = export_df['energy']
        else:
            export_df['Score2'] = export_df['Score1']
    
    # File paths
    csv_path = os.path.join(output_dir, "combined_consensus_docking_results.csv")
    parquet_path = os.path.join(output_dir, "combined_consensus_docking_results.parquet")
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    backup_csv = os.path.join(output_dir, f"backup_consensus_data_{timestamp}.csv")
    summary_path = os.path.join(output_dir, f"data_preparation_summary_{timestamp}.json")
    
    # Export with progress tracking
    with tqdm(total=4, desc="💾 Exporting Files", bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]") as pbar:
        
        # Export CSV
        export_df.to_csv(csv_path, index=False)
        pbar.set_postfix({'File': 'CSV'})
        pbar.update(1)
        
        # Export Parquet
        export_df.to_parquet(parquet_path, index=False)
        pbar.set_postfix({'File': 'Parquet'})
        pbar.update(1)
        
        # Export backup
        export_df.to_csv(backup_csv, index=False)
        pbar.set_postfix({'File': 'Backup'})
        pbar.update(1)
        
        # Export summary
        summary = generate_data_summary(export_df)
        summary['compatibility'] = {
            'analysis_workflow': 'compatible',
            'expected_files_created': [
                'combined_consensus_docking_results.csv',
                'combined_consensus_docking_results.parquet'
            ],
            'required_columns_present': required_columns,
            'column_mappings_applied': column_mapping
        }
        
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        pbar.set_postfix({'File': 'Summary'})
        pbar.update(1)
    
    export_time = time.time() - export_start
    
    print(f"✅ Export completed in {export_time:.2f} seconds")
    print(f"📄 Files created:")
    print(f"  CSV: {csv_path}")
    print(f"  Parquet: {parquet_path}")
    print(f"  Backup: {backup_csv}")
    print(f"  Summary: {summary_path}")
    
    return {
        'csv_path': csv_path,
        'parquet_path': parquet_path,
        'backup_path': backup_csv,
        'summary_path': summary_path,
        'compatible_with_analysis': True
    }

def main():
    """
    Main execution function with comprehensive performance monitoring.
    """
    total_start_time = time.time()
    
    print("🚀 HIGH-PERFORMANCE CONSENSUS DOCKING DATA PARSER")
    print("=" * 60)
    
    system_info = get_system_info()
    print(f"💻 System: {system_info['cpu_count']} cores, {system_info['memory_gb']:.1f}GB RAM")
    print(f"💾 Available: {system_info['available_memory_gb']:.1f}GB")
    
    # Step 1: File Discovery
    result_files = find_consensus_results_files_optimized(
        base_directory="/media/onur/Elements/cavity_space_consensus_docking",
        max_workers=min(multiprocessing.cpu_count(), 12)  # Cap workers for I/O
    )
    
    total_files = sum(len(files) for files in result_files.values())
    if total_files == 0:
        print("❌ No consensus docking files found. Exiting.")
        return
    
    # Step 2: Data Loading
    consensus_data = load_consensus_data_optimized(
        result_files, 
        max_workers=min(8, system_info['cpu_count'])  # Optimal for I/O bound loading
    )
    
    if consensus_data.empty:
        print("❌ No data could be loaded. Exiting.")
        return
    
    # Step 3: Data Cleaning
    cleaned_data = clean_and_prepare_data_optimized(consensus_data)
    
    # Step 4: Data Export
    export_results = export_prepared_data_optimized(cleaned_data, output_dir="./")
    
    # Final Performance Summary
    total_time = time.time() - total_start_time
    
    print(f"\n🎉 PROCESSING COMPLETE")
    print(f"=" * 60)
    print(f"⏱️  Total execution time: {total_time:.2f} seconds")
    print(f"📊 Total records processed: {len(cleaned_data):,}")
    print(f"⚡ Processing rate: {len(cleaned_data)/total_time:,.0f} records/sec")
    print(f"💾 Final memory usage: {cleaned_data.memory_usage(deep=True).sum() / (1024**2):.1f} MB")
    
    if export_results.get('compatible_with_analysis'):
        print(f"\n✅ Files ready for analysis workflow!")
        print(f"🚀 You can now run the Jupyter notebook analysis.")
    
    return cleaned_data, export_results

if __name__ == "__main__":
    try:
        cleaned_data, export_results = main()
    except KeyboardInterrupt:
        print("\n⚠️  Process interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)