#!/usr/bin/env python3
"""
Create IC50 Mapping Table for DrugBank ID - UniProt ID Interactions

This script integrates data from multiple sources to create a comprehensive mapping
table of IC50 values for drug-target interactions using DrugBank IDs and UniProt IDs.

Data Sources:
1. TTD Target Download: Maps TTD Target IDs to UniProt Entry Names
2. UniProt Gene Mapping: Maps UniProt Entry Names to UniProt Entry IDs
3. DrugBank Links: Maps TTD Drug IDs to DrugBank IDs
4. TTD Activity Data: Contains IC50 values for TTD Target-Drug pairs

Output:
A CSV file with columns: drugbank_id, uniprot_id, ttd_target_id, ttd_drug_id, 
                         pubchem_cid, activity, ic50_value, ic50_unit
"""

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict
import polars as pl

# File paths
TTD_TARGET_FILE = "/opt/data/ttd/P1-01-TTD_target_download.txt"
TTD_ACTIVITY_FILE = "/opt/data/ttd/P1-09-Target_compound_activity.txt"
DRUGBANK_LINKS_FILE = "/opt/data/drugbank/small_molecule_drug_links.csv"
UNIPROT_MAPPING_FILE = "/media/onur/Elements/cavity_space_consensus_docking/2025_06_29_batch_dock/uniprot_gene_mapping.csv"
OUTPUT_FILE = "/media/onur/Elements/cavity_space_consensus_docking/2025_06_29_batch_dock/ic50_mapping.csv"


def parse_ttd_targets(filepath):
    """
    Parse TTD target download file to extract TTD Target ID -> UniProt Entry Name mapping.
    
    The file format is:
    TTD_TARGET_ID   FIELD_TYPE   VALUE
    
    Returns:
        dict: {ttd_target_id: uniprot_entry_name}
    """
    print(f"📖 Parsing TTD targets from: {filepath}")
    
    ttd_to_uniprot = {}
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('TTD -') or line.startswith('Title') or \
               line.startswith('Version') or line.startswith('Provided') or \
               line.startswith('Any question') or line.startswith('Dr.') or \
               line.startswith('---') or line.startswith('Abbreviations'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                ttd_id = parts[0]
                field_type = parts[1]
                value = parts[2] if len(parts) > 2 else ''
                
                # Look for UNIPROID field
                if field_type == 'UNIPROID':
                    ttd_to_uniprot[ttd_id] = value
    
    print(f"   Found {len(ttd_to_uniprot):,} TTD targets with UniProt mappings")
    return ttd_to_uniprot


def load_uniprot_mapping(filepath):
    """
    Load UniProt mapping file to convert Entry Names to Entry IDs.
    
    Returns:
        dict: {entry_name: entry_id}
    """
    print(f"📖 Loading UniProt mapping from: {filepath}")
    
    # Use polars for efficient CSV reading
    df = pl.read_csv(filepath)
    
    # Create mapping from Entry Name to Entry ID
    mapping = dict(zip(df['Entry Name'].to_list(), df['Entry'].to_list()))
    
    print(f"   Loaded {len(mapping):,} UniProt mappings")
    return mapping


def load_drugbank_pubchem_mapping(drugbank_file):
    """
    Load mapping from DrugBank file to get PubChem Compound ID -> DrugBank ID mapping.
    """
    pubchem_to_drugbank = {}
    
    with open(drugbank_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            drugbank_id = row.get('DrugBank ID', '').strip()
            pubchem_id = row.get('PubChem Compound ID', '').strip()
            
            if drugbank_id and pubchem_id:
                # Store the mapping (PubChem ID may map to multiple DrugBank IDs, keep first)
                if pubchem_id not in pubchem_to_drugbank:
                    pubchem_to_drugbank[pubchem_id] = drugbank_id
    
    return pubchem_to_drugbank


def parse_ic50_value(activity_string):
    """
    Parse IC50/Ki value from activity string.
    
    Examples:
        'IC50 = 3270000 nM' -> (3270000.0, 'nM', 'IC50', '=')
        'IC50 > 10000 nM' -> (10000.0, 'nM', 'IC50', '>')
        'Ki = 2390 nM' -> (2390.0, 'nM', 'Ki', '=')
    
    Returns:
        tuple: (value, unit, measurement_type, operator) or (None, None, None, None)
    """
    # Pattern to match IC50/Ki values with operators
    # Matches: IC50 = 123 nM, IC50 > 456 uM, Ki = 789 nM, etc.
    pattern = r'(IC50|Ki)\s*([>=<]+)\s*([\d.]+)\s*([nuµμ]?M)'
    
    match = re.search(pattern, activity_string, re.IGNORECASE)
    if match:
        measurement_type = match.group(1)
        operator = match.group(2)
        value = float(match.group(3))
        unit = match.group(4)
        
        return value, unit, measurement_type, operator
    
    return None, None, None, None


def parse_activity_data(filepath, ttd_target_to_uniprot, uniprot_name_to_id, pubchem_to_drugbank):
    """
    Parse the TTD activity file and create final mapping table.
    
    Args:
        filepath: Path to P1-09-Target_compound_activity.txt
        ttd_target_to_uniprot: Mapping from TTD Target ID to UniProt Entry Name
        uniprot_name_to_id: Mapping from UniProt Entry Name to UniProt Entry ID
        pubchem_to_drugbank: Mapping from PubChem Compound ID to DrugBank ID
    
    Returns:
        list of dicts with keys: drugbank_id, uniprot_id, activity, pubchem_cid
    """
    print(f"📖 Parsing activity data from: {filepath}")
    
    results = []
    total_lines = 0
    mapped_lines = 0
    skipped_no_drugbank = 0
    skipped_no_uniprot = 0
    skipped_no_ic50 = 0
    
    with open(filepath, 'r', encoding='utf-8') as f:
        # Skip header line
        next(f)
        
        for line in f:
            total_lines += 1
            parts = line.strip().split('\t')
            
            if len(parts) < 4:
                continue
            
            ttd_target_id = parts[0]
            ttd_drug_id = parts[1]
            pubchem_cid = parts[2] if len(parts) > 2 else ''
            activity = parts[3] if len(parts) > 3 else ''
            
            # Parse IC50 value
            ic50_value, ic50_unit, measurement_type, operator = parse_ic50_value(activity)
            
            if ic50_value is None:
                skipped_no_ic50 += 1
                continue
            
            # Map TTD Target ID to UniProt Entry Name
            uniprot_entry_name = ttd_target_to_uniprot.get(ttd_target_id)
            if not uniprot_entry_name:
                skipped_no_uniprot += 1
                continue
            
            # Map UniProt Entry Name to UniProt Entry ID
            uniprot_id = uniprot_name_to_id.get(uniprot_entry_name)
            if not uniprot_id:
                skipped_no_uniprot += 1
                continue
            
            # Map PubChem CID to DrugBank ID
            drugbank_id = pubchem_to_drugbank.get(pubchem_cid)
            if not drugbank_id:
                skipped_no_drugbank += 1
                continue
            
            # Successfully mapped!
            results.append({
                'drugbank_id': drugbank_id,
                'uniprot_id': uniprot_id,
                'ttd_target_id': ttd_target_id,
                'ttd_drug_id': ttd_drug_id,
                'pubchem_cid': pubchem_cid,
                'activity': activity,
                'measurement_type': measurement_type,
                'operator': operator,
                'ic50_value': ic50_value,
                'ic50_unit': ic50_unit
            })
            mapped_lines += 1
    
    print(f"\n📊 Activity Data Processing Summary:")
    print(f"   Total activity records: {total_lines:,}")
    print(f"   Successfully mapped: {mapped_lines:,} ({mapped_lines/total_lines*100:.1f}%)")
    print(f"   Skipped - no IC50 value: {skipped_no_ic50:,}")
    print(f"   Skipped - no DrugBank mapping: {skipped_no_drugbank:,}")
    print(f"   Skipped - no UniProt mapping: {skipped_no_uniprot:,}")
    
    return results


def write_results(results, output_file):
    """Write results to CSV file."""
    print(f"\n💾 Writing results to: {output_file}")
    
    if not results:
        print("   ⚠️  No results to write!")
        return
    
    # Define column order
    columns = ['drugbank_id', 'uniprot_id', 'ttd_target_id', 'ttd_drug_id', 
               'pubchem_cid', 'activity', 'measurement_type', 'operator', 
               'ic50_value', 'ic50_unit']
    
    # Convert to Polars DataFrame for efficient writing
    df = pl.DataFrame(results)
    
    # Reorder columns
    df = df.select(columns)
    
    # Write CSV
    df.write_csv(output_file)
    
    print(f"   ✅ Wrote {len(results):,} records to {output_file}")
    
    # Print sample records
    print(f"\n📋 Sample records (first 5):")
    print(df.head(5))
    
    # Print statistics
    print(f"\n📊 Dataset Statistics:")
    print(f"   Unique DrugBank IDs: {df['drugbank_id'].n_unique():,}")
    print(f"   Unique UniProt IDs: {df['uniprot_id'].n_unique():,}")
    print(f"   Unique drug-target pairs: {df.select(['drugbank_id', 'uniprot_id']).unique().height:,}")
    
    # IC50 value statistics
    print(f"\n📈 IC50 Value Statistics:")
    ic50_stats = df.select([
        pl.col('ic50_value').min().alias('min'),
        pl.col('ic50_value').max().alias('max'),
        pl.col('ic50_value').mean().alias('mean'),
        pl.col('ic50_value').median().alias('median')
    ])
    print(ic50_stats)
    
    # Measurement type distribution
    print(f"\n🔬 Measurement Type Distribution:")
    print(df.group_by('measurement_type').agg(pl.len()).sort('measurement_type'))
    
    # Operator distribution
    print(f"\n🔢 Operator Distribution:")
    print(df.group_by('operator').agg(pl.len()).sort('operator'))
    
    # Unit distribution
    print(f"\n📏 Unit Distribution:")
    print(df.group_by('ic50_unit').agg(pl.len()).sort('ic50_unit'))


def main():
    """Main execution function."""
    print("=" * 80)
    print("🧬 IC50 MAPPING TABLE CREATION")
    print("=" * 80)
    print("\nIntegrating data from multiple sources to create DrugBank-UniProt IC50 mappings\n")
    
    # Step 1: Parse TTD targets to get TTD ID -> UniProt Entry Name mapping
    ttd_to_uniprot = parse_ttd_targets(TTD_TARGET_FILE)
    
    # Step 2: Load UniProt mapping to get Entry Name -> Entry ID mapping
    uniprot_name_to_id = load_uniprot_mapping(UNIPROT_MAPPING_FILE)
    
    # Step 3: Load DrugBank-PubChem mapping
    pubchem_to_drugbank = load_drugbank_pubchem_mapping(DRUGBANK_LINKS_FILE)
    print(f"   Found {len(pubchem_to_drugbank):,} PubChem IDs with DrugBank mappings")
    
    # Step 4: Parse activity data and create final mapping
    results = parse_activity_data(
        TTD_ACTIVITY_FILE, 
        ttd_to_uniprot, 
        uniprot_name_to_id, 
        pubchem_to_drugbank
    )
    
    # Step 5: Write results to file
    write_results(results, OUTPUT_FILE)
    
    print("\n" + "=" * 80)
    print("✅ IC50 mapping creation complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
