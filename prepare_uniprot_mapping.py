#!/usr/bin/env python3
"""
Script to query UniProt for gene name to UniProt ID mapping.
This script creates a mapping file that can be used by the docking workflow.

Usage:
    python prepare_uniprot_mapping.py

Output:
    - uniprot_gene_mapping.csv: Contains gene names and their corresponding UniProt IDs
    - uniprot_preparation.log: Log file for this process
"""

import os
import sys
import pandas as pd
import math
import time
import logging
from pathlib import Path
from tqdm import tqdm

# --- Configuration ---
DRUG_TO_PROTEIN_TSV = "/opt/data/multiscale_interactome_data/1_drug_to_protein.tsv"
UNIPROT_MAPPING_CSV = "uniprot_gene_mapping.csv"
LOG_FILE = "uniprot_preparation.log"

# Set up logging
logging.basicConfig(filename=LOG_FILE, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logging.getLogger().addHandler(console_handler)

def get_uniprot_data(gene_names_list, out_csv):
    """
    Queries UniProt for UniProt IDs for a list of human gene names.
    Each gene is queried individually to ensure accurate results.
    Saves the mapping to a CSV.
    
    Args:
        gene_names_list (list): List of gene names to query
        out_csv (str): Output CSV file path
        
    Returns:
        pd.DataFrame: DataFrame with UniProt mapping data
    """
    # Remove duplicates and ensure we only query each gene once
    unique_gene_names = list(set(gene_names_list))
    logging.info(f"Fetching UniProt IDs for {len(unique_gene_names)} unique gene names...")
    logging.info(f"(Removed {len(gene_names_list) - len(unique_gene_names)} duplicate gene names)")
    
    try:
        from bioservices import UniProt
    except ImportError:
        logging.critical("bioservices package not found. Please install it with: pip install bioservices")
        sys.exit(1)
    
    service = UniProt()
    
    dfs = []
    successful_queries = 0
    failed_queries = 0
    genes_found = set()
    genes_not_found = set()
    
    # Query each gene individually
    for i, gene_name in enumerate(tqdm(unique_gene_names, desc="Querying UniProt", unit="gene")):
        try:
            # Query with just the gene name (no special formatting)
            df = service.get_df(gene_name, columns="Entry,Entry Name,Gene Names,Organism")
            
            if df is not None and not df.empty:
                # Filter to ensure we only keep human entries
                df_filtered = df[df['Organism'].str.contains("Homo sapiens", na=False)].copy()
                
                if not df_filtered.empty:
                    # Check if this gene is actually mentioned in the Gene Names
                    found_for_this_gene = False
                    for gene_name_entry in df_filtered['Gene Names'].dropna():
                        gene_names_in_entry = str(gene_name_entry).split()
                        if gene_name in gene_names_in_entry:
                            found_for_this_gene = True
                            break
                    
                    if found_for_this_gene:
                        # Map 'Gene Names' column to the queried gene name for easier lookup
                        df_filtered['Mapped_Gene_Name'] = gene_name
                        dfs.append(df_filtered)
                        genes_found.add(gene_name)
                        successful_queries += 1
                        logging.info(f"Gene {i+1}/{len(unique_gene_names)}: Found {len(df_filtered)} entries for '{gene_name}'")
                    else:
                        genes_not_found.add(gene_name)
                        logging.warning(f"Gene {i+1}/{len(unique_gene_names)}: No matching entries for '{gene_name}' (found other genes)")
                        failed_queries += 1
                else:
                    genes_not_found.add(gene_name)
                    logging.warning(f"Gene {i+1}/{len(unique_gene_names)}: No human entries found for '{gene_name}'")
                    failed_queries += 1
            else:
                genes_not_found.add(gene_name)
                logging.warning(f"Gene {i+1}/{len(unique_gene_names)}: No data returned for '{gene_name}'")
                failed_queries += 1
                
            # Be respectful to the UniProt server
            time.sleep(0.5)  # Shorter sleep since we're doing individual queries
            
        except Exception as e:
            genes_not_found.add(gene_name)
            logging.error(f"Error querying UniProt for gene '{gene_name}': {e}")
            tqdm.write(f"Error querying UniProt for gene '{gene_name}': {e}")
            failed_queries += 1
            time.sleep(2)  # Longer sleep for errors

    # Handle case where no data was retrieved
    if not dfs:
        logging.error("No UniProt data retrieved for any gene.")
        # Create empty DataFrame with proper columns
        empty_df = pd.DataFrame(columns=['Entry', 'Entry Name', 'Gene Names', 'Mapped_Gene_Name', 'Organism'])
        # Still save the empty file for consistency
        empty_df.to_csv(out_csv, index=False)
        logging.info(f"Empty mapping file saved to '{out_csv}'")
        return empty_df

    # Combine all dataframes and remove duplicates
    try:
        final_df = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=['Entry']).reset_index(drop=True)
    except ValueError as e:
        logging.error(f"Error combining dataframes: {e}")
        # Create empty DataFrame with proper columns if concat fails
        final_df = pd.DataFrame(columns=['Entry', 'Entry Name', 'Gene Names', 'Mapped_Gene_Name', 'Organism'])
    
    # Save the complete mapping
    final_df.to_csv(out_csv, index=False)
    logging.info(f"UniProt mapping saved to '{out_csv}'")
    logging.info(f"Total unique UniProt entries found: {len(final_df)}")
    logging.info(f"Successful individual gene queries: {successful_queries}")
    logging.info(f"Failed individual gene queries: {failed_queries}")
    
    # Create a detailed summary of the mapping
    logging.info(f"=== MAPPING SUMMARY ===")
    logging.info(f"Total unique genes queried: {len(unique_gene_names)}")
    logging.info(f"Genes successfully mapped: {len(genes_found)}")
    logging.info(f"Genes not found: {len(genes_not_found)}")
    logging.info(f"Success rate: {len(genes_found)/len(unique_gene_names)*100:.1f}%")
    
    if genes_found:
        logging.info(f"Sample of mapped genes: {sorted(list(genes_found))[:10]}{'...' if len(genes_found) > 10 else ''}")
    
    if genes_not_found:
        logging.warning(f"Genes not found in UniProt:")
        # Log all unmapped genes, but limit to avoid excessive output
        genes_not_found_list = sorted(list(genes_not_found))
        for i in range(0, len(genes_not_found_list), 20):
            batch = genes_not_found_list[i:i+20]
            logging.warning(f"  {', '.join(batch)}")
    
    return final_df

def validate_input_files():
    """Validate that required input files exist."""
    if not os.path.exists(DRUG_TO_PROTEIN_TSV):
        logging.critical(f"Error: Drug-to-protein TSV file not found at '{DRUG_TO_PROTEIN_TSV}'. Please correct the path.")
        return False
    return True

def generate_summary_report(uniprot_df, original_gene_count, output_file="uniprot_mapping_summary.txt"):
    """
    Generate a detailed summary report of the UniProt mapping results.
    
    Args:
        uniprot_df (pd.DataFrame): The UniProt mapping DataFrame
        original_gene_count (int): Number of original gene names queried
        output_file (str): Path to save the summary report
    """
    try:
        with open(output_file, 'w') as f:
            f.write("UniProt Mapping Summary Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Query Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total genes queried: {original_gene_count}\n")
            f.write(f"UniProt entries found: {len(uniprot_df)}\n")
            
            if not uniprot_df.empty:
                mapped_genes = uniprot_df['Mapped_Gene_Name'].dropna().unique()
                f.write(f"Unique genes mapped: {len(mapped_genes)}\n")
                f.write(f"Success rate: {len(mapped_genes)/original_gene_count*100:.1f}%\n\n")
                
                f.write("Sample of mapped genes:\n")
                for gene in sorted(mapped_genes)[:20]:
                    f.write(f"  - {gene}\n")
                if len(mapped_genes) > 20:
                    f.write(f"  ... and {len(mapped_genes) - 20} more\n")
            else:
                f.write("Success rate: 0.0%\n")
                f.write("No genes were successfully mapped.\n")
            
            f.write(f"\nDetailed results saved in: {UNIPROT_MAPPING_CSV}\n")
            f.write(f"Log file: {LOG_FILE}\n")
        
        logging.info(f"Summary report saved to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error generating summary report: {e}")

def main():
    """Main execution function."""
    logging.info("Starting UniProt mapping preparation script.")
    logging.info(f"Drug-to-Protein TSV: {DRUG_TO_PROTEIN_TSV}")
    logging.info(f"Output UniProt Mapping CSV: {UNIPROT_MAPPING_CSV}")
    
    # Validate input files
    if not validate_input_files():
        sys.exit(1)
    
    # Check if mapping already exists
    uniprot_map_file = Path(UNIPROT_MAPPING_CSV)
    if uniprot_map_file.exists():
        age_days = (time.time() - uniprot_map_file.stat().st_mtime) / (24 * 3600)
        logging.info(f"Existing UniProt mapping file found (age: {age_days:.1f} days)")
        
        # Ask user if they want to regenerate
        try:
            response = input("Do you want to regenerate the UniProt mapping? (y/N): ").strip().lower()
            if response not in ['y', 'yes']:
                logging.info("Using existing mapping file. Script completed.")
                return
        except KeyboardInterrupt:
            logging.info("Operation cancelled by user.")
            return
    
    try:
        # Extract gene names from drug-protein data
        logging.info("Extracting gene names from drug-protein interaction data...")
        drug_protein_data = pd.read_csv(DRUG_TO_PROTEIN_TSV, sep='\t')
        
        # Get all gene names (including duplicates) for counting purposes
        all_gene_names = drug_protein_data['node_2_name'].tolist()
        unique_gene_names = drug_protein_data['node_2_name'].unique().tolist()
        
        logging.info(f"Found {len(all_gene_names)} total gene entries in '{DRUG_TO_PROTEIN_TSV}'.")
        logging.info(f"Found {len(unique_gene_names)} unique gene names (removed {len(all_gene_names) - len(unique_gene_names)} duplicates).")
        
        # Get UniProt mapping
        logging.info("Starting UniProt query process...")
        uniprot_df = get_uniprot_data(unique_gene_names, UNIPROT_MAPPING_CSV)
        
        if uniprot_df.empty:
            logging.warning("No UniProt data was retrieved for any gene names.")
            logging.warning("This could be due to:")
            logging.warning("  1. Network issues")
            logging.warning("  2. Gene names not being found in UniProt")
            logging.warning("  3. Issues with the UniProt service")
            logging.warning("An empty mapping file has been created. You may want to:")
            logging.warning("  - Check your gene names manually in UniProt")
            logging.warning("  - Try running the script again later")
            logging.warning("  - Check the log file for more details")
        else:
            logging.info("UniProt mapping preparation completed successfully!")
            logging.info(f"Retrieved data for {len(uniprot_df)} UniProt entries")
        
        logging.info(f"Mapping file saved as: {UNIPROT_MAPPING_CSV}")
        logging.info(f"Log file saved as: {LOG_FILE}")
        
        # Generate summary report
        generate_summary_report(uniprot_df, len(unique_gene_names))
        
    except Exception as e:
        logging.critical(f"Error during UniProt mapping preparation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
