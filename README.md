# Consensus Docking Analysis Workflow

This repository contains scripts for running batch consensus docking experiments using AlphaFold structures and cavity data. The workflow processes **small molecule drugs only**.

## Workflow Overview

The docking workflow consists of the following main steps:

1. **UniProt Mapping Generation** (`prepare_uniprot_mapping.py`)
2. **Cavity Extraction** (`extract_cavities.py`)
3. **Required Structures Identification** (`identify_required_structures.py`)
4. **Negative Sample Generation** (`generate_negative_samples.py`) - Optional for ML training
5. **AlphaFold Model Extraction** (`extract_alphafold_models.py`)
6. **PDB Fixing** (`fix_required_pdbs.py`)
7. **PDB to PDBQT Conversion** (`convert_pdb_to_pdbqt.py`)
8. **Batch Docking Execution** (`batch_dock.py`)

## Scripts

### 1. prepare_uniprot_mapping.py
Generates UniProt ID to gene name mappings required for the docking workflow.

**Output**: `uniprot_gene_mapping.csv`

### 2. extract_cavities.py
Extracts receptor and cavity PDB files from AlphaFold cavity tarballs. This is a separate preprocessing step that can take significant time but only needs to be run once.

**Features**:
- Handles multiple cavities per protein (cavity_1, cavity_2, etc.)
- Creates mapping CSV for fast loading by main docking script
- Extracts files to `extracted_cavities/` folder in current working directory
- Can resume from existing extractions

**Configuration**: Edit the `CAVITY_TARBALL_FOLDER` path in the script

**Usage**:
```bash
python extract_cavities.py
```

**Output**: 
- `cavity_mapping.csv` - Mapping file for main docking script
- `extracted_cavities/` - Directory with extracted PDB files
- `cavity_extraction.log` - Extraction log

### 3. identify_required_structures.py
Identifies which AlphaFold structures are actually needed based on drug-protein interactions and ligand availability. This allows selective processing of only the structures that will be used for docking, dramatically improving efficiency.

**Features**:
- **Selective identification**: Only identifies structures needed for docking
- **Drug-protein matching**: Based on interaction data and available ligands
- **Priority scoring**: Ranks targets by interaction count
- **Mapping integration**: Works with all downstream scripts
- **Significant speedup**: Reduces processing time for large datasets

**Configuration**: Update file paths in script for drug-protein data and ligand folders

**Usage**:
```bash
python identify_required_structures.py
```

**Output**:
- `required_structures.csv` - List of required structures for docking
- `required_structures.log` - Processing log

### 4. generate_negative_samples.py (Optional)
Generates negative samples (non-interacting drug-target pairs) to complement positive samples. This is essential for machine learning model training and evaluation, enabling proper assessment of docking quality.

**Purpose**: Creates balanced datasets with both positive (known interactions) and negative (non-interactions) samples for downstream analysis.

**Strategy**: Implements **Balanced Negative Sampling** based on Najm et al.
- **Key principle**: Each drug and target appears in equal numbers of positive and negative samples
- **Prevents spurious learning**: Model won't learn that "Drug X interacts with everything" or "Target Y never interacts"
- **Greedy algorithm**: Efficiently satisfies both drug and target balance constraints simultaneously
- **Example**: A drug with 20 positive targets gets 20 negative targets; a drug with 1 positive target gets 1 negative target

**Features**:
- **Balanced representation**: Maintains equal positive/negative counts per drug and per target
- **Intelligent validation**: Ensures no known interactions are included as negatives
- **Reproducible**: Seeded random generation for consistent results
- **Resource checking**: Validates ligand SDF files and protein structure availability
- **Comprehensive reporting**: Detailed balance statistics and underbalanced entity warnings
- **Augments existing data**: Works with `required_structures.csv` from step 3

**Configuration**: Edit file paths in script:
- `DRUG_TO_PROTEIN_TSV`: Known drug-protein interactions to avoid as negatives
- `SMALL_MOLECULE_DRUGS_CSV`: Available DrugBank drugs
- `CAVITY_MAPPING_CSV`: Available protein structures/pockets
- `PROCESSED_LIGAND_SDF_FOLDER`: Location of ligand SDF files

**Usage**:
```bash
python generate_negative_samples.py
```

**Output**:
- `required_structures_with_negatives.csv` - Augmented structure list with positive and negative samples
- `negative_samples_metadata.json` - Statistics and balance metrics
- `negative_samples.log` - Processing log with balance statistics

**Note**: Downstream scripts (steps 5-8) will automatically use `required_structures_with_negatives.csv` if it exists, otherwise fall back to `required_structures.csv`. This allows seamless integration without modifying existing scripts.

### 5. extract_alphafold_models.py
Extracts full AlphaFold models for required structures only (if step 3 was run). This extracts structures directly from the local AlphaFold database archive. Now optimized to only extract structures that will actually be used.

**Features**:
- **Selective extraction**: Only extracts required structures (from `required_structures.csv` or `required_structures_with_negatives.csv`)
- **Fallback capability**: Uses `cavity_mapping.csv` if required structures not available
- **Offline extraction** from local AlphaFold database archive
- **Parallel processing** with configurable thread count
- **Resume capability**: Skips already extracted files
- **Progress tracking**: Detailed logging and progress bars
- **Complete structures**: Includes all atoms needed for docking
- **Automatic integration**: Works seamlessly with existing workflow

**Configuration**: Edit `ALPHAFOLD_ARCHIVE` path and `NUM_THREADS` for optimal performance

**Usage**:
```bash
python extract_alphafold_models.py
```

**Output**:
- `alphafold_structures/` - Directory with AlphaFold structures organized by UniProt ID
- `alphafold_mapping.csv` - Updated mapping with AlphaFold file paths
- `alphafold_extraction.log` - Extraction log

**Requirements**:
- Local AlphaFold database archive
- `tarfile` and `gzip` Python libraries

### 6. fix_required_pdbs.py

Fixes PDB files by adding hydrogens and performing other corrections needed for accurate docking. This step only processes the structures identified as required, making it much more efficient than fixing all structures.

**Features**:
- **Selective processing** of only required structures
- **Hydrogen addition** at physiological pH (7.0)
- **Missing atom/residue fixing** using PDBFixer
- **Parallel processing** with configurable thread count
- **Resume capability** - skips already fixed structures

**Configuration**: Edit `NUM_THREADS` and `ALPHAFOLD_ARCHIVE` path

**Usage**:
```bash
python fix_required_pdbs.py
```

**Output**:
- `fixed_structures/` - Directory with fixed PDB files
- `fixed_mapping.csv` - Updated mapping with fixed structure paths
- `pdb_fixing.log` - Processing log

**Requirements**:
- `pdbfixer` and `openmm` packages
- Required structures list from step 3 or 4

### 7. convert_pdb_to_pdbqt.py
Converts receptor PDB files to PDBQT format with parallel processing. This step is optional but highly recommended as it significantly speeds up docking by pre-converting files. Uses fixed structures if available from step 6.

**Features**:
- **Parallel processing** with configurable thread count
- **Smart deduplication**: Same UniProt ID receptors are only converted once
- **Resume capability**: Skips already converted files
- **File verification**: Checks for identical structures using MD5 hashes
- **Progress tracking**: Detailed logging and progress bars

**Configuration**: Edit `NUM_THREADS` based on your CPU cores

**Usage**:
```bash
python convert_pdb_to_pdbqt.py
```

**Output**:
- `converted_pdbqt/` - Directory with PDBQT files organized by UniProt ID
- `pdbqt_mapping.csv` - Enhanced mapping with PDBQT file paths
- `pdb_to_pdbqt_conversion.log` - Conversion log

**Requirements**:
- OpenBabel Python bindings (`openbabel-python` or `conda install openbabel`)

### 8. batch_dock.py
Main docking execution script that uses pre-generated mappings to run consensus docking jobs for small molecule drugs only.

**Features**:
- **Small molecule filtering**: Only processes small molecule drugs from DrugBank
- **Automatic negative sample support**: Uses `required_structures_with_negatives.csv` if available, falls back to `DRUG_TO_PROTEIN_TSV` otherwise
- **Explicit cavity targeting**: When using negative samples file, docks only specified drug-target-pocket combinations (more efficient)
- Uses pre-extracted cavity mappings (faster startup)
- **Automatic PDBQT detection**: Uses pre-converted PDBQT files if available for faster docking
- **Smart mapping priority**: Uses fixed structures > AlphaFold structures > cavity structures
- Supports multiple cavities per protein
- Test mode for validation
- Resume capability (skips completed jobs)

**Docking Modes**:
1. **With negative samples** (`required_structures_with_negatives.csv` present):
   - Docks exact drug-target-pocket combinations specified in the file
   - Includes both positive (known interactions) and negative (non-interactions) samples
   - Only docks to specified cavities (more efficient, avoids unnecessary work)
   - Ideal for ML training datasets with balanced positive/negative samples

2. **Without negative samples** (fallback to `DRUG_TO_PROTEIN_TSV`):
   - Docks all known drug-protein interactions
   - Docks to all available cavities for each target
   - Original behavior for standard docking experiments

**Prerequisites**:
- `uniprot_gene_mapping.csv` (from step 1)
- `cavity_mapping.csv` (from step 2)
- `required_structures_with_negatives.csv` (from step 4, preferred) OR `DRUG_TO_PROTEIN_TSV` (fallback)
- `alphafold_mapping.csv` (from step 5, recommended for better results)
- `fixed_mapping.csv` (from step 6, optional but recommended)
- `pdbqt_mapping.csv` (from step 7, optional but recommended)

**Usage**:
```bash
python batch_dock.py
```

**Output**:
- `consensus_docking_results/` - Docking results for small molecule drugs only
- `docking_automation.log` - Docking log

## Key Improvements

### Multiple Cavity Support
All scripts now handle multiple cavities per protein:
- Detects `_vacant_1.pdb`, `_vacant_2.pdb`, etc.
- Detects `_cavity_1.pdb`, `_cavity_2.pdb`, etc.
- Creates separate docking jobs for each cavity

### AlphaFold Integration
- **Full structures**: Uses complete AlphaFold models instead of cavity-extracted "vacant" files
- **Complete atom sets**: AlphaFold models include all atoms needed for docking (no fixing required)
- **Better docking**: Complete protein structures provide better binding site context
- **Quality confidence**: Extracts confidence scores for structure quality assessment
- **Automated extraction**: Parallel extraction from local AlphaFold database archive

### Parallel PDBQT Conversion
- **Significant speed improvement**: Pre-converting PDB to PDBQT eliminates conversion time during docking
- **Intelligent deduplication**: Identical receptors (same UniProt ID) only converted once
- **Parallel processing**: Utilizes multiple CPU cores for faster conversion
- **Smart copying**: If multiple cavities use the same receptor, files are copied instead of re-converted

### Enhanced Performance
- **Faster docking startup**: Main script starts immediately with cached data
- **Optimized file handling**: PDBQT files are passed directly to consensus_docker.py
- **Better resource utilization**: Parallel conversion maximizes CPU usage

## Recommended Workflow

### Option 1: Manual Step-by-Step (Positive Samples Only)
```bash
# Step 1: Generate UniProt mappings (if not already done)
python prepare_uniprot_mapping.py

# Step 2: Extract cavities (can take a while, run once)
python extract_cavities.py

# Step 3: Identify required structures (filters for needed structures only)
python identify_required_structures.py

# Step 4: Extract AlphaFold structures (extracts only required structures)
python extract_alphafold_models.py

# Step 5: Fix PDB files (adds hydrogens to required structures)
python fix_required_pdbs.py

# Step 6: Convert PDB to PDBQT (saves significant time during docking)
python convert_pdb_to_pdbqt.py

# Step 7: Run batch docking (now much faster with pre-converted files)
python batch_dock.py
```

### Option 1b: Manual Step-by-Step (With Negative Samples for ML)
```bash
# Steps 1-3: Same as above
python prepare_uniprot_mapping.py
python extract_cavities.py
python identify_required_structures.py

# Step 4: Generate negative samples (optional, for ML training)
python generate_negative_samples.py

# Steps 5-8: Continue as normal (will process both positive and negative samples)
python extract_alphafold_models.py
python fix_required_pdbs.py
python convert_pdb_to_pdbqt.py
python batch_dock.py
```

### Option 2: Automated Workflow Runner
```bash
# Run the complete workflow automatically
python run_full_workflow.py

# Skip steps you've already completed
python run_full_workflow.py --skip-extract --skip-download

# Run in test mode
python run_full_workflow.py --test-mode
```

The workflow runner automatically:
- Runs each step in the correct order
- Checks for required files before proceeding  
- Handles errors gracefully (optional steps won't stop the workflow)
- Provides detailed logging and timing information
- Summarizes all generated files at the end

## Configuration

Before running the scripts, update the configuration paths in each script:

- **extract_cavities.py**: `CAVITY_TARBALL_FOLDER`
- **identify_required_structures.py**: `DRUG_TO_PROTEIN_TSV`, `SMALL_MOLECULE_DRUGS_CSV`, `PROCESSED_LIGAND_SDF_FOLDER`
- **generate_negative_samples.py**: `DRUG_TO_PROTEIN_TSV`, `SMALL_MOLECULE_DRUGS_CSV`, `CAVITY_MAPPING_CSV`, `PROCESSED_LIGAND_SDF_FOLDER`
- **extract_alphafold_models.py**: `NUM_THREADS` (adjust for extraction performance)
- **fix_required_pdbs.py**: `NUM_THREADS` (adjust for fixing performance)
- **convert_pdb_to_pdbqt.py**: `NUM_THREADS` (adjust for your CPU)
- **batch_dock.py**: `CONSENSUS_DOCKER_SCRIPT`, `PROCESSED_LIGAND_SDF_FOLDER`, `DRUG_TO_PROTEIN_TSV`, `SMALL_MOLECULE_DRUGS_CSV`

## File Organization

```
your_working_directory/
├── prepare_uniprot_mapping.py
├── extract_cavities.py
├── identify_required_structures.py
├── generate_negative_samples.py        # Optional: for negative sample generation
├── extract_alphafold_models.py
├── fix_required_pdbs.py
├── convert_pdb_to_pdbqt.py
├── batch_dock.py
├── run_full_workflow.py                # Automated workflow runner
├── uniprot_gene_mapping.csv            # Generated by step 1
├── cavity_mapping.csv                  # Generated by step 2
├── required_structures.csv             # Generated by step 3 (positive samples only)
├── required_structures_with_negatives.csv  # Generated by step 4 (optional, with negatives)
├── negative_samples_metadata.json      # Generated by step 4 (optional)
├── alphafold_mapping.csv               # Generated by step 5
├── fixed_mapping.csv                   # Generated by step 6
├── pdbqt_mapping.csv                   # Generated by step 7
├── extracted_cavities/                 # Generated by step 2
│   ├── AF-P12345-F1-model_v1/
│   │   ├── ...vacant_1.pdb
│   │   ├── ...cavity_1.pdb
│   │   ├── ...vacant_2.pdb
│   │   ├── ...cavity_2.pdb
│   └── ...
├── alphafold_structures/               # Generated by step 5
│   ├── P12345/
│   │   ├── AF-P12345-F1-model_v4.pdb
│   │   └── AF-P12345-F1-predicted_aligned_error_v4.json
│   └── ...
├── fixed_structures/                   # Generated by step 6
│   ├── P12345/
│   │   ├── AF-P12345-F1-model_v4_fixed.pdb
│   │   └── ...
│   └── ...
├── converted_pdbqt/                    # Generated by step 7
│   ├── P12345/
│   │   ├── AF-P12345-F1-model_v4.pdbqt
│   │   └── ...
│   └── ...
├── consensus_docking_results/          # Generated by step 8
│   ├── DB00001_GENE1_P12345_cavity_1/
│   ├── DB00001_GENE1_P12345_cavity_2/
│   └── ...
├── cavity_extraction.log
├── required_structures.log             # Generated by step 3
├── negative_samples.log                # Generated by step 4 (optional)
├── alphafold_extraction.log            # Generated by step 5
├── pdb_fixing.log                      # Generated by step 6
├── pdb_to_pdbqt_conversion.log         # Generated by step 7
└── docking_automation.log              # Generated by step 8
```

## Performance Benefits

### With AlphaFold + PDBQT Pre-conversion:
- **Better structures**: Complete AlphaFold models instead of cavity-extracted "vacant" files
- **Complete atom sets**: AlphaFold models include all necessary atoms (no fixing required)
- **Faster docking**: No PDB→PDBQT conversion during docking jobs
- **Reduced redundancy**: Each unique receptor converted only once
- **Better parallelization**: Conversion happens separately with optimal threading
- **Resume capability**: Can restart conversion where it left off
- **Quality confidence**: AlphaFold confidence scores help assess structure reliability

### Time Savings Example:
- **Without AlphaFold/PDBQT**: Each docking job includes ~30-60s conversion time + potential failures from malformed PDBs + poor binding site context
- **With AlphaFold + PDBQT**: Docking jobs start immediately, conversion done once upfront, complete structures, better results
- **For 1000 jobs**: Save ~8-16 hours of total processing time + improved docking quality

## Notes

- Steps 4-7 are optional but highly recommended for significant time savings and better docking quality
- Step 3 (identify required structures) is especially important as it optimizes the workflow to only process needed structures
- Step 4 (negative sample generation) is optional but essential for machine learning model training and evaluation
  - Negative samples allow proper assessment of model performance and prevent overfitting
  - Can be skipped if only running basic docking experiments
  - Downstream scripts automatically detect and use negatives if available
- Step 5 (AlphaFold extraction) replaces cavity-extracted "vacant" files with complete protein structures
- Run `extract_cavities.py` in test mode first to verify your cavity tarball folder is correct
- The structure identification, negative sample generation, AlphaFold extraction, PDB fixing, and PDBQT conversion steps only need to be run once unless you get new cavity data
- Use TEST_MODE=True in `batch_dock.py` to validate your setup before running full batch
- All scripts support resumption - they will skip work that's already been completed
- Ensure OpenBabel and requests are installed before running the respective scripts
