# Consensus## Workflow Overview

The docking workflow consists of seven main steps and processes **small molecule drugs only**:

1. **UniProt Mapping Generation** (`prepare_uniprot_mapping.py`)
2. **Cavity Extraction** (`extract_cavities.py`)
3. **Required Structures Identification** (`identify_required_structures.py`) - **NEW** (filters for small molecules)
4. **AlphaFold Model Extraction** (`extract_alphafold_models.py`) - **UPDATED** (extracts only required structures)
5. **PDB Fixing** (`fix_required_pdbs.py`) - **RECOMMENDED**
6. **PDB to PDBQT Conversion** (`convert_pdb_to_pdbqt.py`) - **OPTIONAL BUT RECOMMENDED**
7. **Batch Docking Execution** (`batch_dock.py`) - **UPDATED** (processes small molecules only)alysis Workflow

This repository contains scripts for running batch consensus doc**Prerequisites**:
- `uniprot_gene_mapping.csv` (from step 1)
- `cavity_mapping.csv` (from step 2)
- `cavityspace_mapping.csv` (from step 3, recom- **extract_alphafold_models.py**: `NUM_THREADS` (adjust for extraction performance)
- **identify_required_structures.py**: File paths for drug-protein data and ligands
- **fix_required_pdbs.py**: `NUM_THREADS` (adjust for fixing performance)
- **convert_pdb_to_pdbqt.py**: `NUM_THREADS` (adjust for your CPU cores)nded)
- `pdbqt_mapping.csv` (from step 4, optional but recommended)xperiments using AlphaFold structures and cavity data.

## Workflow Overview

The docking workflow consists of seven main steps:

1. **UniProt Mapping Generation** (`prepare_uniprot_mapping.py`)
2. **Cavity Extraction** (`extract_cavities.py`)
3. **Required Structures Identification** (`identify_required_structures.py`) - **NEW & RECOMMENDED**
4. **AlphaFold Model Extraction** (`extract_alphafold_models.py`) - **RECOMMENDED** (now extracts only required structures)
5. **PDB Fixing** (`fix_required_pdbs.py`) - **RECOMMENDED**
6. **PDB to PDBQT Conversion** (`convert_pdb_to_pdbqt.py`) - **OPTIONAL BUT RECOMMENDED**
7. **Batch Docking Execution** (`batch_dock.py`)

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

### 3. identify_required_structures.py ⭐ **NEW - OPTIMIZATION**
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

### 4. extract_alphafold_models.py ⭐ **UPDATED - RECOMMENDED**
Downloads full AlphaFold models for required structures only (if step 3 was run). This extracts structures directly from the local AlphaFold database archive, avoiding network limitations. Now optimized to only extract structures that will actually be used.

**Features**:
- **Selective extraction**: Only extracts required structures (from `required_structures.csv`)
- **Fallback capability**: Uses `cavity_mapping.csv` if required structures not available
- **Offline extraction** from local AlphaFold database
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

Identifies which AlphaFold structures are actually needed based on drug-protein interactions. This optimization step ensures we only process structures that will be used in docking, making the workflow more efficient.

**Features**:
- **Smart filtering** based on drug-protein interaction data
- **Ligand availability checking** to ensure corresponding SDF files exist
- **Priority ranking** by interaction count
- **Detailed reporting** of required vs available structures

**Configuration**: Edit paths for drug-protein data and ligand folder

**Usage**:
```bash
python identify_required_structures.py
```

**Output**:
- `required_structures.csv` - List of structures needed for docking
- `required_structures.log` - Processing log

**Requirements**:
- Drug-protein interaction TSV file
- UniProt mapping and cavity mapping files
- Processed ligand SDF folder

### 5. fix_required_pdbs.py ⭐ **NEW - RECOMMENDED**

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
- Required structures list from step 4

### 6. convert_pdb_to_pdbqt.py ⭐ **OPTIONAL BUT RECOMMENDED**
Converts receptor PDB files to PDBQT format with parallel processing. This step is optional but highly recommended as it significantly speeds up docking by pre-converting files. **Uses fixed structures if available from step 5.**

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

### 7. batch_dock.py ⭐ **UPDATED**
Main docking execution script that uses pre-generated mappings to run consensus docking jobs for **small molecule drugs only**.

**Features**:
- **Small molecule filtering**: Only processes small molecule drugs from DrugBank
- Uses pre-extracted cavity mappings (faster startup)
- **Automatic PDBQT detection**: Uses pre-converted PDBQT files if available for faster docking
- **Smart mapping priority**: Uses fixed structures > AlphaFold structures > cavity structures
- Supports multiple cavities per protein
- Test mode for validation
- Resume capability (skips completed jobs)

**Prerequisites**:
- `uniprot_gene_mapping.csv` (from step 1)
- `cavity_mapping.csv` (from step 2)
- `required_structures.csv` (from step 3, recommended for efficiency)
- `alphafold_mapping.csv` (from step 4, recommended for better results)
- `fixed_mapping.csv` (from step 5, optional but recommended)
- `pdbqt_mapping.csv` (from step 6, optional but recommended)

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
- **Quality confidence**: Downloads confidence scores for structure quality assessment
- **Automated download**: Parallel downloading from AlphaFold database

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

### Option 1: Manual Step-by-Step
```bash
# Step 1: Generate UniProt mappings (if not already done)
python prepare_uniprot_mapping.py

# Step 2: Extract cavities (can take a while, run once)
python extract_cavities.py

# Step 3: Extract AlphaFold structures (RECOMMENDED - better than vacant files)
python extract_alphafold_models.py

# Step 4: Convert PDB to PDBQT (RECOMMENDED - saves significant time later)
python convert_pdb_to_pdbqt.py

# Step 5: Run batch docking (now much faster with pre-converted files)
python batch_dock.py
```

### Option 2: Automated Workflow Runner ⭐ **NEW**
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
- **extract_alphafold_models.py**: `NUM_THREADS` (adjust for extraction performance)
- **fix_required_pdbs.py**: `NUM_THREADS` (adjust for fixing performance)
- **convert_pdb_to_pdbqt.py**: `NUM_THREADS` (adjust for your CPU)
- **batch_dock.py**: `CONSENSUS_DOCKER_SCRIPT`, `PROCESSED_LIGAND_SDF_FOLDER`, `DRUG_TO_PROTEIN_TSV`, `SMALL_MOLECULE_DRUGS_CSV`

## File Organization

```
your_working_directory/
├── prepare_uniprot_mapping.py
├── extract_cavities.py
├── extract_alphafold_models.py         # NEW
├── identify_required_structures.py     # NEW
├── fix_required_pdbs.py               # NEW
├── convert_pdb_to_pdbqt.py             # NEW
├── batch_dock.py
├── run_full_workflow.py                # NEW - Automated workflow runner
├── uniprot_gene_mapping.csv            # Generated by step 1
├── cavity_mapping.csv                  # Generated by step 2
├── alphafold_mapping.csv               # Generated by step 3 (recommended)
├── required_structures.csv             # Generated by step 4 (optimization)
├── fixed_mapping.csv                   # Generated by step 5 (recommended)
├── pdbqt_mapping.csv                   # Generated by step 6 (optional)
├── extracted_cavities/               # Generated by step 2
│   ├── AF-P12345-F1-model_v1/
│   │   ├── ...vacant_1.pdb
│   │   ├── ...cavity_1.pdb
│   │   ├── ...vacant_2.pdb
│   │   └── ...cavity_2.pdb
│   └── ...
├── cavityspace_structures/             # Generated by step 3 (recommended)
│   ├── P12345/
│   │   ├── AF-P12345-F1-model_v4.pdb
│   │   └── AF-P12345-F1-predicted_aligned_error_v4.json
│   └── ...
├── converted_pdbqt/                    # Generated by step 4 (optional)
│   ├── P12345/
│   │   ├── AF-P12345-F1-model_v4.pdbqt
│   │   └── ...
│   └── ...
├── consensus_docking_results/          # Generated by step 5
│   ├── DB00001_GENE1_P12345_cavity_1/
│   ├── DB00001_GENE1_P12345_cavity_2/
│   └── ...
├── cavity_extraction.log
├── cavityspace_download.log            # NEW
├── pdb_to_pdbqt_conversion.log         # NEW
└── docking_automation.log
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

- **Steps 3-4 are optional but highly recommended** for significant time savings and better docking quality
- **Step 3 (AlphaFold download) is especially important** as it replaces cavity-extracted "vacant" files with complete protein structures
- Run `extract_cavities.py` in test mode first to verify your cavity tarball folder is correct
- The AlphaFold download and PDBQT conversion steps only need to be run once unless you get new cavity data
- Use TEST_MODE=True in `batch_dock.py` to validate your setup before running full batch
- All scripts support resumption - they will skip work that's already been completed
- Ensure OpenBabel and requests are installed before running the respective scripts
