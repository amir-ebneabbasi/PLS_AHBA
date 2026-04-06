# PLS Pipeline for Imaging Transcriptomics


## 📦 Data Requirements

### AHBA Gene Expression
- Rows: Brain regions
- Columns: Genes

### Phenotype Data
- Rows: Brain regions
- Must include: phenotype column (`pheno_name`)

### Spins (Spatial Null Model)
- Spatial permutation indices
- Each column = one permutation

### Pathways
- Gene sets in GMX (tab-separated)


## 🔧 Pipeline Overview

### AHBA Processing (abagen)
Generates regional gene expression matrices using the *abagen* toolbox.

- **Script:** `abagen.sh`
- **Output:** ROI × gene expression matrix (input to PLS)

### Spin Permutations (Spatial Null Model)
Generates spatially constrained permutations to control for spatial autocorrelation.

- **Script:** `generate_spins.R`
- **Inputs:** ROI coordinates
- **Methods:** `"hungarian"`, `"vasa"`
- **Output:** ROI × nrot permutation index matrix

## 🧠 PLS + FGSEA Pipeline
Links AHBA gene expression to brain phenotypes using PLS and explores pathway enrichment.
### Step 1 — Data Preprocessing
- Removes genes with excessive missingness
- Ensures ROI alignment across datasets
- Merges phenotype + gene expression

### Step 2 — PLS Regression
- Extracts explained variance (R²)
- Aligns component signs for interpretability

### Step 3 — Spin Permutation Test
- Applies spatial permutations
- Computes spin p-values for components

### Step 4 — Bootstrapping
- Estimates variability of gene weights
- Computes corrected weights (Z-like scores)

### Step 5 — FGSEA
- Uses PLS1 component
- Performs preranked gene set enrichment analysis (GSEA)


## 📤 Outputs

### Results_PLS_Spin.csv
- Variance explained (R²)
- Spin-test p-values

### Results_PLS_cWeights.csv
- Bootstrapped gene weights
- Used for biological interpretation

### Results_FGSEA.csv
- Enriched pathways
- Normalized Enrichment Score (NES)
- Leading-edge genes
