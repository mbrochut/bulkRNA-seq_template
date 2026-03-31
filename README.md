# RNA-seq Differential Expression & Pathway Analysis Pipeline

This repository provides a full pipeline to go from **raw count data** to:

* Differential expression analysis (DESeq2 via PyDESeq2)
* Pathway enrichment (GSEA / ORA)
* Automated Quarto reports

---
## Setup

### 1. Clone the repository

```bash
git clone <your-repo-url>
cd <your-repo-name>
```

---

### 2. Create and activate environment

```bash
python -m venv env

# bash / zsh
source env/bin/activate

# fish
source env/bin/activate.fish
```

---

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

---

### 4. Enable Plotly image export

```bash
plotly_get_chrome
```

---

## Project Initialization

```bash
python 00_init_repo.py --organism Mouse
# or
python 00_init_repo.py --organism Human
```

This will:

* Create the full project structure
* Download pathway databases for the selected organism into:

  ```
  data/DB/<organism>/
  ```

### Project Structure

```
data/
  DB/<organism>/
meta/
QUARTO/
results/
  contrasts/
  models/
  pathway/
    concat/
    GSEA/
    ORA/
    GSEA_object/
  QC/
```

## Input Data Format (IMPORTANT ⚠️)

### Expression matrix (`.tsv`)

Must be **wide format**:

| gene_id    | gene_name | OC1  | OC2  | ... |
| ---------- | --------- | -------- | -------- | --- |
| ENSMUSG... | Gnai3     | 6591     | 6228     | ... |

* Must include:

  * `gene_id`
  * `gene_name`

### Metadata (`.xlsx`)

| id | treatment | stimulation | Age |
| --------- | --------- | ----------- | --- |
| OC1       | Control   | NO          | Old |

* `id` must match column names in expression matrix


---

## Pipeline Steps

Run notebooks **in order**:

---

### **01_create_anndata_object.ipynb**

Creates the main object:

```
adata_filter_genes.h5ad
```

#### Parameters:

```python
path_to_your_data = "./data/salmon.merged.gene_counts.tsv"
path_to_your_metadata = "./meta/meta.xlsx"

gene_id_col = 'gene_id'

condition_columns = ['treatment', 'Age']
# → merged into a single "condition" column: treatment_Age. You can put only one column if needed

filtering_sum = 10
# → keeps genes with total counts ≥ 10
```

Notes:

* Additional filtering uses default PyDESeq2 settings
* The `condition` column is used for all downstream analysis

---

### **02_data_exploration.ipynb**

* Quality control plots
* PCA / distributions
* General dataset exploration

Free exploration — no strict requirements

---

### **03_create_contrast.ipynb**

Define comparisons of interest:

```python
paired_contrast = [
    ("Control_Young", "Trained_Young"),
    ("Control_Young", "Control_Old"),
]
```

#### Model design:

```python
design = '~condition'
```

Outputs:

* Contrast files saved in `results/contrasts/`

---

### **04_compute_pathway_per_contrast.ipynb**

Runs pathway enrichment for each contrast:

* GSEA
* ORA

#### Parameter:

```python
organism = "Mouse"  # or "Human"
```

Uses databases from:

```
data/DB/<organism>/
```

Outputs saved in:

```
results/pathway/
```

---

### **05_generate_quarto_contrast_file.py**

```bash
python 05_generate_quarto_contrast_file.py
```

This script:

* Generates one `.qmd` per contrast from a template
* Copies:

  * `QC.qmd`
  * `utils.py`
* Creates `_quarto.yml` automatically

---

## Quarto Report

### Install Quarto

Not sure if necessary
https://quarto.org/docs/download/

### Generate report

```bash
cd QUARTO
quarto render
```

Outputs a full website with:

* QC page
* One page per contrast
* Interactive plots

---

## 🎉 Summary

This pipeline provides:

* Clean DESeq2 workflow
* Automated pathway analysis
* Fully reproducible reporting

👉 From raw counts → publication-ready report in a few steps
