# RNA-seq Differential Expression & Pathway Analysis Pipeline

This repository provides a full pipeline to go from **raw count data** to:

* Differential expression analysis (DESeq2 via PyDESeq2)
* Pathway enrichment (GSEA / ORA)
* Automated Quarto reports

---
## Setup

### 1. Clone the repository

```bash
git clone https://github.com/mbrochut/bulkRNA-seq_template.git
cd bulkRNA-seq_template
```

---

### 2. Create and activate environment

```bash
python3 -m venv env

# bash
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
python 00_init_repo.py --organism Mouse # TO USE WITH TESTING DATA
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

## Quarto Report Generation (`05_generate_quarto_files.py`)


### Install Quarto

Not sure if necessary
https://quarto.org/docs/download/

Quarto report generation is now fully driven by the `config.yaml` file.  
This step is mandatory and no longer relies on hardcoded parameters.

### Configuration (`config.yaml`)

The section `05_generate_quarto` controls how the report is built:

```yaml
05_generate_quarto:
  contrast_folder: "results/contrasts/"
  output_folder: "QUARTO/"
  template_dir: "Quarto_template/"

  project_title: "YOUR TITLE"
  author_name: "AUTHOR NAME"
  split_to_remove: 1 # for now this parameter is used to isolate the contrast name in generation of files.

  modules:
    QC:
      name: "Quality Control"
      type: "single"
      file: "QC.qmd"
      include: true

    DE:
      name: "Differential Analysis"
      type: "menu"
      template: "Differential_analysis_template.qmd"
      prefix: ""
      model: "condition"
      include: true
```

### Modules system
By default, all the modules are include.
Modules define what appears in the final Quarto report.

Each module has:

| Parameter | Description |
|----------|------------|
| `name` | Display name in the navbar |
| `type` | `"single"` or `"menu"` |
| `include` | Enable or disable the module |
| `file` | (single) static `.qmd` file |
| `template` | (menu) template used to generate multiple pages |
| `prefix` | Prefix added to generated files |
| `model` | Column used to group contrasts |

### Module types

#### Single

- One static `.qmd` file
- Used for global sections (QC, Venn, summary)

#### Menu

- Generates one page per contrast
- Uses a template located in `Quarto_template/`
- 

The script automatically:
- loops over contrasts
- fills the template
- creates one `.qmd` per contrast
- adds them to the Quarto navbar


### Run Quarto generation

```
python 05_generate_quarto_files.py
```

### Generate report

```bash
cd QUARTO
quarto render
```

---

# Snakemake Pipeline

A Snakemake pipeline is provided to automate the full workflow.

to lunch the pipeline:
```
snakemake --cores N # use N as the number of core you want to use
```


## Required configuration

The pipeline relies on `config.yaml`.  
You must define at least:

```yaml
general:
  organism: Mouse # CHOSE YOUR ORGANISM: Human or Mouse

paths:
  counts: "data/salmon.merged.gene_counts.tsv"
  metadata: "meta/meta.xlsx"

01_create_anndata_object:
  gene_id_col: "gene_id" # Colname of gene id (gene_id from the NF-core pieline)
  condition_columns: # Define which metadata columns are combined to create the `condition` column 
    - "treatment"
    - "Age"

03_create_contrast:
  paired_contrast:
    - ["Control_Young","Trained_Young"] # always: [reference, test] in that order. example: ["control", "treated"]
```


## Key parameters

| Parameter | Description |
|----------|------------|
| `organism` | Used for pathway databases (GSEA / ORA) |
| `counts` | Gene count matrix |
| `metadata` | Sample metadata |
| `gene_id_col` | Column containing gene identifiers |
| `condition_columns` | Columns combined to build the `condition` variable |
| `paired_contrast` | List of comparisons `[reference, test]` |



## Pipeline steps

The pipeline executes the following steps:

1. Initialization (folder structure and database download)
2. AnnData creation
3. Data exploration
4. Differential expression (DESeq2)
5. Pathway analysis (GSEA / ORA)
6. Quarto file generation
7. Quarto rendering



## Initialization behavior

The pipeline automatically handles initialization:

- If the project is not initialized, it will create the structure and download databases
- If already initialized, the step is skipped
- If databases already exist, they are reused

If you change the `organism`, make sure the corresponding databases exist or rerun initialization.



## Run full pipeline

```
snakemake --cores 8
```



## Run a specific step

```
snakemake quarto_render
```


## Notes

- Notebooks can still be run independently without Snakemake
- The YAML file ensures reproducibility and consistency across the pipeline
- The module system allows easy extension (e.g. TF inference or additional analyses)


## 🎉 Summary

This pipeline provides:

* Clean DESeq2 workflow
* Automated pathway analysis
* Fully reproducible reporting

👉 From raw counts → publication-ready report in a few steps


## In Progress

The following features are currently under development:

* **Heatmaps**
  → Improved visualization of gene expression and pathway activity for specific genes of interest

---

Contributions, suggestions, and feedback are welcome!
