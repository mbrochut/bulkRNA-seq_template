configfile: "config.yaml"
import yaml

def get_contrast_outputs(config):
    contrasts = config["03_create_contrast"]["paired_contrast"]
    return [
        f"results/contrasts/condition_{c2}_VS_{c1}.csv"
        for c1, c2 in contrasts
    ]


def get_pathway_outputs(config):
    contrasts = config["03_create_contrast"]["paired_contrast"]

    # Build contrast names
    contrast_names = [
        f"{c2}_VS_{c1}"
        for c1, c2 in contrasts
    ]

    gsea_dbs = ["GO", "KEGG", "Hallmark"]
    ora_types = ["ALL", "UP", "DOWN"]

    outputs = []

    for contrast in contrast_names:
        # GSEA files
        for db in gsea_dbs:
            outputs.append(
                f"results/pathway/GSEA/GSEA_{db}_condition_{contrast}.csv"
            )

        # ORA files (only GO)
        for t in ora_types:
            outputs.append(
                f"results/pathway/ORA/ORA_GO_condition_{contrast}_{t}.csv"
            )

    return outputs

CONTRAST_OUTPUTS = get_contrast_outputs(config)
PATHWAY_OUTPUTS = get_pathway_outputs(config)
print(PATHWAY_OUTPUTS)

rule all:
    input:
        "results/papermill/02_data_exploration.ipynb",
        CONTRAST_OUTPUTS,
        PATHWAY_OUTPUTS,
        "QUARTO/_quarto.yml"

rule init_repo:
    output:
        touch("results/.init_done")
    params:
        organism=config["general"]["organism"]
    shell:
        """
        python 00_init_repo.py --organism {params.organism}
        touch {output}
        """

rule create_anndata:
    input:
        nb="01_create_anndata_object.ipynb",
        init="results/.init_done"
    output:
        nb="results/papermill/01_create_anndata_object.ipynb",
        adata=config["paths"]["adata_filtered"]
    params:
        yaml_params=lambda wc: yaml.dump({
            "path_to_your_data": config["paths"]["counts"],
            "path_to_your_metadata": config["paths"]["metadata"],
            "gene_id_col": config["01_create_anndata_object"]["gene_id_col"],
            "condition_columns": config["01_create_anndata_object"]["condition_columns"],
            "filtering_sum": config["01_create_anndata_object"]["filtering_sum"],
            "adata_output": config["paths"]["adata_filtered"],
        }),
        adata_out=config["paths"]["adata_filtered"]
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml_params}'
        """



rule data_exploration:
    input:
        nb="02_data_exploration.ipynb",
        adata=config["paths"]["adata_filtered"]
    output:
        nb="results/papermill/02_data_exploration.ipynb"
    params:
        yaml_params=lambda wc: yaml.dump({
            "path_adata_filter": config["paths"]["adata_filtered"],
            "condition": config["general"]["condition"],
            "design": config["02_data_exploration"]["design"],
        })
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml_params}'
        """


rule create_contrasts:
    input:
        nb="03_create_contrast.ipynb",
        adata=config["paths"]["adata_filtered"]
    output:
        nb="results/papermill/03_create_contrast.ipynb",
        csvs=CONTRAST_OUTPUTS
    params:
        yaml_params=lambda wc: yaml.dump({
            "paired_contrast": config["03_create_contrast"]["paired_contrast"],
            "condition": config["general"]["condition"],
            "design": config["03_create_contrast"]["design"],
            "path_adata_filter": config["paths"]["adata_filtered"],
        })
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml_params}'
        """


rule compute_pathway:
    input:
        nb="04_compute_pathway_per_contrast.ipynb",
        contrasts=CONTRAST_OUTPUTS
    output:
        nb="results/papermill/04_compute_pathway_per_contrast.ipynb",
        results=PATHWAY_OUTPUTS
    params:
        yaml_params=lambda wc: yaml.dump({
            "organism": config["general"]["organism"]
        })
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y '{params.yaml_params}'
        """


rule generate_quarto:
    input:
        PATHWAY_OUTPUTS
    output:
        "QUARTO/_quarto.yml"
    shell:
        """
        python 05_generate_quarto_files.py
        """