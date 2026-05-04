configfile: "config.yaml"

def get_contrast_outputs(config):
    contrasts = config["03_create_contrast"]["paired_contrast"]
    return [
        f"results/contrasts/condition_{c1}_VS_{c2}.csv"
        for c1, c2 in contrasts
    ]

CONTRAST_OUTPUTS = get_contrast_outputs(config)

rule all:
    input:
        "results/papermill/02_data_exploration.ipynb",
        CONTRAST_OUTPUTS

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
        counts=config["paths"]["counts"],
        metadata=config["paths"]["metadata"],
        gene_id=config["01_create_anndata_object"]["gene_id_col"],
        condition_columns=config["01_create_anndata_object"]["condition_columns"],
        filtering_sum=config["01_create_anndata_object"]["filtering_sum"],
        adata_out=config["paths"]["adata_filtered"]
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p path_to_your_data "{params.counts}" \
            -p path_to_your_metadata "{params.metadata}" \
            -p gene_id_col "{params.gene_id}" \
            -p condition_columns '{params.condition_columns}' \
            -p filtering_sum {params.filtering_sum} \
            -p adata_output "{params.adata_out}"
        """


rule data_exploration:
    input:
        nb="02_data_exploration.ipynb",
        adata=config["paths"]["adata_filtered"]
    output:
        nb="results/papermill/02_data_exploration.ipynb"
    params:
        condition=config["general"]["condition"],
        design=config["02_data_exploration"]["design"]
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p path_adata_filter "{input.adata}" \
            -p condition "{params.condition}" \
            -p design "{params.design}"
        """


rule create_contrasts:
    input:
        nb="03_create_contrast.ipynb",
        adata=config["paths"]["adata_filtered"]
    output:
        nb="results/papermill/03_create_contrast.ipynb",
        csvs=CONTRAST_OUTPUTS
    params:
        paired_contrast=config["03_create_contrast"]["paired_contrast"],
        condition=config["general"]["condition"],
        design=config["03_create_contrast"]["design"]
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p paired_contrast '{params.paired_contrast}' \
            -p condition "{params.condition}" \
            -p design "{params.design}" \
            -p path_adata_filter "{input.adata}"
        """