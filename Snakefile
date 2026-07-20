configfile: "config.yaml"

ANNOTATION_DATA = config["data"]["annotation_data"]
NEURONS = config["data"]["neurons_subclass"]

rule all:
    input:
        "results/composition",
        "results/neurons",
        "results/ppi",
        "results/neuronchat",
        "results/C0"


rule composition_tsne:
    input:
        annotation_data=ANNOTATION_DATA,
        neurons=NEURONS
    output:
        directory("results/composition")
    conda:
        "envs/scrna.yaml"
    script:
        "scripts/01_composition_tsne.R"


rule neuronal_deg_enrichment:
    input:
        neurons=NEURONS
    output:
        directory("results/neurons")
    conda:
        "envs/scrna.yaml"
    script:
        "scripts/02_neuronal_deg_enrichment.R"


rule ppi:
    input:
        neuronal_results=rules.neuronal_deg_enrichment.output
    output:
        directory("results/ppi")
    conda:
        "envs/scrna.yaml"
    script:
        "scripts/03_ppi_analysis.R"


rule neuronchat:
    input:
        neurons=NEURONS
    output:
        directory("results/neuronchat")
    conda:
        "envs/neuronchat.yaml"
    script:
        "scripts/04_neuronchat_analysis.R"


rule C0_deg_enrichment:
    input:
        neurons=NEURONS
    output:
        directory("results/C0")
    conda:
        "envs/scrna.yaml"
    script:
        "scripts/05_C0_deg_enrichment.R"
