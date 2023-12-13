"""Lemcke et al., 2023"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("Run", drop=False)
resolutions = [0.001]
bprj = "PRJNA938499"
prj  = "lemcke2023-lha-infl"

def plots_doublets_raw(wildcards):
    x = "output/figures/{wildcards.run}_raw/doublets_call".format(wildcards=wildcards)
    return x.replace("\.", "_")


def get_mem_mb(wildcards, attempt):
    return attempt * 500000


##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand("cellbender/{run}/{run}_output_filtered.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("scrublet/{run}/{run}_initial_annotation.h5ad",
                run=samples["Run"]),
        # expand(["output/figures/combined-top5_logreg-umap-whole_dataset-fpr_{res}.pdf",
        #     "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_{res}.pdf",
        #     "output/tables/01A-eda-whole_dataset-fpr_{res}/parameters.json"], res=resolutions),
        # expand(["data/{bprj}-whole_dataset-fpr_{res}-clusters.h5Seurat",
        #     "data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad",
        #     "data/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.h5ad"], bprj=bprj, res=resolutions),
        # expand(["output/tables/01A-eda-whole_dataset-fpr_{res}/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_{res}.csv",
        #     "output/tables/01A-eda-whole_dataset-fpr_{res}/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_{res}.csv"], prj=prj, res=resolutions),
        # expand(["data/{bprj}-whole_dataset-nc-clusters.h5Seurat",
        #     "data/{bprj}-whole_dataset-nc-clusters.h5ad",
        #     "data/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.h5ad"], bprj=bprj),
        # expand(["output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
        #     "output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-logreg_sct-combined-whole_dataset-nc.csv"], prj=prj),
        # "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
        # "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
        # "output/tables/01-eda-whole_dataset-nc/parameters.json",

##### load rules #####

CELLRANGER="cd cellranger && source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    input:
        sample=directory("fastq"),
        idx=directory("/home/etretiakov/src/scRNA-seq-references/mm10_optimized")
    output:
        raw="cellranger/{run}/outs/raw_feature_bc_matrix.h5",
        filtered="cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
        summary="cellranger/{run}/outs/web_summary.html",
        bam="cellranger/{run}/outs/possorted_genome_bam.bam",
    params:
        ids="cellranger/{run}",
        sample="{run}"
    threads: 32
    resources:
        mem_mb=64000
    shell:
        ("{CELLRANGER} count --include-introns true \
            --id={params.sample} \
            --sample={params.sample} \
            --transcriptome={input.idx} \
            --fastqs={input.sample} \
            --jobmode=local \
            --localcores={threads} ")

rule cellbender:
    input:
        "cellranger/{run}/outs/raw_feature_bc_matrix.h5"
    output:
        expand(["cellbender/{{run}}/{{run}}_output.h5", "cellbender/{{run}}/{{run}}_output_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender/{run}/{run}_output.h5"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    threads: 4
    resources:
        nvidia_gpu=1,
        mem_mb=10000
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {params.h5} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.001 \
            --epochs 150")

rule doublets_call:
    input:
        filt_h5="cellbender/{run}/{run}_output_filtered.h5"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls.tsv",
        dr="cellbender/{run}/{run}_latent_gene_expression.csv",
        h5ad="scrublet/{run}/{run}_initial_annotation.h5ad"
    params:
        expected_dblt=lambda wildcards: samples["NExpectedDoubletRate"][wildcards.run],
        sample_run_name="{run}",
        plots=plots_doublets_raw
    container:
        "docker://etretiakov/scrna-seq-cellbender:2023-10-08_2"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "code/scrublet_cb-z.py"





# rule exploratory_data_analysis_0_001:
#     input:
#         rmd="analysis/01A-eda-whole_dataset-fpr_0.001.Rmd"
#     output:
#         "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.001.pdf",
#         "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.001.pdf",
#         f"output/tables/01A-eda-whole_dataset-fpr_0.001/{prj}_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.001.csv",
#         f"output/tables/01A-eda-whole_dataset-fpr_0.001/{prj}_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.001.csv",
#         "output/tables/01A-eda-whole_dataset-fpr_0.001/parameters.json",
#         f"data/{bprj}-whole_dataset-fpr_0.001-clusters.h5Seurat"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 32
#     resources:
#         mem_mb=get_mem_mb
#     shell:
#         ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


# rule exploratory_data_analysis:
#     input:
#         rmd="analysis/01-eda-whole_dataset-nc.Rmd",
#         raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
#         scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.001.tsv", run=samples["Run"])
#     output:
#         "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
#         "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
#         f"output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
#         f"output/tables/01-eda-whole_dataset-nc/{prj}_all_mrk-logreg_sct-combined-whole_dataset-nc.csv",
#         "output/tables/01-eda-whole_dataset-nc/parameters.json",
#         f"data/{bprj}-whole_dataset-nc-clusters.h5Seurat"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 32
#     resources:
#         mem_mb=get_mem_mb
#     shell:
#         ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


# rule convert_seurat_to_h5ad:
#     input:
#         h5ad="data/{bprj}-whole_dataset-fpr_{res}-clusters.h5Seurat"
#     output:
#         h5ad="data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 2
#     resources:
#         mem_mb=30000
#     script:
#         "../code/convert_h5ad.R"

# rule convert_nc_seurat_to_h5ad:
#     input:
#         h5ad=f"data/{bprj}-whole_dataset-nc-clusters.h5Seurat"
#     output:
#         h5ad=f"data/{bprj}-whole_dataset-nc-clusters.h5ad"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 2
#     resources:
#         mem_mb=30000
#     script:
#         "../code/convert_h5ad.R"


# rule subset_astrocytes:
#     input:
#         "data/{bprj}-whole_dataset-fpr_{res}-clusters.h5ad"
#     output:
#         h5ad_annotations_all="data/class_cello/{bprj}-whole_dataset-{res}-cello_annotation.h5ad",
#         tables_annotations_all="output/tables/class_cello/{bprj}-whole_dataset-{res}-CellO_output.tsv",
#         h5ad_annotations_astrocytes="data/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.h5ad",
#         tables_annotations_astrocytes="output/tables/class_cello/{bprj}-astrocytes_dataset-{res}-initial_selection.tsv"
#     params:
#         prj=prj,
#         bioprj=bprj,
#         res="{res}"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 2
#     resources:
#         mem_mb=30000
#     script:
#         "../code/class_cello.py"


# rule subset_nc_astrocytes:
#     input:
#         "data/{bprj}-whole_dataset-nc-clusters.h5ad"
#     output:
#         h5ad_annotations_all="data/class_cello/{bprj}-whole_dataset-nc-cello_annotation.h5ad",
#         tables_annotations_all="output/tables/class_cello/{bprj}-whole_dataset-nc-CellO_output.tsv",
#         h5ad_annotations_astrocytes="data/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.h5ad",
#         tables_annotations_astrocytes="output/tables/class_cello/{bprj}-astrocytes_dataset-nc-initial_selection.tsv"
#     params:
#         prj=prj,
#         bioprj=bprj,
#         res="nc"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 2
#     resources:
#         mem_mb=30000
#     script:
#         "../code/class_cello.py"
