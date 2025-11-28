######################################################
### Reprohackathon Project - 2025                  ###
### Authors : Saïda MOUSSAEVA,                     ###
### Camille LE CORRE, Laura LEFEVRE, Muyao GUO     ###
### Description : Snakefile with rule all to       ###
### run the pipeline                               ###
######################################################

# Loading configuration file
configfile: "config.yaml"

# Import of variables from config file
THREADS = config["THREADS"]
SRR_LIST = config["SRR_LIST"]
OUTDIR = "output/fastq"
IMAGEDIR = "output/images"
SIF_LIST = config["SIF_LIST"]
ZENODO_URL = config["ZENODO_URL"]

# Including rules from separate files for modularity
include: "rules/downloading_fastQ.rules"
include: "rules/trimming.rules"
include: "rules/downloading_ref_genome.rules"
include: "rules/downloading_ref_annot.rules"
include: "rules/bowtie_build_index.rules"
include: "rules/downloading_gtf.rules"
include: "rules/mapping.rules"
include: "rules/counting.rules"
include: "rules/deseq2.rules"
include: "rules/enrichment.rules"
include: "rules/plots.rules"

# Main rule
rule all:
    input:
        "output/counts/all_samples_counts.txt",

        "output/deseq2/res_df.csv",
        "output/deseq2/deg.csv",
        "output/deseq2/deg_up.csv",
        "output/deseq2/deg_down.csv",
        "output/vsd.rds",

        "output/enrichment/gene_ids.csv",
        "output/enrichment/gene_pathways.csv",
        "output/enrichment/annot_long.csv",
        "output/enrichment/df_enrich.csv",
        "output/enrichment/df_enrich_top25.csv",

        "output/plots/plots.done"