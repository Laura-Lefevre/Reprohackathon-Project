###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

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



rule all:
    input:
        "output/counts/all_samples_counts.txt"


###############################################
# 0 — Download images from Zenodo
###############################################

###############################################
# 1 — Download FASTQ
###############################################

#rule download_all_fastq:
#    input:
#        expand(f"{OUTDIR}/{{srr}}.fastq", srr=SRR_LIST)

###############################################
# 2 — Trimming
###############################################

#rule trimming_all:
#    input:
#        expand("output/trimmed/{srr}_trimmed.fq", srr=SRR_LIST)

###############################################
# 3 — Mapping
###############################################

#rule mapping_all:
#    input:
#        expand("output/mapping/{srr}.bam", srr=SRR_LIST),
#        expand("output/mapping/{srr}.bam.bai", srr=SRR_LIST)

###############################################
# 4 — Counting
###############################################

#rule counting_all:
#    input:
#        "output/counts/all_samples_counts.txt"

###############################################
# 5 — Statistical Analysis (DESeq2)
###############################################
