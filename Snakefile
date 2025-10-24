###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

# Chargement du fichier de configuration
configfile: "config.yaml"


###############################################
# RULE 1 — Download FASTQ
###############################################




###############################################
# RULE 2 — Trimming
###############################################




###############################################
# RULE 3 — Mapping
###############################################




###############################################
# RULE 4 — Counting
###############################################




###############################################
# RULE 5 — Statistical Analysis (DESeq2)
###############################################
rule deseq2_test:
    input:
        "output/counts.tsv",
    output:
        "results/deseq2_results.csv"
    container:
        "containers/DESeq2.sif"
    shell:
        "Rscript scripts/DESeq2R"



