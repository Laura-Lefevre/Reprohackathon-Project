###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

# Chargement du fichier de configuration
configfile: "config.yaml"

THREADS = config["THREADS"]

SRR_LIST=["SRR10510434", "SRR10510435", "SRR10510436", "SRR10510437", "SRR10510438", "SRR10510439"]
OUTDIR = "output/fastq"

#include: "rules/downloading_fastq.rules"


###############################################
# RULE 1 — Download FASTQ
###############################################

rule download_fastq:
    output:
        "output/fastq/{srr}.fastq"
    threads: THREADS
    shell:
        """
        echo "Downloading {wildcards.srr}..."
        mkdir -p {OUTDIR}

        if [ -f {output} ]; then
            echo "Existing file detected: deleting {output}"
            rm {output}
        fi

        fasterq-dump --threads {threads} --progress --outdir {OUTDIR} {wildcards.srr}

        TMPDIR=$(find {OUTDIR} -maxdepth 1 -type d -name "fasterq.tmp.*" | head -n 1)
        if [ -n "$TMPDIR" ]; then
            echo "Deleting temporary folder: $TMPDIR"
            rm -rf "$TMPDIR"
        fi

        echo "Download complete for {wildcards.srr}"
        """

rule download_all_fastq:
    input:
        expand(f"{OUTDIR}/{{srr}}.fastq", srr=SRR_LIST)

# Pour lancer jusqu'à cette partie : snakemake --cores 4 download_all_fastq

###############################################
# RULE 2 — Trimming
###############################################

rule trimming:
    input:
        fastq="output/fastq/{srr}.fastq"
    output:
        trimmed_fastq="output/trimmed/{srr}_trimmed.fq"
    threads: THREADS
    params:
        outdir="output/trimmed",
        quality=20,       # par exemple
        phred=33,         # 33 ou 64 selon ton jeu de données
        minlen=25         # longueur minimale des reads
    shell:
        """
        mkdir -p {params.outdir}

        echo "Trimming {input.fastq} with cutadapt..."
        
        cutadapt -m 25 -q {params.quality} -o {output.trimmed_fastq} {input.fastq}
        
        echo "Trimming complete for {wildcards.srr}"
        """

rule trimming_all:
    input:
        expand("output/trimmed/{srr}_trimmed.fq", srr=SRR_LIST)


# Pour lancer jusqu'à cette partie : snakemake --cores 4 trimming_all


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



