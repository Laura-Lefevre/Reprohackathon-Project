###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

# Chargement du fichier de configuration
configfile: "config.yaml"

THREADS = config["THREADS"]

SRR_LIST=["SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726"]

OUTDIR = "output/fastq"
IMAGEDIR = "output/images"

#include: "rules/downloading_fastq.rules"


###############################################
# RULE 1 — Download FASTQ
###############################################

rule build_sif_download_fastq:
    input:
        definition = config["def_download_fastq"]
    output:
        image = f"{IMAGEDIR}/downloading_fastq.sif"
    params:
        build_cmd = config.get("build_cmd", "apptainer build")
    shell:
        """
        mkdir -p {IMAGEDIR}
        {params.build_cmd} {output.image} {input.definition}
        """

rule download_fastq:
    input:
        sif="output/images/downloading_fastq.sif"
    output:
        "output/fastq/{srr}.fastq"
    threads: THREADS
    container:
        "output/images/downloading_fastq.sif"
    shell:
        """
        echo "Downloading {wildcards.srr}..."
        mkdir -p output/fastq

        if [ -f output/fastq/{wildcards.srr}.fastq ]; then
            rm output/fastq/{wildcards.srr}.fastq
        fi

        fasterq-dump --threads 4 --progress --outdir output/fastq {wildcards.srr}

        TMPDIR=$(find output/fastq -maxdepth 1 -type d -name "fasterq.tmp.*" | head -n 1)
        if [ -n "$TMPDIR" ]; then
            rm -rf "$TMPDIR"
        fi

        echo "Download complete for {wildcards.srr}"
        """


rule download_all_fastq:
    input:
        expand(f"{OUTDIR}/{{srr}}.fastq", srr=SRR_LIST)

# Pour lancer jusqu'à cette partie : snakemake --use-singularity --singularity-args "-B $(pwd)" --cores 4 download_all_fastq

###############################################
# RULE 2 — Trimming
###############################################

rule build_sif_trimming:
    input:
        definition = config["def_trimming"]
    output:
        image = f"{IMAGEDIR}/trimming.sif"
    params:
        build_cmd = config.get("build_cmd", "apptainer build")
    shell:
        """
        mkdir -p {IMAGEDIR}
        {params.build_cmd} {output.image} {input.definition}
        """

rule trimming:
    input:
        sif="output/images/trimming.sif",
        fastq="output/fastq/{srr}.fastq"
    output:
        trimmed_fastq="output/trimmed/{srr}_trimmed.fq"
    threads: THREADS
    container:
        "output/images/trimming.sif"
    params:
        outdir="output/trimmed",
        quality=20,
        phred=33,
        minlen=25
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


# Pour lancer jusqu'à cette partie : snakemake --use-singularity --singularity-args "-B $(pwd)" --cores 4 trimming_all



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