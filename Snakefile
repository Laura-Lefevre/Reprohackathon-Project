###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

# Chargement du fichier de configuration
configfile: "config.yaml"

THREADS = config["THREADS"]

SRR_LIST=["SRR10510434", "SRR10510435", "SRR10510436", "SRR10510437", "SRR10510438", "SRR10510439"]

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

# Build SIF for DESeq2
rule build_sif_deseq2:
    input:
        definition = config["def_deseq2"]  # ex: containers/DESeq2.def
    output:
        image = "containers/DESeq2.sif"
    params:
        build_cmd = config.get("build_cmd", "apptainer build")
    shell:
        """
        mkdir -p containers
        {params.build_cmd} {output.image} {input.definition}
        """

# Differential expression with DESeq2
# Attend un tableau de comptes issu de featureCounts (gènes en lignes, samples en colonnes),
# et une table d'échantillons (colonnes: sample, condition).
rule deseq2:
    input:
        sif     = "containers/DESeq2.sif",
        counts  = "output/counts/featureCounts.tsv",   # output rule 4
        samples = "metadata/samples.tsv",              # mapping sample->condition
        script  = "scripts/deseq2.R"                   # ton script R
    output:
        results = "results/deseq2/deseq2_results.csv",
        rds     = "results/deseq2/dds.rds"
    container:
        "containers/DESeq2.sif"
    params:
        outdir  = "results/deseq2"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {input.script} \
            --counts  {input.counts} \
            --samples {input.samples} \
            --outdir  {params.outdir}
        """




