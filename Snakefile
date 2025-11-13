###############################################
# Snakefile — Reprohackathon Project (version monolithique)
###############################################

# Chargement du fichier de configuration
configfile: "config.yaml"

#include : rules/sif_images.rules 

THREADS = config["THREADS"]

#SRR_LIST=["SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726"]
SRR_LIST=["SRR10379721"]

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

rule build_sif_mapping:
    input:
        definition = config["def_mapping"]
    output:
        image = f"{IMAGEDIR}/mapping.sif"
    params:
        build_cmd = config.get("build_cmd", "apptainer build")
    shell:
        """
        mkdir -p {IMAGEDIR}
        {params.build_cmd} {output.image} {input.definition}
        """

# Règle pour télécharger le génome de référence
rule download_ref_genome:
    input:
        sif = "output/images/mapping.sif"
    output:
        "output/genome/reference.fasta"
    container:
        "output/images/mapping.sif"
    shell:
        "wget -q -O {output} 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta'"

# Règle pour télécharger les annotations
rule download_ref_annot:
    input:
        sif = "output/images/mapping.sif"
    output:
        "output/genome/reference.gff"
    container:
        "output/images/mapping.sif"
    shell:
        "wget -O {output} 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1'"

# Règle pour créer l'index Bowtie
rule bowtie_build_index:
    input:
        sif = "output/images/mapping.sif",
        fasta = "output/genome/reference.fasta"
    output:
        touch("output/genome/index.done")
    params:
        prefix="output/genome/reference_index"
    container:
        "output/images/mapping.sif"
    shell:
        "bowtie-build {input.fasta} {params.prefix} && touch {output}"

# Règle pour le mapping
rule bowtie_map:
    input:
        sif = "output/images/mapping.sif",
        fastq="output/trimmed/{srr}_trimmed.fq",
        index_done="output/genome/index.done"
    output:
        bam="output/mapping/{srr}.bam",
        bai="output/mapping/{srr}.bam.bai"
    log:
        "logs/mapping/{srr}.log"
    params:
        index_prefix="output/genome/reference_index",
        outdir="output/mapping"
    threads: THREADS
    container:
        "output/images/mapping.sif"
    shell:
        """
        mkdir -p {params.outdir}
        bowtie -S -p {threads} {params.index_prefix} {input.fastq} 2> {log} |
        samtools sort -@ {threads} > {output.bam}
        samtools index {output.bam}
        """

rule mapping_all:
    input:
        expand("output/mapping/{srr}.bam", srr=SRR_LIST),
        expand("output/mapping/{srr}.bam.bai", srr=SRR_LIST)

#Pour lancer jusqu'à cette partie : snakemake --use-singularity --cores 4 mapping_all

###############################################
# RULE 4 — Counting
###############################################
# Build SIF for featureCounts
rule build_sif_counting:
    input:
        definition = config["def_counting"]  # ex: containers/featurecounts.def
    output:
        image = "containers/featureCounts.sif"
    params:
        build_cmd = config.get("build_cmd", "apptainer build")
    shell:
        """
        mkdir -p containers
        {params.build_cmd} {output.image} {input.definition}
        """

# Counting with featureCounts
rule counting:
    input:
        sif  = "output/images/featurecounts.sif",
        bams = expand("output/mapping/{srr}.bam", srr=SRR_LIST),
        gff  = "output/genome/reference.gff"
    output:
        "output/counts/all_samples_counts.txt"
    container:
        "output/images/featurecounts.sif"
    params:
        outdir = "output/counts"
    threads: THREADS
    shell:
        """
        echo "Running featureCounts for all samples..."
        mkdir -p {params.outdir}
        featureCounts -t gene -g ID -s 1 \
                      -a {input.gff} \
                      -o {output} \
                      {input.bams}
        echo "Combined featureCounts results saved to {output}"
        """

#rule counting_all:
#    input:
 #       "output/counts/all_samples_counts.txt"

# Pour lancer jusqu'à cette partie : 
# snakemake --use-singularity --cores 4 mapping_all
# snakemake --cores 4 build_sif_counting   # si l'image n'est pas encore construite
# snakemake --use-singularity --cores 4 counting




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

# include : "rules/DESeq2.rules"


