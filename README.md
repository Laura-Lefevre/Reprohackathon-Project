# Reprohackathon-Project

## Overview

This repository provides a fully reproducible computational workflow designed to replicate the RNA-seq differential-expression results presented in:

Peyrusson, F., Varet, H., Nguyen, T.K. et al.   
**Intracellular Staphylococcus aureus persisters upon antibiotic exposure.**  
Nat Commun 11, 2200 (2020).   
https://doi.org/10.1038/s41467-020-15966-7  


The workflow reconstructs the preprocessing steps applied to the raw sequencing data, using a controlled and modular Snakemake pipeline in which every rule executes inside a dedicated Singularity container. This architecture ensures strict reproducibility, environment isolation, and long-term traceability.  
The referenced study investigates dormancy induction in *Staphylococcus aureus* under metabolic stress, using RNA-seq to characterize transcriptional shifts associated with persistence. Replicating such results requires exact control over software versions, trimming parameters, and preprocessing logic.

## Design choices for reproducibility  

This repository is structured to make the analysis as deterministic and auditable as possible.

### Containers  
Each major step (download, trimming, mapping, counting, DESeq2, enrichment) runs in its own Singularity image.

Tool versions are pinned inside each image.

Parameters are passed from Snakemake/config rather than hard-coded in the container.

Re-running the pipeline in the future reuses exactly the same software stack as in this ReproHackathon, maximising proximity to the original article’s results.

### Rules

Snakemake rules are split across multiple .rules files under `rules/` (e.g. download, trimming, mapping, counting, DESeq2, enrichment).

Improves readability and makes it clear which piece of code is responsible for each stage.

Allows local modification or extension of a single module (e.g. adding a new mapper or an alternative trimming strategy) without touching the rest of the workflow.

Facilitates review: each file can be checked against the methods section of the paper.

### Centralised configuration

All user-tunable parameters are stored in `config.yaml` (quality threshold, PHRED encoding, minimal read length, number of threads, paths to reference genome/annotation, container image list, etc.).

The logic in rules/ does not depend on hard-coded paths or thresholds.

Reproducing or adapting the analysis to a new dataset mainly requires editing config.yaml, not the workflow code.


## Execution

To run the workflow, execute the provided script:

`run.sh` 

The script positions itself in the project directory, creates the base logs/ folder, and launches Snakemake with Singularity enabled. Snakemake then resolves dependencies, builds the directed acyclic graph (DAG), and executes each rule in the correct order inside the corresponding Singularity container.


## Repository structure 

├── containers/            # Singularity images used by rules  
├── enrichment_data/       # data used for enrichment analyses (pathways, gene IDs)  
├── metadata/              # sample metadata (e.g. sample list, annotation)  
├── rules/                 # Snakemake rule definitions (e.g. trimming, mapping)    
├── scripts/               # scripts   (mainly .R)  
├── config.yaml            # configuration file (genomes, quality threshold, min length, etc.)  
├── Snakefile              # main workflow definition    
├── .gitignore    
└── README.md              # this file  

## Requirements

`Snakemake v.7.32.4`  

`Singularity/Apptainer v.1.4.3`

## Outputs

This pipeline generates the folowwing output structure : 

├── outputs  
  ├── fastq               # Raw FASTQ files downloaded for all samples  
  ├── mapping             # Aligned reads: .bam files and index files  
  ├── trimmed             # FASTQ files after Cutadapt trimming (quality/length filtered)  
  ├── counts              # Gene-level count tables (featureCounts outputs)  
  ├── deseq2              # DESeq2 results: DEG tables, .csv files  
  ├── enrichment          # Pathway/GO enrichment results (Fisher tests, filtered sets)  
  ├── plots               # All generated figures (MA plot, volcano, PCA, heatmaps, enrichment plots)  
  ├── vsd.rds             # Variance-stabilized data object (VST/VSD)  
├── logs                  # Logs for trimming, mapping and counting rules  

