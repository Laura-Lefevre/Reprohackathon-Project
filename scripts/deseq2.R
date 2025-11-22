#!/usr/bin/env Rscript

########################
# Parse CLI arguments
########################
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- args[i]
    val <- args[i + 1]
    if (key == "--counts")  out$counts  <- val
    if (key == "--samples") out$samples <- val
    if (key == "--outdir")  out$outdir  <- val
  }
  out
}

# This function parses Snakemake script : 
# Rscript scripts/deseq2.R \
# --counts output/counts/all_samples_counts.txt \
# --samples metadata/samples.tsv \
# --outdir results/deseq2
# and creates R list with input/samples/output : 

opt <- parse_args(args)
counts_file  <- opt$counts
# counts_file <- "output/counts/all_samples_counts.txt"
samples_file <- opt$samples
# samples_file <- metadata/samples.tsv
outdir       <- opt$outdir
# outdir <- "results/deseq2"


# Creates DESeq2 result directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

########################
# Load required packages
########################
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(rentrez)
library(EnrichmentBrowser)
library(tidyr)
library(forcats)
library(ggrepel)

########################
# Import featureCounts output
########################
counts <- read.delim(counts_file, comment.char="#", check.names = FALSE)

# extract raw counts
cts <- as.matrix(counts[ , 7:ncol(counts)])
rownames(cts) <- counts$Geneid

########################
# Import sample metadata
########################
# Expected format:
# sample   condition
# SRR001   Normal
# SRR002   Persistent
samples <- read.delim(samples_file, sep="\t", stringsAsFactors = FALSE)

if (!all(samples$sample %in% colnames(cts))) {
  stop("Sample names in samples.tsv do not match count matrix column names.")
}

# reorder columns to match counts file
cts <- cts[ , samples$sample]

# build colData
coldata <- data.frame(condition = factor(samples$condition))
rownames(coldata) <- samples$sample

# set reference level
coldata$condition <- relevel(coldata$condition, ref = "Normal")

########################
# DESeq2 analysis
########################
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  # Model the gene expression **as a function of** the variable 'condition'.
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

########################
# Convert DESeqResults → data.frame
########################

res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)
res_df$signif <- !is.na(res_df$padj) & res_df$padj < 0.05
