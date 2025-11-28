#!/usr/bin/env Rscript

######################################################
### Reprohackathon Project - 2025                  ###
### Author : Saïda MOUSSAEVA                       ###
### Description : R file for Differential          ###
### expression analysis with DESeq2s               ###
######################################################

###############
# Parse args
###############
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- args[i]
    val <- args[i+1]
    if (key == "--counts")  out$counts  <- val
    if (key == "--samples") out$samples <- val
    if (key == "--outdir")  out$outdir  <- val
    if (key == "--vsd")     out$vsd     <- val
  }
  out
}

opt <- parse_args(args)
counts_file  <- opt$counts
samples_file <- opt$samples
outdir       <- opt$outdir
vsd_file     <- opt$vsd

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###############
# Libraries
###############
library(DESeq2)

###############
# Load counts
###############
counts <- read.delim(counts_file, comment.char="#", check.names = FALSE)

cts <- as.matrix(counts[,7:ncol(counts)])
rownames(cts) <- counts$Geneid

samples <- read.delim(samples_file, sep="\t", stringsAsFactors = FALSE)
cts <- cts[, samples$sample]

coldata <- data.frame(condition = factor(samples$condition))
rownames(coldata) <- samples$sample
coldata$condition <- relevel(coldata$condition, ref = "Normal")

###############
# DESeq2
###############
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

###############
# Dataframes
###############
res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)
res_df$signif <- !is.na(res_df$padj) & res_df$padj < 0.05

deg_df      <- res_df[res_df$signif, ]
deg_up_df   <- deg_df[deg_df$log2FoldChange > 0, ]
deg_down_df <- deg_df[deg_df$log2FoldChange < 0, ]

###############
# VST
###############
vsd <- vst(dds, blind = FALSE)

###############
# Write outputs — ALWAYS row.names = FALSE
###############
write.csv(res_df,      file.path(outdir, "res_df.csv"),      row.names = FALSE)
write.csv(deg_df,      file.path(outdir, "deg.csv"),         row.names = FALSE)
write.csv(deg_up_df,   file.path(outdir, "deg_up.csv"),      row.names = FALSE)
write.csv(deg_down_df, file.path(outdir, "deg_down.csv"),    row.names = FALSE)

saveRDS(vsd, vsd_file)



