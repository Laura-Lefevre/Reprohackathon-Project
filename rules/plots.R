#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

resdf_file      <- args[1]
deg_file        <- args[2]
deg_up_file     <- args[3]
deg_down_file   <- args[4]
vsd_file        <- args[5]
gene_pathways_file <- args[6]
annot_long_file <- args[7]
df_enrich_file  <- args[8]
df_top25_file   <- args[9]
outdir          <- args[10]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(forcats)
library(ggrepel)
library(RColorBrewer)

############################################################
# Utility to load CSV with rownames restored
############################################################
load_df <- function(file) {
  df <- read.csv(file, check.names = FALSE)
  
  # remettre les GeneID comme rownames si colonne Geneid existe
  if ("Geneid" %in% colnames(df)) {
    rownames(df) <- df$Geneid
  }
  
  # nettoyer colonnes parasites type X ou Unnamed:0
  bad_cols <- c("X", "Unnamed: 0")
  df <- df[ , !(colnames(df) %in% bad_cols), drop = FALSE]
  
  return(df)
}


############################################################
# LOAD DATA
############################################################
res_df   <- load_df(resdf_file)
deg      <- load_df(deg_file)
deg_up   <- load_df(deg_up_file)
deg_down <- load_df(deg_down_file)

vsd <- readRDS(vsd_file)
mat <- assay(vsd)

gene_pathways <- read.csv(gene_pathways_file, check.names = FALSE)
annot_long    <- read.csv(annot_long_file,    check.names = FALSE)

df_enrich      <- read.csv(df_enrich_file,     check.names = FALSE)
df_top25       <- read.csv(df_top25_file,      check.names = FALSE)

df_top25$pathway <- as.character(df_top25$pathway)

############################################################
# PREP DESEQ2 RESULTS
############################################################
res_df$signif <- !is.na(res_df$padj) & res_df$padj < 0.05
res_df$neglog10padj <- -log10(res_df$padj)

############################################################
# PLOT 1 — MA PLOT
############################################################
p1 <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = signif), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("FALSE"="black","TRUE"="red")) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw()

ggsave(file.path(outdir, "MAplot_global.png"), p1, width = 6, height = 5)

############################################################
# PLOT 2 — VOLCANO
############################################################
p2 <- ggplot(res_df %>% filter(!is.na(padj)),
             aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = signif), size=1, alpha=0.6) +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("FALSE"="black","TRUE"="red")) +
  theme_bw()

ggsave(file.path(outdir, "Volcano.png"), p2, width = 6, height = 5)

############################################################
# PLOT 3 — HEATMAP DEGS
############################################################
genes_deg <- intersect(rownames(mat), rownames(deg))

mat_deg <- mat[genes_deg, , drop = FALSE]
mat_scaled <- t(scale(t(mat_deg)))

pheatmap(
  mat_scaled,
  color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
  show_rownames=FALSE,
  show_colnames=FALSE,
  filename=file.path(outdir,"Heatmap_DEGs.png")
)

############################################################
# PLOT 4 — PATHWAY BARPLOT (≥20 genes)
############################################################
df_counts <- gene_pathways %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep=", ") %>%
  count(pathways, name="gene_count") %>%
  filter(gene_count >= 20)

p3 <- ggplot(df_counts,
             aes(x=reorder(pathways, gene_count), y=gene_count)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  theme_bw()

ggsave(file.path(outdir, "Pathways_barplot.png"), p3, width = 6, height = 7)

############################################################
# PLOT 5 — ENRICHMENT TOP 25 (READ-ONLY)
############################################################
df_top25 <- df_top25 %>% filter(!is.na(pathway))

p_enrich <- ggplot(df_top25,
                   aes(x=odds_ratio,
                       y=fct_reorder(pathway, odds_ratio),
                       color=FDR)) +
  geom_point(aes(size=n_genes)) +
  scale_color_gradient(low="#104E8B", high="#1E90FF") +
  scale_size(range=c(2,8)) +
  theme_bw()

ggsave(file.path(outdir, "enrichment_top25.png"), p_enrich,
       width=7, height=7)

############################################################
# PLOT 6 — MA-PLOT TRANSLATION GENES
############################################################
res_df$locus_tag <- sub(".*?(SAOUHSC_\\d+)", "\\1", res_df$Geneid)
############################################################
# Translation pathways — classification KEGG
############################################################

translation_pathways <- c(
  "sao00970_Aminoacyl-tRNA_biosynthesis",
  "sao03010_Ribosome",
  "sao01240_Biosynthesis_of_cofactors",
  "sao01230_Biosynthesis_of_amino_acids",
  "sao00470_D-Amino_acid_metabolism",
  "sao00250_Alanine,_aspartate_and_glutamate_metabolism",
  "sao00400_Phenylalanine,_tyrosine_and_tryptophan_biosynthesis",
  "sao00350_Tyrosine_metabolism",
  "sao00360_Phenylalanine_metabolism",
  "sao00380_Tryptophan_metabolism",
  "sao00260_Glycine,_serine_and_threonine_metabolism",
  "sao00330_Arginine_and_proline_metabolism",
  "sao00290_Valine,_leucine_and_isoleucine_biosynthesis",
  "sao00280_Valine,_leucine_and_isoleucine_degradation",
  "sao00270_Cysteine_and_methionine_metabolism",
  "sao00340_Histidine_metabolism",
  "sao00310_Lysine_degradation",
  "sao00300_Lysine_biosynthesis",
  "sao00220_Arginine_biosynthesis"
)

############################################################
# Extract translation-related genes from annot_long
############################################################

translation_genes <- annot_long %>%
  filter(pathways %in% translation_pathways) %>%
  distinct(locus_tag) %>%
  pull(locus_tag)

res_translation <- res_df %>%
  filter(locus_tag %in% translation_genes)

############################################################
# Forced labels (tsf, frr, infA/B/C)
############################################################

forced_translation_genes <- tibble::tribble(
  ~locus_tag, ~ids,
  "SAOUHSC_01234","tsf",
  "SAOUHSC_01236","frr",
  "SAOUHSC_02489","infA",
  "SAOUHSC_01246","infB",
  "SAOUHSC_01786","infC"
)

res_translation <- res_translation %>%
  left_join(forced_translation_genes, by="locus_tag")

############################################################
# AA–tRNA synthetase genes (automatic from KEGG annotation)
############################################################

aatRNA_synthetase <- gene_pathways %>%
  filter(grepl("Aminoacyl-tRNA_biosynthesis", pathways)) %>%
  pull(locus_tag) %>%
  unique()

############################################################
# Build plot — translation-only MA-plot with labels
############################################################

p_translation <- ggplot(res_translation,
                        aes(x = log2(baseMean + 1),
                            y = log2FoldChange)) +
  geom_point(aes(color = signif),
             alpha = 0.6, size = 1.3) +
  
  geom_point(data = res_translation %>%
               filter(locus_tag %in% aatRNA_synthetase),
             shape = 21, fill = NA, color = "black",
             stroke = 1.1, size = 3.5) +
  
  geom_point(data = res_translation %>%
               filter(!is.na(ids)),
             color = "black", size = 3) +
  
  geom_label_repel(data = res_translation %>%
                     filter(!is.na(ids)),
                   aes(label = ids),
                   color = "black", size = 4,
                   box.padding = 0.4,
                   max.overlaps = Inf) +
  
  scale_color_manual(values = c("FALSE"="grey60","TRUE"="red")) +
  
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(
    title = "MA-plot – Translation-related genes",
    x = "log2(mean normalized counts)",
    y = "log2 Fold Change"
  ) +
  theme_bw()

ggsave(file.path(outdir, "MAplot_translation.png"),
       p_translation, width = 7, height = 6)
