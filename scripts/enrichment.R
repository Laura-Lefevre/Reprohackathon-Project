
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
resdf_file  <- args[1]
outdir      <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)

###########################################################
# LOAD INPUT
###########################################################
res_df <- read.csv(resdf_file, stringsAsFactors = FALSE)

# Extraction locus tags
extract_tag <- function(x) sub(".*?(SAOUHSC_\\d+)", "\\1", x)
clean_tags  <- extract_tag(res_df$Geneid)

###########################################################
# LOAD OFFLINE DATA (pré-générée une fois)
###########################################################
gene_ids_df      <- read.csv("gene_ids.csv", stringsAsFactors = FALSE)
gene_pathways_df <- read.csv("gene_pathways.csv", stringsAsFactors = FALSE)
pathways_list    <- readRDS("kegg_pathways.rds")   # structure KEGG complète

###########################################################
# WRITE OFFLINE GENE IDS (direct copy)
###########################################################
# On restreint aux gènes présents dans res_df
sub_gene_ids <- gene_ids_df %>%
  filter(locus_tag %in% clean_tags)

write.csv(
  sub_gene_ids,
  file.path(outdir, "gene_ids.csv"),
  row.names = FALSE
)

###########################################################
# WRITE OFFLINE GENE → PATHWAY (direct copy)
###########################################################
sub_gene_pathways <- gene_pathways_df %>%
  filter(locus_tag %in% clean_tags)

write.csv(
  sub_gene_pathways,
  file.path(outdir, "gene_pathways.csv"),
  row.names = FALSE
)

###########################################################
# LONG FORMAT
###########################################################
annot_long <- sub_gene_pathways %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep = ", ") %>%
  distinct()

write.csv(
  annot_long,
  file.path(outdir, "annot_long.csv"),
  row.names = FALSE
)

###########################################################
# ENRICHMENT (FISHER)
###########################################################
deg <- res_df %>% filter(!is.na(padj) & padj < 0.05)
deg_genes <- extract_tag(deg$Geneid)
all_vec   <- clean_tags

enrich_one <- function(path){
  genes_in <- annot_long$locus_tag[annot_long$pathways == path]
  
  a = sum(deg_genes %in% genes_in)
  b = sum(!(deg_genes %in% genes_in))
  c = sum(all_vec %in% genes_in) - a
  d = length(all_vec) - (a + b + c)
  
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
  
  data.frame(
    pathway = path,
    odds_ratio = ft$estimate,
    pvalue     = ft$p.value,
    n_genes    = a,
    stringsAsFactors = FALSE
  )
}

all_paths <- unique(annot_long$pathways)

if (length(all_paths) == 0) {
  df_enrich <- data.frame(
    pathway=character(0),
    odds_ratio=numeric(0),
    pvalue=numeric(0),
    n_genes=integer(0),
    FDR=numeric(0)
  )
} else {
  df_enrich <- bind_rows(lapply(all_paths, enrich_one))
  df_enrich$FDR <- p.adjust(df_enrich$pvalue, "BH")
}

write.csv(
  df_enrich,
  file.path(outdir, "df_enrich.csv"),
  row.names = FALSE
)

df_top25 <- df_enrich %>% arrange(desc(odds_ratio)) %>% head(25)

write.csv(
  df_top25,
  file.path(outdir, "df_enrich_top25.csv"),
  row.names = FALSE
)
