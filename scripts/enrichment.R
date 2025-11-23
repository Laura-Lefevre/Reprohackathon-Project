#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
resdf_file  <- args[1]
outdir      <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(rentrez)
library(EnrichmentBrowser)

###########################################################
# LOAD RES_DF
###########################################################
res_df <- read.csv(resdf_file, stringsAsFactors = FALSE)

all_genes <- res_df$Geneid

###########################################################
# GENE → NCBI ID
###########################################################
extract_tag <- function(x) sub(".*?(SAOUHSC_\\d+)", "\\1", x)

get_gene_ids <- function(tags){
  sapply(tags, function(tag){
    out <- entrez_search(db="gene", term=paste0(tag,"[locus_tag]"))
    if(length(out$ids)==0) return(NA)
    out$ids[1]
  })
}

clean_tags <- extract_tag(all_genes)
ncbi_ids   <- get_gene_ids(clean_tags)

gene_ids_df <- data.frame(
  locus_tag = clean_tags,
  gene_id   = ncbi_ids,
  stringsAsFactors = FALSE
)

write.csv(gene_ids_df,
          file.path(outdir,"gene_ids.csv"),
          row.names = FALSE)

###########################################################
# GET PATHWAYS
###########################################################
pathways <- getGenesets(org="sao", db="kegg")

get_pw <- function(id){
  if(is.na(id)) return(NA)
  pw <- names(Filter(function(gs) id %in% gs, pathways))
  if(length(pw)==0) return(NA)
  paste(pw, collapse=", ")
}

pw_list <- sapply(ncbi_ids, get_pw)

gene_pathways <- data.frame(
  locus_tag = clean_tags,
  gene_id   = ncbi_ids,
  pathways  = pw_list,
  stringsAsFactors = FALSE
)

write.csv(gene_pathways,
          file.path(outdir,"gene_pathways.csv"),
          row.names = FALSE)

###########################################################
# LONG FORMAT
###########################################################
annot_long <- gene_pathways %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep=", ") %>%
  distinct()

write.csv(annot_long,
          file.path(outdir,"annot_long.csv"),
          row.names = FALSE)

###########################################################
# ENRICHMENT (FISHER)
###########################################################
deg <- res_df %>% filter(!is.na(padj) & padj < 0.05)
deg_genes <- extract_tag(deg$Geneid)
all_vec   <- extract_tag(res_df$Geneid)

enrich_one <- function(path){
  genes_in <- annot_long$locus_tag[annot_long$pathways == path]
  
  a = sum(deg_genes %in% genes_in)
  b = sum(!(deg_genes %in% genes_in))
  c = sum(all_vec %in% genes_in) - a
  d = length(all_vec) - (a+b+c)
  
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

df_enrich <- bind_rows(lapply(all_paths, enrich_one))
df_enrich$FDR <- p.adjust(df_enrich$pvalue, "BH")

write.csv(df_enrich,
          file.path(outdir,"df_enrich.csv"),
          row.names = FALSE)

df_top25 <- df_enrich %>%
  arrange(desc(odds_ratio)) %>%
  head(25)

write.csv(df_top25,
          file.path(outdir,"df_enrich_top25.csv"),
          row.names = FALSE)

