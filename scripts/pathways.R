#!/usr/bin/env Rscript

######################################################
### Reprohackathon Project - 2025                  ###
### Author : Saïda MOUSSAEVA, Camille LE CORRE     ###
### Description : R file to get pathways of        ###
### each gene                                      ###
######################################################

library(dplyr)
library(tidyr)
library(rentrez)
library(EnrichmentBrowser)

extract_tag <- function(x) sub(".*?(SAOUHSC_\\d+)", "\\1", x)
get_gene_ids <- function(tags){
  sapply(tags, function(tag){
    out <- entrez_search(db="gene", term=paste0(tag,"[locus_tag]"))
    if(length(out$ids)==0) return(NA)
    out$ids[1]
  })
}

# full res_df available
res_df <- read.csv("~/Bureau/reprohack/res_df.csv", stringsAsFactors = FALSE)
clean_tags <- extract_tag(res_df$Geneid)
ncbi_ids   <- get_gene_ids(clean_tags)

gene_ids_df <- data.frame(locus_tag=clean_tags, gene_id=ncbi_ids)
write.csv(gene_ids_df, "gene_ids.csv", row.names=FALSE)

pathways <- getGenesets(org="sao", db="kegg")
saveRDS(pathways, "kegg_pathways.rds")

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
  pathways  = pw_list
)
write.csv(gene_pathways, "gene_pathways.csv", row.names=FALSE)