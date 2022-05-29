#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("GSEA, need prepare kegg_info.RData(containing ko2gene pathway2gene and pathway2name) and org.*.eg.db_1.0.tar.gz in $work_dir/db,")

# Add command line arguments
p <- add_argument(p, "--de_result", help = "input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--enrich_pvalue", help = "pvalue cutoff for enrichment", type="numeric", default = 0.05)
p <- add_argument(p, "--orgdb", help = "orgdb package name, must have been installed in ./R_library", type="character")
p <- add_argument(p, "--species", help = "character, either the kegg code, scientific name or the common name of the target species.", type="character", default = "ko")
p <- add_argument(p, "--drawPdf", help = "draw pdf image or not", type = "logical", default = "FALSE")

# Parse the command line arguments
argv <- parse_args(p)

########################################################################
## test
test <- F
if (test) {
  argv <- list()
  argv$de_result <- "./DE/genes.counts.matrix.F64_11_vs_F64_22.DESeq2.DE_results"
  argv$enrich_pvalue <- 0.05
  argv$orgdb <- "org.My.eg.db"
  argv$species <- "ko"
  argv$drawPdf <- FALSE
}

########################################################################

filename <- basename(argv$de_result)
orgdb <- argv$orgdb
enrich_pvalue <- argv$enrich_pvalue
draw_pdf <- argv$drawPdf
species <- argv$species

## load R packages
library(clusterProfiler)
library(enrichplot)
library(tidyverse)

if (!dir.exists("../db/R_Library/")) {
  dir.create('../db/R_Library', recursive = T)
}
if (!requireNamespace(orgdb, lib.loc = "../db/R_Library/", quietly = TRUE)) {
  orgpkg <- list.files(path = "../db/", pattern = "org.*eg.db_.*tar\\.gz")[[1]]
  install.packages(paste("../db/", orgpkg, sep = ""), 
                   repos = NULL, #从本地安装
                   lib = '../db/R_Library') # 安装文件夹
}
library(orgdb, lib.loc = "../db/R_Library", character.only = TRUE)

# read DE result
de_result <- read.table(argv$de_result)
genes <- de_result$logFC
names(genes) <- rownames(de_result)
genes <- sort(genes, decreasing = T)

## GO
# GSEA
gseGO_res <- gseGO(geneList = genes,
                   OrgDb = get(orgdb),
                   ont = "ALL",
                   keyType = "GID",
                   minGSSize = 5,
                   maxGSSize = 1000,
                   pvalueCutoff = enrich_pvalue,
                   pAdjustMethod = "BH")

# write GO result
write_tsv(gseGO_res@result, paste(filename, "gseGO.txt", sep = "."))
write_csv(gseGO_res@result, paste(filename, "gseGO.csv", sep = "."))

# create GSEA GO dictionary and plot
if (!dir.exists(paste(filename, "GO", sep = "/"))) {
  dir.create(paste(filename, "GO", sep = "/"), recursive = T)
}
for (i in 1:length(gseGO_res@result$ID)) {
  cat("Plot ", gseGO_res@result$ID[i], " ", gseGO_res@result$Description[i], "\n", sep = "")
    prefix <- str_replace(gseGO_res@result$ID[i], ":", "")
  p <- gseaplot2(gseGO_res, geneSetID = gseGO_res@result$ID[i], title = gseGO_res@result$Description[i])
  if(draw_pdf){
    pdf(file = paste(filename, "/GO/", prefix, ".pdf", sep = ""), width = 6, height = 4.5)
    print(p)
    dev.off()
  }
  png(file = paste(filename, "/GO/", prefix, ".png", sep = ""), width = 6, height = 4.5, units = "in", res = 500)
  print(p)
  dev.off()
}

#gseaplot2(gseGO_res, geneSetID = 4, title = gseGO_res@result$Description[4])

## KEGG
load("../db/kegg_info.RData")
# GSEA
gseKEGG_res <- GSEA(geneList = genes,
                    TERM2GENE = pathway2gene,
                    TERM2NAME = pathway2name,
                    minGSSize = 5,
                    maxGSSize = 1000,
                    pvalueCutoff = enrich_pvalue,
                    pAdjustMethod = "BH")
# write KEGG result
write_tsv(gseKEGG_res@result, paste(filename, "gseKEGG.txt", sep = "."))
write_csv(gseKEGG_res@result, paste(filename, "gseKEGG.csv", sep = "."))

# create GSEA KEGG dictionary and plot
if (!dir.exists(paste(filename, "KEGG", sep = "/"))) {
  dir.create(paste(filename, "KEGG", sep = "/"), recursive = T)
}
for (i in 1:length(gseKEGG_res@result$ID)) {
  cat("Plot ", gseKEGG_res@result$ID[i], " ", gseKEGG_res@result$Description[i], "\n", sep = "")
  p <- gseaplot2(gseKEGG_res, geneSetID = gseKEGG_res@result$ID[i], title = gseKEGG_res@result$Description[i])
  if(draw_pdf){
    pdf(file = paste(filename, "/KEGG/", gseKEGG_res@result$ID[i], ".pdf", sep = ""), width = 6, height = 4.5)
    print(p)
    dev.off()
  }
  png(file = paste(filename, "/KEGG/", gseKEGG_res@result$ID[i], ".png", sep = ""), width = 6, height = 4.5, units = "in", res = 500)
  print(p)
  dev.off()
}

