#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("do GO/KEGG enrichment, need prepare kegg_info.RData(containing ko2gene pathway2gene and pathway2name) and org.*.eg.db_1.0.tar.gz in $work_dir/db")

# Add command line arguments
p <- add_argument(p, "--de_result", help="input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--de_log2FoldChange", help="log2FoldChange cutoff", type="numeric", default = 1)
p <- add_argument(p, "--de_fdr", help="FDR cutoff", type="numeric", default = 0.05)
p <- add_argument(p, "--enrich_pvalue", help="pvalue cutoff for enrichment", type="numeric", default = 0.05)
p <- add_argument(p, "--enrich_qvalue", help="qvalue cutoff for enrichment", type="numeric", default = 0.05)
p <- add_argument(p, "--orgdb", help="orgdb package name, must have been installed in ./R_library", type="character")
p <- add_argument(p, "--species", help="character, either the kegg code, scientific name or the common name of the target species.", type="character", default = "ko")


# Parse the command line arguments
argv <- parse_args(p)

########################################################################
## test
test <- FALSE
if (test) {
  argv <- list()
  argv$de_result <- "genes.counts.matrix.PG_vs_PP.edgeR.DE_results"
  argv$de_log2FoldChange <- 1
  argv$de_fdr <- 0.05
  argv$enrich_pvalue <- 0.05
  argv$enrich_qvalue <- 0.05
  argv$orgdb <- "org.My.eg.db"
  argv$species <- "ko"
}

########################################################################

out_prefix <- basename(argv$de_result)
orgdb <- argv$orgdb
# load library ------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)

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

de_result <- read.table(file = argv$de_result)
gene <- filter(de_result,
               abs(logFC) > argv$de_log2FoldChange & FDR < argv$de_fdr) %>%
  rownames(id)

geneList <- de_result$logFC
names(geneList) <- rownames(de_result)
geneList <- sort(geneList, decreasing = TRUE)

de_ego <- enrichGO(gene = gene,
                   OrgDb = orgdb,
                   keyType = 'GID',
                   ont = 'ALL',
                   qvalueCutoff = argv$enrich_qvalue,
                   pvalueCutoff = argv$enrich_pvalue)
de_ego_df <- as.data.frame(de_ego)
head(de_ego_df)
write_csv(de_ego_df, paste(out_prefix, "GO_result", "csv", sep = "."))
write_tsv(de_ego_df, paste(out_prefix, "GO_result", "txt", sep = "."))

# GO绘图
pdf(file = paste(out_prefix, "GO_barplot.pdf", sep = "."))
barplot(de_ego, showCategory=15, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")
dev.off()

png(filename = paste(out_prefix, "GO_barplot.png", sep = "."))
barplot(de_ego, showCategory=15, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")
dev.off()

## KEGG
load("../db/kegg_info.RData")
de_ekp <- enricher(gene,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = argv$enrich_pvalue,
                   qvalueCutoff = argv$enrich_qvalue)
de_ekp_df <- as.data.frame(de_ekp)
head(de_ekp_df)
write_csv(de_ekp_df, paste(out_prefix, "KEGG_result", "csv", sep = "."))
write_tsv(de_ekp_df, paste(out_prefix, "KEGG_result", "txt", sep = "."))

# KEGG绘图
pdf(file = paste(out_prefix, "KEGG_dotplot.pdf", sep = "."))
dotplot(de_ekp, showCategory=20, x = "GeneRatio")
dev.off()
png(filename = paste(out_prefix, "KEGG_dotplot.png", sep = "."))
dotplot(de_ekp, showCategory=20, x = "GeneRatio")
dev.off()

pdf(file = paste(out_prefix, "KEGG_emapplot.pdf", sep = "."))
emapplot(de_ekp)
dev.off()
png(filename = paste(out_prefix, "KEGG_emapplot.png", sep = "."))
emapplot(de_ekp)
dev.off()

pdf(file = paste(out_prefix, "KEGG_cnetplot.pdf", sep = "."))
cnetplot(de_ekp, 
         foldChange = geneList, 
         showCategory = 5,
         node_label = "category", # category | gene | all | none
         #circular = TRUE, 
         colorEdge = TRUE)
dev.off()
png(file = paste(out_prefix, "KEGG_cnetplot.png", sep = "."))
cnetplot(de_ekp, 
         foldChange = geneList, 
         showCategory = 5,
         node_label = "category", # category | gene | all | none
         #circular = TRUE, 
         colorEdge = TRUE)
dev.off()

# pathway view ------------------------------------------------------------
# 此处并不完美，并不是对应物种专用的pathway，有改善空间
#id.map <- select(org.My.eg.db, keys = names(geneList), columns = "Ko")
#gene.ko <- mol.sum(mol.data = geneList, id.map = ko2gene)
gene.ko <- geneList
id <- tibble(GID = names(gene.ko)) %>% left_join(ko2gene, by = "GID")
names(gene.ko) <- id$KO
gene.ko <- as.matrix(gene.ko[!is.na(names(gene.ko))])

sig.pathway <- as.character(filter(de_ekp_df, p.adjust < argv$enrich_qvalue)$ID)

work_dir <- getwd()
pathview_dir <- paste(out_prefix, 'pathwiew', sep = "_")
try(dir.create(pathview_dir, recursive=T))
setwd(pathview_dir)

result <- tryCatch(pathview(gene.data  = gene.ko,
         pathway.id = str_replace(sig.pathway, "bna", "ko"),
         species    = "ko"), error = function(e){"Can't download pathway map, maybe network is not available!\n"})

setwd(work_dir)
