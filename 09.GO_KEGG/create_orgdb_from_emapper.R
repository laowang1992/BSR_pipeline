#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(tidyverse)

# Create a parser
p <- arg_parser("make OrgDB from emapper")

# Add command line arguments
p <- add_argument(p, "annotation", help="emapper annotation result", type="character")

# Parse the command line arguments
argv <- parse_args(p)

emapper <- read_delim(argv$annotation, 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X6, 
                GO = X7, 
                KO = X9, 
                Pathway = X10, 
                OG = X21, 
                Gene_Name = X22)

gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))

gene2go <- dplyr::select(emapper, GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  filter(!is.na(GO)) %>%
  mutate(EVIDENCE = 'IEA') 

AnnotationForge::makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               maintainer='zhangsan <zhangsan@genek.tv>',
               author='zhangsan',
               outputDir="./",
               tax_id=0000,
               genus='M',
               species='y',
               goTable="go",
               version="1.0")

pkgbuild::build('.//org.My.eg.db', dest_path = ".")

## 准备GO数据库
dir.create('R_Library', recursive = T)
install.packages('org.My.eg.db_1.0.tar.gz', 
                 repos = NULL, #从本地安装
                 lib = 'R_Library') # 安装文件夹

## 准备 TERM2GENE

emapper <- read_delim(argv$annotation, 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X9, 
                Pathway = X10)

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))


## 准备 TERM2NAME

library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()

## 准备 ko2gene
ko2gene <- dplyr::select(emapper, GID, KO) %>%
  separate_rows(KO, sep = ',', convert = TRUE) %>%
  filter(str_detect(KO, 'ko')) %>%
  mutate(KO = str_remove(KO, 'ko:'))

## 
save(pathway2gene, pathway2name, ko2gene, file = "kegg_info.RData")

