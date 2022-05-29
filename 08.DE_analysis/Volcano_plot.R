library(argparser, quietly=TRUE)

p <- arg_parser("create volcano plot")

# Add command line arguments
p <- add_argument(p, "--de_result", help = "DE_result obtained form DEseq2", type = "character")
p <- add_argument(p, "--select_genes", help = "a file contain GeneID to label on the plot, one gene one row", type = "character")
p <- add_argument(p, "--FDR_cutoff", help = "significance level", type = "numeric", default = 0.05)
p <- add_argument(p, "--log2FC_cutoff", help = "if absolute value of log2FC is larger than log2FC_cutoff, It's a DEG, default: 1", type = "numeric", default = 1)
p <- add_argument(p, "--height", help="volcano plot height, units: 'in'", type="numeric", default = 5)
p <- add_argument(p, "--width", help="volcano plot width, units: 'in'", type="numeric", default = 6)

# Parse the command line arguments
argv <- parse_args(p)

filename <- argv$de_result
select_geneID <- argv$select_genes
FDR_cutoff <- argv$FDR_cutoff
log2FC_cutoff <- argv$log2FC_cutoff
width <- argv$width
height <- argv$height

#filename <- "./genes.counts.matrix.XiaoYun_4d_vs_Darmor_4d.DESeq2.DE_results"
#select_geneID <- "./select.txt"
#FDR_cutoff <- 0.05
#log2FC_cutoff <- 1
#width <- 6
#height <- 5

library(tidyverse)
library(ggsci)
library(cowplot)

de_result <- read.table(file = filename) %>% rownames_to_column("GeneID") %>%
  select(GeneID, logFC, FDR) %>%
  mutate(direction = if_else(FDR > FDR_cutoff | abs(logFC) < log2FC_cutoff, "NS",
                             if_else(logFC >= log2FC_cutoff, "UP", "DOWN")))
stat <- de_result %>% group_by(direction) %>% count()
up <- stat %>% filter(direction == "UP") %>% pull(n)
down <- stat %>% filter(direction == "DOWN") %>% pull(n)

if (is.na(select_geneID)) {
  Pvolcano <- ggplot(data = de_result, aes(x = logFC, y = -log10(FDR))) +
    geom_point(size = 2, aes(color = direction), show.legend = F) +
    #geom_text(data = res_selected, aes(label = GENE_NAME)) +
    #geom_point(data = res_selected, size = 3, shape = 21, stroke = 1) +
    #geom_text_repel(data = res_selected, aes(label = GENE_NAME)) +
    geom_hline(yintercept = -log10(FDR_cutoff), 
               linetype = 'dotdash', color = 'grey30') +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), 
               linetype = 'dotdash', color = 'grey30') +
    annotate("text", x = max(de_result$logFC), y = max(-log10(de_result$FDR)), 
             label = paste("Up ", up, "\n", "Down ", down, sep = ""), 
             vjust=1, hjust=1, colour="black", size=4) +
    scale_color_manual(values = c('#1613BF', '#727272', '#D5161A')) +
    theme_half_open()
} else {
  library(ggrepel)
  selected_genes <- read_tsv(file = select_geneID, col_names = F) %>% pull(X1)
  res_selected <- filter(de_result, GeneID %in% selected_genes)
  Pvolcano <- ggplot(data = de_result, aes(x = logFC, y = -log10(FDR))) +
    geom_point(size = 2, aes(color = direction), show.legend = F) +
    #geom_text(data = res_selected, aes(label = GeneID)) +
    geom_point(data = res_selected, size = 2, shape = 21, stroke = 1.5) +
    geom_text_repel(data = res_selected, aes(label = GeneID)) +
    geom_hline(yintercept = -log10(FDR_cutoff), 
               linetype = 'dotdash', color = 'grey30') +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), 
               linetype = 'dotdash', color = 'grey30') +
    annotate("text", x = max(de_result$logFC), y = max(-log10(de_result$FDR)), 
             label = paste("Up ", up, "\n", "Down ", down, sep = ""), 
             vjust=1, hjust=1, colour="black", size=4) +
    scale_color_manual(values = c('#1613BF', '#727272', '#D5161A')) +
    theme_half_open()

}

ggsave(Pvolcano, filename = str_replace(filename, "DE_results", "pdf"), width = width, height = height)
ggsave(Pvolcano, filename = str_replace(filename, "DE_results", "png"), width = width, height = height, dpi = 500)

