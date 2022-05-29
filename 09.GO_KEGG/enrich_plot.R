#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("enrich results plot")
# Add command line arguments
p <- add_argument(p, "--enrich_result", help="prefix of GO and KEGG result, eg.: genes.counts.matrix.YY_9_vs_ZZ_9.DESeq2.DE_results", type="character")
p <- add_argument(p, "--GO_height", help="GO barplot height, units: 'in'", type="numeric", default = 9)
p <- add_argument(p, "--GO_width", help="GO barplot width, units: 'in'", type="numeric", default = 9)
p <- add_argument(p, "--KEGG_height", help="KEGG botplot height, units: 'in'", type="numeric", default = 9)
p <- add_argument(p, "--KEGG_width", help="KEGG botplot width, units: 'in'", type="numeric", default = 9)
# Parse the command line arguments
argv <- parse_args(p)

## test
#argv <- list()
#argv$enrich_result <- ""

# 
out_prefix <- basename(argv$enrich_result)
GO_file <- paste(argv$enrich_result, "GO_result.csv", sep = ".")
EKGG_file <- paste(argv$enrich_result, "KEGG_result.csv", sep = ".")


library(tidyverse)
library(ggsci)
library(cowplot)

mytheme <- theme(panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.grid.major.y = element_line(colour = "grey",linetype = 2),
                 panel.grid.minor.y = element_line(colour = "grey",linetype = 1),
                 axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 10,
                                            face = "plain", color = "black"),
                 axis.text.y = element_text(face = "plain", color = "black", size = 10),
                 axis.title = element_text(face = "bold", size = 12),
                 legend.text = element_text(size = 9, face = "bold"),
                 legend.title = element_text(size = 12))


## GO绘图

de_ego_df <- read_csv(file = GO_file)
de_ego_df

de_ego_df <- de_ego_df %>% 
  separate(GeneRatio, c("term_deg", "ann_deg"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio = term_deg / ann_deg)

p1 <- ggplot(de_ego_df, aes(x = reorder(Description, p.adjust), y = GeneRatio)) +
  geom_bar(aes(fill = -log10(p.adjust)), color="black", stat= 'identity',width = 0.7) +
  geom_text(aes(x = Description, y = GeneRatio), label = de_ego_df$term_deg, vjust=-1,size=3.5) +
  scale_fill_gradient(low = "blue", high = "red") +
  ylim(0, max(de_ego_df$GeneRatio) * 1.2) +
  labs(x = NULL, y = "Gene Ratio", fill = "-log10 P.ajust") +
  facet_grid(. ~ ONTOLOGY, scales = "free_x",space = "free_x") +
  theme_half_open() +
  mytheme
p1
ggsave(p1, filename = paste(out_prefix, "GO_barplot1.pdf", sep = "."), height = argv$GO_height, width = argv$GO_width)
ggsave(p1, filename = paste(out_prefix, "GO_barplot1.png", sep = "."), height = argv$GO_height, width = argv$GO_width, dpi = 500)

p2 <- ggplot(de_ego_df, aes(x = reorder(Description, p.adjust), y = GeneRatio)) +
  geom_bar(aes(fill = ONTOLOGY), color="black", stat= 'identity',width = 0.7) +
  geom_text(aes(x = Description, y = GeneRatio), label = de_ego_df$term_deg, vjust=-1,size=3.5) +
  scale_fill_aaas() +
  ylim(0, max(de_ego_df$GeneRatio) * 1.2) +
  labs(x = NULL, y = "Gene Ratio") +
  facet_grid(. ~ ONTOLOGY, scales = "free_x",space = "free_x") +
  theme_half_open() +
  mytheme

#argv$GO_height <- 9
#argv$GO_width <- 9
ggsave(p2, filename = paste(out_prefix, "GO_barplot2.pdf", sep = "."), height = argv$GO_height, width = argv$GO_width)
ggsave(p2, filename = paste(out_prefix, "GO_barplot2.png", sep = "."), height = argv$GO_height, width = argv$GO_width, dpi = 500)

## KEGG绘图
de_ekp_df <- read_csv(file = KEGG_file)
de_ekp_df

de_ekp_df <- de_ekp_df %>% 
  separate(GeneRatio, c("term_deg", "ann_deg"), sep = "/", convert = TRUE) %>% 
  mutate(GeneRatio = term_deg / ann_deg)

p3 <- ggplot(de_ekp_df, aes(x = reorder(Description, GeneRatio), y = GeneRatio)) +
#  geom_segment(aes(xend=Description, yend=0), color = "orange") +
  geom_point(aes(size = term_deg, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = NULL, y = "Gene Ratio", color = "-log10 P.ajust", size = "Count") +
  coord_flip() +
  theme_bw() +
  mytheme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
#argv$KEGG_height <- 9
#argv$KEGG_width <- 9
ggsave(p3, filename = paste(out_prefix, "KEGG_dotplot.pdf", sep = "."), height = argv$KEGG_height, width = argv$KEGG_width)
ggsave(p3, filename = paste(out_prefix, "KEGG_dotplot.png", sep = "."), height = argv$KEGG_height, width = argv$KEGG_width, dpi = 500)
