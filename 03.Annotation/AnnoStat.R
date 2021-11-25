#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for annovar results statistic")

## Add command line arguments
#
p <- add_argument(p, "--snpvar", help = "<filename>.SNPs.anno.variant_function", type = "character")
p <- add_argument(p, "--indelvar", help = "<filename>.INDELs.anno.variant_function", type = "character")

p <- add_argument(p, "--snpex", help = "<filename>.SNPs.anno.exonic_variant_function", type = "character")
p <- add_argument(p, "--indelex", help = "<filename>.INDELs.anno.exonic_variant_function", type = "character")

# Parse the command line arguments
argv <- parse_args(p)

snpvarFilename <- argv$snpvar
indelvarFilename <- argv$indelvar
snpexFilename <- argv$snpex
indelexFilename <- argv$indelex

test <- FALSE
if(test){
  snpvarFilename <- "./Tae_BSR.filter.SNPs.anno.variant_function"
  indelvarFilename <- "./Tae_BSR.filter.INDELs.anno.variant_function"
  snpexFilename <- "./Tae_BSR.filter.SNPs.anno.exonic_variant_function"
  indelexFilename <- "./Tae_BSR.filter.INDELs.anno.exonic_variant_function"
}

library(tidyverse)
library(cowplot)
library(RColorBrewer)

## 
# snp stat
snpvar <- read_tsv(file = snpvarFilename, col_names = F) %>%
  mutate(varType = "SNP")

snpvar_stat <- snpvar %>% 
  group_by(X1, varType) %>%
  summarise(number = n()) %>%
  rename(Distribution = X1)

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
Psnpvar <- ggplot(snpvar_stat, aes(x = Distribution, y = number)) +
  geom_bar(aes(fill = Distribution), stat="identity", alpha = 0.6) +
  scale_fill_manual(values = getPalette(length(unique(snpvar$X1)))) +
  labs(x = NULL, y = "Number") +
  theme_classic(base_line_size = 1) +
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "NA"
  )
ggsave(Psnpvar, filename = str_replace(basename(snpvarFilename), "anno\\.variant_function", "distribution.pdf"), width = 6, height = 5)
ggsave(Psnpvar, filename = str_replace(basename(snpvarFilename), "anno\\.variant_function", "distribution.png"), width = 6, height = 5, dpi = 500)

# indel stat
indelvar <- read_tsv(file = indelvarFilename, col_names = F) %>%
  mutate(varType = "INDEL")

indelvar_stat <- indelvar %>% 
  group_by(X1, varType) %>%
  summarise(number = n()) %>% 
  rename(Distribution = X1)

getPalette = colorRampPalette(brewer.pal(8, "Set1"))
Pindelvar <- ggplot(indelvar_stat, aes(x = Distribution, y = number)) +
  geom_bar(aes(fill = Distribution), stat="identity", alpha = 1) +
  scale_fill_manual(values = getPalette(length(unique(indelvar$X1)))) +
  labs(x = NULL, y = "Number") +
  theme_classic(base_line_size = 1) +
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "NA"
  )
ggsave(Pindelvar, filename = str_replace(basename(indelvarFilename), "anno\\.variant_function", "distribution.pdf"), width = 6, height = 5)
ggsave(Pindelvar, filename = str_replace(basename(indelvarFilename), "anno\\.variant_function", "distribution.png"), width = 6, height = 5, dpi = 500)

# 
var_stat <- rbind(snpvar_stat, indelvar_stat)
Pvar <- ggplot(var_stat, aes(x = Distribution, y = number)) +
  geom_bar(aes(fill = Distribution, alpha = varType), stat="identity", position = "dodge") +
  scale_fill_manual(values = getPalette(length(unique(var_stat$Distribution)))) +
  scale_alpha_manual(values = c(0.6, 1), breaks = c("SNP", "INDEL")) +
  labs(x = NULL, y = "Number", fill = "Distribution") +
  theme_classic(base_line_size = 1) +
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0)
  )
ggsave(Pvar, filename = str_replace(basename(snpvarFilename), "SNPs\\.anno\\.variant_function", "variations.distribution.pdf"), width = 10, height = 6)
ggsave(Pvar, filename = str_replace(basename(snpvarFilename), "SNPs\\.anno\\.variant_function", "variations.distribution.png"), width = 10, height = 6, dpi = 500)

## 
snpex <- read_tsv(file = snpexFilename, col_names = F)
snpex_stat <- snpex %>% 
  group_by(X2) %>% 
  summarise(Number = n()) %>%
  rename(Type = X2) %>%
  mutate(Type = fct_reorder(Type, desc(Number)))
snpex_stat 
Psnpex <- ggplot(snpex_stat, aes(x = Type, y = Number)) +
  geom_bar(aes(fill = Type), stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Number of SNPs") +
  theme_classic(base_line_size = 1) +
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "NULL"
  )
ggsave(Psnpex, filename = str_replace(basename(snpexFilename), "$", ".pdf"), width = 4, height = 4.5)
ggsave(Psnpex, filename = str_replace(basename(snpexFilename), "$", ".png"), width = 4, height = 4.5, dpi = 500)


#
indelex <- read_tsv(file = indelexFilename, col_names = F)
indelex_stat <- indelex %>% 
  group_by(X2) %>% 
  summarise(Number = n()) %>%
  rename(Type = X2) %>%
  mutate(Type = fct_reorder(Type, desc(Number)))
indelex_stat 
Pindelex <- ggplot(indelex_stat, aes(x = Type, y = Number)) +
  geom_bar(aes(fill = Type), stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = NULL, y = "Number of INDELs") +
  theme_classic(base_line_size = 1) +
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "NULL"
  )
ggsave(Pindelex, filename = str_replace(basename(indelexFilename), "$", ".pdf"), width = 5, height = 4.5)
ggsave(Pindelex, filename = str_replace(basename(indelexFilename), "$", ".png"), width = 5, height = 4.5, dpi = 500)


## 输出表格
write_csv(x = snpvar_stat, file = str_replace(basename(snpvarFilename), "anno\\.variant_function", "distribution.csv"))
write_csv(x = indelvar_stat, file = str_replace(basename(indelvarFilename), "anno\\.variant_function", "distribution.csv"))
write_csv(x = snpex_stat, file = str_replace(basename(snpexFilename), "$", ".csv"))
write_csv(x = indelex_stat, file = str_replace(basename(indelexFilename), "$", ".csv"))

