# 加载扩展???
library(tidyverse)
library(ggplot2)
library(ggsci)
library(cowplot)
library(reshape2)
# 读取数据并转换格???
df <- read.table(file="genes.TMM.EXPR.matrix", header=T)
head(df)
d <- melt(df)

mytheme<-theme(axis.title = element_text(face = "bold",size = 10),
               axis.text = element_text(face = "bold",size = 9),
               panel.background = element_rect(fill = "white",color = "darkblue"),
               panel.grid.major.y = element_line(color = "grey",linetype = 1),
               panel.grid.minor.y = element_line(color = "grey",linetype = 2),
               panel.grid.minor.x = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 1))

p1 <- ggplot(data=d, aes(x=variable, y=log10(value+1), fill=variable)) +
  geom_violin(size=.5) +
  geom_boxplot(width=.3, size=.5) +
  xlab("Samples") +
  ylab("log10(TMM+1)") +
  labs(fill="Sample") +
  #scale_fill_aaas() +
  mytheme
ggsave(p1, filename="genes.TMM.EXPR.pdf", height = 5, width = 13)
ggsave(p1, filename="genes.TMM.EXPR.png", height = 5, width = 13, dpi = 500)

p2 <-ggplot(data=d, aes(x=log10(value+1), color=variable)) +
  geom_density(size=1) +
  xlab("log10(TMM+1)") +
  ylab("Density")+
  labs(fill="Sample") +
  mytheme
ggsave(p2, filename="genes.TMM.EXPR.density.pdf", height = 4, width = 8.5)
ggsave(p2, filename="genes.TMM.EXPR.density.png", height = 4, width = 8.5, dpi = 500)


gene_exp <- read.table("genes.TMM.EXPR.matrix", 
                       header = T, row.names=1)
metadata <- read.table("./metadata.txt", 
                       header = T, row.names = 1)
head(metadata)

# PCA
library(PCAtools)
p <- pca(log10(gene_exp+1), metadata = metadata)
pdf(file = "screeplot.pdf", width = 10, height = 5)
screeplot(p)
dev.off()
png(file = "screeplot.png", width = 10, height = 5, units = "in", res = 500)
screeplot(p)
dev.off()
#biplot(p, 
#       x = 'PC1',                 # x ???
#       y = 'PC2',                 # y ???
#       colby = 'Strain',          # 颜色映射
#       shape = 'Stage',
#       legendPosition = 'right',  # 图例位置
#       lab = rownames(metadata)                    # 样本名称显示
#)

pca_rotated_plus <- rownames_to_column(p$rotated, 
                                       var = 'sample_name') %>%
  left_join(rownames_to_column(metadata, var = 'sample_name'), 
            by = 'sample_name')


p_pca <- ggplot(data = pca_rotated_plus, aes(x = PC1, y = PC2)) +
  geom_point(size = 6, 
             aes(shape = Strain, fill = Stage)) +
  #stat_ellipse(aes(color = Stage)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_aaas() +
  scale_color_npg() +
  labs(x = paste("PC1 (", round(p[["variance"]][["PC1"]], 2), "% variance explained)", sep = ""),
       y = paste("PC1 (", round(p[["variance"]][["PC2"]], 2), "% variance explained)", sep = "")) +
  theme_half_open() + 
  theme(#legend.position = c(0.2, 0.9), 
    legend.background = element_rect(fill = NA)
  ) +
  guides(fill = guide_legend(override.aes=list(shape=21)))

ggsave("PCAplot_PC1&PC2.pdf", height = 5, width = 7)
ggsave("PCAplot_PC1&PC2.png", height = 5, width = 7, dpi = 500)

# 相关???
sample_cor <- cor(log10(gene_exp+1), 
                  method = "pearson") # 算法??? pearson | kendall | spearman
write.table(x = as.table(sample_cor), file = "./correlation.txt")
write.csv(x = as.table(sample_cor), file = "./correlation.csv")

#annotation_colors = list( treatment = c(ck = "red", PW = "orange"))

pdf("correlation_withClust.pdf",width = 12,height = 11)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = T,       # 是否对行（基因）聚类
                   cluster_cols = T,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   #cutree_col = 3,            # 热图按列分为3???
                   #cutree_rows = 2,           # 热图按行分为2???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   #annotation_colors = annotation_colors,
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = data.frame(strain=metadata$Strain,row.names = rownames(metadata)), # 添加列注释，col_metadata.csv中信???
                   annotation_row = data.frame(stage=metadata$Stage,row.names = rownames(metadata)), # 添加行注释，row_metadata.csv中信???
                   cellwidth = 12,
                   cellheight = 12,
                   color = colorRampPalette(c("green","white","red"))(100),
                   #annotation_legend = F,
                   #legend = F
                   #main = "example heatmap"   # 添加标题
)
dev.off()
png("correlation_withClust.png", width = 12, height = 11, units = "in", res = 500)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = T,       # 是否对行（基因）聚类
                   cluster_cols = T,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   #cutree_col = 3,            # 热图按列分为3???
                   #cutree_rows = 2,           # 热图按行分为2???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   #annotation_colors = annotation_colors,
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = data.frame(strain=metadata$Strain,row.names = rownames(metadata)), # 添加列注释，col_metadata.csv中信???
                   annotation_row = data.frame(stage=metadata$Stage,row.names = rownames(metadata)), # 添加行注释，row_metadata.csv中信???
                   cellwidth = 12,
                   cellheight = 12,
                   color = colorRampPalette(c("green","white","red"))(100),
                   #annotation_legend = F,
                   #legend = F
                   #main = "example heatmap"   # 添加标题
)
dev.off()


pdf("correlation_withoutClust.pdf",width = 12,height = 11)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = F,       # 是否对行（基因）聚类
                   cluster_cols = F,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   #cutree_col = 3,            # 热图按列分为3???
                   #cutree_rows = 2,           # 热图按行分为2???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   #annotation_colors = annotation_colors,
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = data.frame(strain=metadata$Strain,row.names = rownames(metadata)), # 添加列注释，col_metadata.csv中信???
                   annotation_row = data.frame(stage=metadata$Stage,row.names = rownames(metadata)), # 添加行注释，row_metadata.csv中信???
                   cellwidth = 12,
                   cellheight = 12,
                   color = colorRampPalette(c("green","white","red"))(100),
                   #annotation_legend = F,
                   #legend = F
                   #main = "example heatmap"   # 添加标题
)
dev.off()
png("correlation_withoutClust.png", width = 12, height = 11, units = "in", res = 500)
pheatmap::pheatmap(sample_cor,                # 绘图数据，应进行log转换，log2或log10，并且防止表达量???0应该df+1
                   show_rownames = TRUE,      # 是否显示行名（基因名），基因太多时选FALSE
                   show_colnames = TRUE,      # 是否显示列明（样本名），同上
                   cluster_rows = F,       # 是否对行（基因）聚类
                   cluster_cols = F,       # 是否对列（样本）聚类
                   border_color = "NA",       # 网格分割线颜???
                   #cutree_col = 3,            # 热图按列分为3???
                   #cutree_rows = 2,           # 热图按行分为2???
                   display_numbers = FALSE,   # 是否显示数值，数值为上面转换过的大小
                   number_format = "%.2f",    # 数值保留小数点???2位数，当display_numbers=TRUE有效
                   #annotation_colors = annotation_colors,
                   scale = "none",            # ???"row", "column"???"none"，是否要进行标准化，如果要突出所有基因在处理组和对照组之间表达有差异，建议改???"row"，如果要突出在所有样本中基因的分组情况应改为"column"
                   annotation_col = data.frame(strain=metadata$Strain,row.names = rownames(metadata)), # 添加列注释，col_metadata.csv中信???
                   annotation_row = data.frame(stage=metadata$Stage,row.names = rownames(metadata)), # 添加行注释，row_metadata.csv中信???
                   cellwidth = 12,
                   cellheight = 12,
                   color = colorRampPalette(c("green","white","red"))(100),
                   #annotation_legend = F,
                   #legend = F
                   #main = "example heatmap"   # 添加标题
)
dev.off()

