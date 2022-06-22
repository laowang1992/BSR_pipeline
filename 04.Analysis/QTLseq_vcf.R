#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for QTLseq, a vcf file as input, 2 parent is needed")

## Add command line arguments
#
p <- add_argument(p, "--list", short = "-l", help = "list table header", flag = TRUE)
#
p <- add_argument(p, "--input", help = "input vcf file for QTLseq", type = "character")
p <- add_argument(p, "--out", help = "A prefix for output file", type = "character")
p <- add_argument(p, "--highP", help = "The parent name with high phenotype", type = "character")
p <- add_argument(p, "--lowP", help = "The parent name with low phenotype", type = "character")
p <- add_argument(p, "--highB", help = "The bulk name with high phenotype", type = "character")
p <- add_argument(p, "--lowB", help = "The bulk name with low phenotype", type = "character")
#p <- add_argument(p, "--popType", help = "Population type, 'F2' or 'RIL'", type = "character")
#
p <- add_argument(p, "--minQ", help = "Minimum QUAL for variation", type = "numeric", default = 100)
#
p <- add_argument(p, "--minHPdp", help = "Minimum depth for high parent", type = "numeric", default = 10)
p <- add_argument(p, "--maxHPdp", help = "Maxmum depth for high parent", type = "numeric", default = 10000)
p <- add_argument(p, "--minLPdp", help = "Minimum depth for low parent", type = "numeric", default = 10)
p <- add_argument(p, "--maxLPdp", help = "Maxmum depth for low parent", type = "numeric", default = 10000)
p <- add_argument(p, "--minHBdp", help = "Minimum depth for high bulk", type = "numeric", default = 10)
p <- add_argument(p, "--maxHBdp", help = "Maxmum depth for high bulk", type = "numeric", default = 10000)
p <- add_argument(p, "--minLBdp", help = "Minimum depth for low bulk", type = "numeric", default = 10)
p <- add_argument(p, "--maxLBdp", help = "Maxmum depth for low bulk", type = "numeric", default = 10000)
#
p <- add_argument(p, "--winSize", help = "Window size for sliding window statistics", type = "numeric", default = 500000)
p <- add_argument(p, "--winStep", help = "Step width for sliding window statistics", type = "numeric", default = 250000)
#
p <- add_argument(p, "--minN", help = "Minimum SNP number in a window for plot", type = "numeric", default = 10)
#
p <- add_argument(p, "--width", help = "Delta SNP index plot width", type = "numeric", default = 15)
p <- add_argument(p, "--height", help = "Delta SNP index plot height", type = "numeric", default = 5)

# Parse the command line arguments
argv <- parse_args(p)
#
if (argv$list) {
  df <- read.table(argv$input, nrows = 1, header = T)
  col_name <- colnames(df)
  print(col_name)
  quit(save = "no")
}

filename <- argv$input
outPrefix <- argv$out
minQ <- argv$minQ
highP <- argv$highP
lowP <- argv$lowP
highB <- argv$highB
lowB <- argv$lowB
minHPdp <- argv$minHPdp
maxHPdp <- argv$maxHPdp
minLPdp <- argv$minLPdp
maxLPdp <- argv$maxLPdp
minHBdp <- argv$minHBdp
maxHBdp <- argv$maxHBdp
minLBdp <- argv$minLBdp
maxLBdp <- argv$maxLBdp
winSize <- argv$winSize
winStep <- argv$winStep
minN <- argv$minN
width <- argv$width
height <- argv$height

if (FALSE) {
  ## file name
  filename <- "./purpleLeaf.filter.SNPs.txt.gz"
  outPrefix <- "purpleLeaf"
  ## filter SNP Quality
  minQ <- 50
  
  ## sample name
  # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s8417	FM195BML	M1F	N1F
  highP <- "PP"
  lowP <- "PG"
  highB <- "F2P"
  lowB <- "F2G"
  
  ## filter depth parameter
  minHPdp <- 6
  maxHPdp <- 60
  
  minLPdp <- 6
  maxLPdp <- 60
  
  minHBdp <- 6
  maxHBdp <- 60
  
  minLBdp <- 6
  maxLBdp <- 60
  
  ## sliding window parameter
  winSize <-1000000
  winStep <- 100000
  
  ## plot parameter
  minN <- 10
  width <-15
  height <- 4
}

library(tidyverse)
library(cowplot)
library(ggsci)
library(windowscanr)
library(RColorBrewer)
library(CMplot)

df <- read_tsv(file = filename)
df <- df %>% dplyr::select(CHROM, POS, REF, ALT, QUAL, HP = all_of(highP), LP = all_of(lowP), HB = all_of(highB), LB = all_of(lowB))

# 先过滤亲本基因型，纯合且不同
df <- df %>% filter((str_detect(HP, "^0/0") & str_detect(LP, "^1/1")) | (str_detect(LP, "^0/0") & str_detect(HP, "1/1"))) %>%
  filter(!(str_detect(HB, "^\\./\\.") | str_detect(LB, "^\\./\\."))) %>%
  filter(ALT == "A" | ALT == "C" | ALT == "G" | ALT == "T")

chromColor <- read_tsv("./chrom.txt")
dd <- df %>% filter(QUAL > minQ) #过滤QUAL小于等于50的行
l1<-nchar(dd$REF)	#求REF字符串长度
l2<-nchar(dd$ALT)	#求ALT字符串长度

#过滤indel行，保留SNP信息
lse <- l1 == 1 & l2 == 1
dd <- dd[lse, ]

dd1 <- dd %>% separate(HP, c("HP.geno", "HP"), sep = ":") %>%
  separate(LP, c("LP.geno", "LP"), sep = ":") %>%
  separate(HB, c("HB.geno", "HB"), sep = ":") %>%
  separate(LB, c("LB.geno", "LB"), sep = ":")

#得到各sample每种基因型的覆盖度
dd2 <- dd1 %>% 
  separate(HP, c("HP.ref.dp", "HP.alt.dp"), sep = ",", convert = T) %>%
  separate(LP, c("LP.ref.dp", "LP.alt.dp"), sep = ",", convert = T) %>%
  separate(HB, c("HB.ref.dp", "HB.alt.dp"), sep = ",", convert = T) %>%
  separate(LB, c("LB.ref.dp", "LB.alt.dp"), sep = ",", convert = T)

# 计算每个样本总深度
hh <- dd2 %>% 
  mutate(HP.sum = HP.ref.dp + HP.alt.dp,
         LP.sum = LP.ref.dp + LP.alt.dp,
         HB.sum = HB.ref.dp + HB.alt.dp,
         LB.sum = LB.ref.dp + LB.alt.dp)

hh <- hh %>% 
  mutate(HB.HP.dp = if_else(HP.geno == "1/1", HB.alt.dp, HB.ref.dp),
         HB.LP.dp = if_else(HP.geno == "1/1", HB.ref.dp, HB.alt.dp),
         LB.HP.dp = if_else(HP.geno == "1/1", LB.alt.dp, LB.ref.dp),
         LB.LP.dp = if_else(HP.geno == "1/1", LB.ref.dp, LB.alt.dp))

table(paste(hh$HP.geno, hh$LP.geno, sep = "_"))

# 过滤index都大于0.7或都小于0.3的位点，这一步对计算每个样本的index结果影响很显著，对Δindex并不显著
hh <- hh %>% filter(!((HB.alt.dp/HB.sum > 0.7 & LB.alt.dp/LB.sum > 0.7) | (HB.alt.dp/HB.sum < 0.3 & LB.alt.dp/LB.sum < 0.3)))


dp <- hh %>% dplyr::select(HP = HP.sum, LP = LP.sum, HB = HB.sum, LB = LB.sum) %>%
  gather(key = "sample", value = "depth")
dp$sample <- factor(dp$sample, levels = c("HP", "LP", "HB", "LB"), labels = c(highP, lowP, highB, lowB))
P_dp <- ggplot(dp, aes(x = depth)) + 
  geom_histogram(aes(y = ..density.., fill = sample), binwidth = 2) +
  geom_density() + 
  theme_half_open() +
  facet_wrap(~sample, nrow = 1)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.pdf", sep = "."), height = 3.5, width = 10)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.png", sep = "."), height = 3.5, width = 10, dpi = 500)
remove(P_dp)

hh <- hh %>%
  filter(HP.sum > minHPdp, HP.sum < maxHPdp,
         LP.sum > minLPdp, LP.sum < maxLPdp,
         HB.sum > minHBdp, HB.sum < maxHBdp,
         LB.sum > minLBdp, LB.sum < maxLBdp
         )

##
SNPnumber <- hh %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))
##
options(scipen = 200)
colourCount = dim(chromColor)[[1]]
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
#df$LABEL <- factor(df$LABEL, levels = chromColor$LABEL)
chromColor$LABEL <- factor(chromColor$LABEL, levels = chromColor$LABEL)
Phist <- chromColor %>% left_join(hh, by = "CHROM") %>% ggplot(aes(x = POS)) +
  geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +
  labs(x = NULL, y = "SNP Count / 1Mb") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_half_open() +
  theme(strip.text = element_text(color = NA, size = 0.1),
        strip.background = element_rect(color = NA, fill = NA)) +
  facet_grid(LABEL ~ .)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.pdf", sep = "."), width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.png", sep = "."), width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5, dpi = 500)
remove(Phist)

#
hh$chr <- hh$CHROM
w <- hh %>% dplyr::select(CHROM, POS, HB.HP.dp, HB.LP.dp)
b <- hh %>% dplyr::select(CHROM, POS, LB.HP.dp, LB.LP.dp)

#############################
#ww$R.p1.depth<-ww$R.p1.depth+0.5
#bb$S.p1.depth<-bb$S.p1.depth+0.5
w1 <- winScan(x = hh,
              groups = "CHROM",
              position = "POS",
              values = c("HB.HP.dp", "HB.LP.dp"),
              win_size = winSize,
              win_step = winStep,
              funs = c("sum"))
b1 <- winScan(x = hh,
              groups = "CHROM",
              position = "POS",
              values = c("LB.HP.dp", "LB.LP.dp"),
              win_size = winSize,
              win_step = winStep,
              funs = c("sum"))
#w1$w<-w1$R.p2.depth_sum/(w1$R.p1.depth_sum)
#b1$b<-b1$S.p2.depth_sum/(b1$S.p1.depth_sum)
#w1$w<-w1$R.p2.depth_sum/(w1$R.p1.depth_sum+w1$R.p2.depth_sum)
#b1$b<-b1$S.p2.depth_sum/(b1$S.p1.depth_sum+b1$S.p2.depth_sum)
###
w1 <- w1 %>% mutate(index = HB.HP.dp_sum / (HB.HP.dp_sum + HB.LP.dp_sum))
b1 <- b1 %>% mutate(index = LB.HP.dp_sum / (LB.HP.dp_sum + LB.LP.dp_sum))
#w1$w <- w1[, paste(highB, highP, "depth_sum", sep = ".")] / (w1[, paste(highB, highP, "depth_sum", sep = ".")] + w1[, paste(highB, lowP, "depth_sum", sep = ".")])
#b1$b <- b1[, paste(lowB, highP, "depth_sum", sep = ".")] / (b1[, paste(lowB, highP, "depth_sum", sep = ".")] + b1[, paste(lowB, lowP, "depth_sum", sep = ".")])
#plot(w1$w - b1$b,col=factor(w1$CHROM),pch=20,cex=0.5,xaxt="n",ylab="SNP reads ratio")
d <- tibble(CHROM = w1$CHROM, w1$win_start, w1$win_mid, w1$win_end, 
            SNP_N = w1$HB.HP.dp_n, INDEX_W = w1$index, INDEX_B = b1$index,
            delta_INDEX = w1$index - b1$index)

colnames(d) <- c("CHROM", "win_start", "win_mid", "win_end", "SNP_N", 
                 paste(highB, "SNP_INDEX", sep = "."), 
                 paste(lowB, "SNP_INDEX", sep = "."), "delta_SNP_INDEX")

write_tsv(x = d, file = paste(outPrefix, "SlidingWindow_SNPindex.txt", sep = "."))
write_csv(x = d, file = paste(outPrefix, "SlidingWindow_SNPindex.csv", sep = "."))

# 加上绘图颜色信息
d <- chromColor %>% left_join(d, by = "CHROM")

#mytheme<-theme(axis.title = element_text(face = "bold",size = 10),
#               axis.text = element_text(face = "bold",size = 9),
#               panel.background = element_rect(fill = "white",color = "darkblue"),
#               panel.grid.major.y = element_line(color = "grey",linetype = 1),
#               panel.grid.minor.y = element_line(color = "grey",linetype = 2),
#               panel.grid.minor.x = element_blank(),
#               legend.position = "none")
d$LABEL <- factor(d$LABEL, levels = chromColor$LABEL)
#Pdelta <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = delta_SNP_INDEX)) +
#  geom_point(aes(color = COLOR), size = 0.7) +
#  ylim(-1, 1) +
#  labs(x="", y="delta SNP index") +
#  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
#  scale_color_aaas() +
#  theme_minimal_grid() +
#  theme(legend.position = "none") +
#  #mytheme +
#  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#ggsave(Pdelta, filename = paste(outPrefix, "delta_SNP_index.pdf", sep = "."), width = width, height = height)
#ggsave(Pdelta, filename = paste(outPrefix, "delta_SNP_index.png", sep = "."), width = width, height = height, dpi = 500)
#plot(b1$b/w1$w,pch=20,xaxt="n",ylab="SNP reads ratio")

#Phigh <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = d[d$SNP_N > minN, paste(highB, "SNP_INDEX", sep = ".")][[1]])) +
#  geom_point(aes(color = COLOR), size = 0.7) +
#  ylim(0, 1) +
#  labs(x="", y="SNP index") +
#  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
#  scale_color_aaas() +
#  theme_minimal_grid() +
#  theme(legend.position = "none") +
#  #mytheme +
#  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#ggsave(Phigh, filename = paste(outPrefix, highB, "SNP_index.pdf", sep = "."), width = width, height = height)
#ggsave(Phigh, filename = paste(outPrefix, highB, "SNP_index.png", sep = "."), width = width, height = height, dpi = 500)

#Plow <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = d[d$SNP_N > minN, paste(lowB, "SNP_INDEX", sep = ".")][[1]])) +
#  geom_point(aes(color = COLOR), size = 0.7) +
#  ylim(0, 1) +
#  labs(x="", y="SNP index") +
#  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
#  scale_color_aaas() +
#  theme_minimal_grid() +
#  theme(legend.position = "none") +
#  #mytheme +
#  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#ggsave(Plow, filename = paste(outPrefix, lowB, "SNP_index.pdf", sep = "."), width = width, height = height)
#ggsave(Plow, filename = paste(outPrefix, lowB, "SNP_index.png", sep = "."), width = width, height = height, dpi = 500)

#
pdf(file = paste(outPrefix, "delta_SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, delta_index = delta_SNP_INDEX), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "delta SNP index", ylim = c(-1, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, "delta_SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, delta_index = delta_SNP_INDEX), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "delta SNP index", ylim = c(-1, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

#
pdf(file = paste(outPrefix, highB, "SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB_index = paste(highB, "SNP_INDEX", sep = ".")), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, highB, "SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB_index = paste(highB, "SNP_INDEX", sep = ".")), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

#
pdf(file = paste(outPrefix, lowB, "SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB_index = paste(lowB, "SNP_INDEX", sep = ".")), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, lowB, "SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(d, SNP_N > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB_index = paste(lowB, "SNP_INDEX", sep = ".")), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chromColor$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

