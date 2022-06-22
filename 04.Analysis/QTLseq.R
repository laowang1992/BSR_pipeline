#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for QTLseq, TABEL file from GATK as input, 2 parent is needed")

## Add command line arguments
#
p <- add_argument(p, "--list", short = "-l", help = "list table header", flag = TRUE)
#
p <- add_argument(p, "--input", help = "TABEL file from GATK", type = "character")
p <- add_argument(p, "--out", help = "A prefix for output file", type = "character")
p <- add_argument(p, "--highP", help = "The parent name with high phenotype", type = "character")
p <- add_argument(p, "--lowP", help = "The parent name with low phenotype", type = "character")
p <- add_argument(p, "--highB", help = "The bulk name with high phenotype", type = "character")
p <- add_argument(p, "--lowB", help = "The bulk name with low phenotype", type = "character")
p <- add_argument(p, "--popType", help = "Population type, 'F2' or 'RIL'", type = "character")
#
p <- add_argument(p, "--minQ", help = "Minimum QUAL for variation", type = "numeric", default = 100)
#
p <- add_argument(p, "--bulkSizeH", help = "Bulk size with high phenotype", type = "numeric")
p <- add_argument(p, "--bulkSizeL", help = "Bulk size with low phenotype", type = "numeric")

p <- add_argument(p, "--minGQ", help = "Minimum GQ for parent genotype", type = "numeric", default = 30)
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
p <- add_argument(p, "--winSize", help = "Window size for sliding window statistics", type = "numeric", default = 1000000)
p <- add_argument(p, "--winStep", help = "Window step for sliding window statistics", type = "numeric", default = 200000)

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
write.table(argv, paste(argv$out, "QTLseqrparameter.txt", sep = "."), quote = F, sep = "\t", row.names = F)

library(tidyverse)
library(QTLseqr)
library(windowscanr)
library(cowplot)
library(CMplot)
library(ggsci)
library(RColorBrewer)

filename <- argv$input
outPrefix <- argv$out
minQ <- argv$minQ
bulkSizeH <- argv$bulkSizeH
bulkSizeL <- argv$bulkSizeL
highP <- argv$highP
lowP <- argv$lowP
highB <- argv$highB
lowB <- argv$lowB
popType <- argv$popType
minGQ <- argv$minGQ
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
  filename <- "./purpleLeaf.filter.SNPs.table.gz"
  outPrefix <- "purpleLeaf"
  ## filter SNP Quality
  minQ <- 50
  ##
  bulkSizeH <- 30
  bulkSizeL <- 30
  ## sample name
  highP <- "PP"
  lowP <- "PG"
  highB <- "F2P"
  lowB <- "F2G"
  ##
  popType <- "F2"  # F2 ro RIL
  ##
  minGQ <- 30
  ## filter depth parameter
  minHPdp <- 6
  maxHPdp <- 1000
  
  minLPdp <- 6
  maxLPdp <- 1000
  
  minHBdp <- 6
  maxHBdp <- 1000
  
  minLBdp <- 6
  maxLBdp <- 1000
  
  ## sliding window parameter
  winSize <- 1000000
  winStep <- 200000
  #
  minN <- 5
  ##
  width <- 15
  height <- 5
}

## 读取数据
chr <- read_tsv("./chrom.txt")
len <- read_tsv(file = "./ref.len", col_names = c("CHROM", "Len"))

df <- read_tsv(file = filename) %>%
  select(CHROM, POS, REF, ALT, 
         starts_with(all_of(highP)), 
         starts_with(all_of(lowP)), 
         starts_with(all_of(highB)), 
         starts_with(all_of(lowB))) %>% 
  na.omit() %>% rename(highParent.GT = paste(highP, "GT", sep = "."),
                       lowParent.GT = paste(lowP, "GT", sep = "."),
                       highBulk.GT = paste(highB, "GT", sep = "."),
                       lowBulk.GT = paste(lowB, "GT", sep = "."),
                       highParent.AD = paste(highP, "AD", sep = "."),
                       lowParent.AD = paste(lowP, "AD", sep = "."),
                       highBulk.AD = paste(highB, "AD", sep = "."),
                       lowBulk.AD = paste(lowB, "AD", sep = "."),
                       highParent.DP = paste(highP, "DP", sep = "."),
                       lowParent.DP = paste(lowP, "DP", sep = "."),
                       highBulk.DP = paste(highB, "DP", sep = "."),
                       lowBulk.DP = paste(lowB, "DP", sep = "."),
                       highParent.GQ = paste(highP, "GQ", sep = "."),
                       lowParent.GQ = paste(lowP, "GQ", sep = "."),
                       highBulk.GQ = paste(highB, "GQ", sep = "."),
                       lowBulk.GQ = paste(lowB, "GQ", sep = "."),
                       highParent.PL = paste(highP, "PL", sep = "."),
                       lowParent.PL = paste(lowP, "PL", sep = "."),
                       highBulk.PL = paste(highB, "PL", sep = "."),
                       lowBulk.PL = paste(lowB, "PL", sep = "."))

dd1 <- df %>% filter(highParent.GQ >= minGQ, lowParent.GQ >= minGQ) %>%
  filter(highParent.GT == paste(REF, REF, sep = "/") | highParent.GT == paste(ALT, ALT, sep = "/")) %>%
  filter(lowParent.GT == paste(REF, REF, sep = "/") | lowParent.GT == paste(ALT, ALT, sep = "/")) %>%
  filter(highParent.GT != lowParent.GT)

dd2 <- dd1 %>% separate(highBulk.AD, c("highBulk.AD_0", "highBulk.AD_1"), sep = ",", convert = TRUE) %>%
  separate(lowBulk.AD, c("lowBulk.AD_0", "lowBulk.AD_1"), sep = ",", convert = TRUE) %>%
  separate(highBulk.PL, c("highBulk.PL_00", "highBulk.PL_01", "highBulk.PL_11"), sep = ",", convert = TRUE) %>%
  separate(lowBulk.PL, c("lowBulk.PL_00", "lowBulk.PL_01", "lowBulk.PL_11"), sep = ",", convert = TRUE)

## 
dd3 <- dd2 %>% filter(!(((highBulk.AD_0 / (highBulk.AD_0 + highBulk.AD_1) < 0.3) & 
                          (lowBulk.AD_0 / (lowBulk.AD_0 + lowBulk.AD_1) < 0.3)) | 
                        ((highBulk.AD_0 / (highBulk.AD_0 + highBulk.AD_1) > 0.7) & 
                           (lowBulk.AD_0 / (lowBulk.AD_0 + lowBulk.AD_1) > 0.7))))
## 
dp <- dd3 %>% dplyr::select(HP = highParent.DP, LP = lowParent.DP, HB = highBulk.DP, LB = lowBulk.DP) %>%
  gather(key = "sample", value = "depth")
dp$sample <- factor(dp$sample, levels = c("HP", "LP", "HB", "LB"), labels = c(highP, lowP, highB, lowB))
P_dp <- ggplot(dp, aes(x = depth)) + 
  geom_histogram(aes(y = ..density.., fill = sample), binwidth = 2) +
  geom_density() + 
  scale_x_continuous(limits = c(0, 100)) +
  theme_half_open() +
  facet_wrap(~sample, nrow = 1)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.pdf", sep = "."), height = 3.5, width = 10)
ggsave(P_dp, filename = paste(outPrefix, "depth_desity.png", sep = "."), height = 3.5, width = 10, dpi = 500)
remove(P_dp)

## 根据深度过滤
dd4 <- dd3 %>% filter(highParent.DP > minHPdp, highParent.DP < maxHPdp,
         lowParent.DP > minLPdp, lowParent.DP < maxLPdp,
         highBulk.DP > minHBdp, highBulk.DP < maxHBdp,
         lowBulk.DP > minLBdp, lowBulk.DP < maxLBdp)


##
SNPnumber <- dd4 %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))
##
options(scipen = 200)
colourCount = dim(chr)[[1]]
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
#df$LABEL <- factor(df$LABEL, levels = chromColor$LABEL)
chr$LABEL <- factor(chr$LABEL, levels = chr$LABEL)
Phist <- chr %>% left_join(dd4, by = "CHROM") %>% ggplot(aes(x = POS)) +
  geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +
  labs(x = NULL, y = "SNP Count / 1Mb", fill = "Chrom") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_half_open() +
  theme(strip.text = element_text(color = NA, size = 0.1),
        strip.background = element_rect(color = NA, fill = NA)) +
  facet_grid(LABEL ~ .)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.pdf", sep = "."), width = 9, height = dim(chr)[[1]] * 0.6 + 0.5)
ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.png", sep = "."), width = 9, height = dim(chr)[[1]] * 0.6 + 0.5, dpi = 500)
remove(Phist)

## 导出用于QTLseqr的数据
dd5  <- dd4 %>% mutate(highBulk.AD = if_else(highParent.GT == paste(ALT, ALT, sep = "/"), 
                                             paste(highBulk.AD_0, highBulk.AD_1, sep = ","), 
                                             paste(highBulk.AD_1, highBulk.AD_0, sep = ",")),
                       lowBulk.AD = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
                                            paste(lowBulk.AD_0, lowBulk.AD_1, sep = ","),
                                            paste(lowBulk.AD_1, lowBulk.AD_0, sep = ",")),
                       highBulk.PL = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
                                             paste(highBulk.PL_00, highBulk.PL_01, highBulk.PL_11, sep = ","),
                                             paste(highBulk.PL_11, highBulk.PL_01, highBulk.PL_00, sep = ",")),
                       lowBulk.PL = if_else(highParent.GT == paste(ALT, ALT, sep = "/"),
                                            paste(lowBulk.PL_00, lowBulk.PL_01, lowBulk.PL_11, sep = ","),
                                            paste(lowBulk.PL_11, lowBulk.PL_01, lowBulk.PL_00, sep = ","))) %>%
  select(CHROM, POS, REF, ALT, 
         highBulk.AD, highBulk.DP, highBulk.GQ, highBulk.PL,
         lowBulk.AD, lowBulk.DP, lowBulk.GQ, lowBulk.PL)
write_tsv(dd5, paste(outPrefix, "reform.txt", sep = "."))

##
dd6 <- dd4 %>% mutate(HB.HP.DP = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         highBulk.AD_0, highBulk.AD_1),
                      HB.LP.DP = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         highBulk.AD_1, highBulk.AD_0),
                      LB.HP.DP = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         lowBulk.AD_0, lowBulk.AD_1),
                      LB.LP.DP = if_else(highParent.GT == paste(REF, REF, sep = "/"),
                                         lowBulk.AD_1, lowBulk.AD_0)) %>%
  dplyr::select(CHROM, POS, HB.HP.DP, HB.LP.DP, LB.HP.DP, LB.LP.DP)
## Sliding window statistical
slidwin <- winScan(x = dd6,
              groups = "CHROM",
              position = "POS",
              values = c("HB.HP.DP", "HB.LP.DP", "LB.HP.DP", "LB.LP.DP"),
              win_size = winSize,
              win_step = winStep,
              funs = c("sum")) %>% as_tibble() %>% 
  mutate(HB.index = HB.HP.DP_sum / (HB.HP.DP_sum + HB.LP.DP_sum),
         LB.index = LB.HP.DP_sum / (LB.HP.DP_sum + LB.LP.DP_sum),
         delta.index = HB.index - LB.index,
         ED = sqrt((HB.index - LB.index)^2 + ((1-HB.index) - (1-LB.index))^2),
         ED4 = ED^4) %>% 
  dplyr::select(CHROM, win_start, win_end, win_mid, SNPn = LB.LP.DP_n, 
                       HB.HP.DP_sum, HB.LP.DP_sum, HB.index, 
                       LB.HP.DP_sum, LB.LP.DP_sum, LB.index, 
                       delta.index, ED, ED4)
outtb <- slidwin
colnames(outtb) <- c("CHROM", "win_start", "win_end", "win_mid", "SNPn",
                     paste(highB, highP, "depth", sep = "."),
                     paste(highB, lowP, "depth", sep = "."),
                     paste(highB, "index", sep = "."),
                     paste(lowB, highB, "depth", sep = "."),
                     paste(lowB, lowP, "depth", sep = "."),
                     paste(lowB, "index", sep = "."),
                     "delta.index", "ED", "ED4")
write_tsv(x = outtb, file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
write_csv(x = outtb, file = paste(outPrefix, "SlidingWindow.csv", sep = "."))



## 我们开始画图了
slidwin <- chr %>% left_join(slidwin, by = "CHROM")
slidwin$CHROM <- factor(slidwin$CHROM, levels = chr$CHROM)

## SNP index
# delta index
pdf(file = paste(outPrefix, "delta_SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, delta_index = delta.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "delta SNP index", ylim = c(-1, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, "delta_SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, delta.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "delta SNP index", ylim = c(-1, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

# high bulk index
pdf(file = paste(outPrefix, highB, "SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, highB, "SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, HB.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

# low bulk index
pdf(file = paste(outPrefix, lowB, "SNP_index.pdf", sep = "."), width = width, height = height)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, LB.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, lowB, "SNP_index.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, LB.index), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "SNP index", ylim = c(0, 1), cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()

## ED and ED^4
# ED
pdf(file = paste(outPrefix, "ED.pdf", sep = "."), width = width, height = height)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, ED), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "ED", cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, "ED.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, ED), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "ED", cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
# ED4
pdf(file = paste(outPrefix, "ED4.pdf", sep = "."), width = width, height = height)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, ED4), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "ED^4", cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()
png(filename = paste(outPrefix, "ED4.png", sep = "."), width = width, height = height, units = "in", res = 500)
CMplot(filter(slidwin, SNPn > minN) %>% select(SNP = CHROM, Chromosome = CHROM, Postion = win_mid, ED4), 
       type = "p", plot.type = c("m"), band = 0.5, LOG10 = FALSE, chr.labels = chr$LABEL,
       ylab = "ED^4", cex = 0.5, signal.cex = 0.8, chr.labels.angle = 45,
       chr.den.col = NULL, ylab.pos = 2.7, amplify = FALSE,
       file.output = FALSE)
dev.off()


#### 
df <- importFromGATK(file = paste(outPrefix, "reform.txt", sep = "."),
                     highBulk = "highBulk",
                     lowBulk = "lowBulk",
                     #chr=c('C1','C2',C3)
                     )
df <- chr %>% left_join(df, by = "CHROM")

#df <- runGprimeAnalysis(
#  SNPset = df,
#  windowSize = winSize,    # window size
#  outlierFilter = "deltaSNP")

df <- runQTLseqAnalysis(
  SNPset = df,
  windowSize = winSize, #window size
  popStruc = popType,  #F2 or RIL
  bulkSize = c(bulkSizeH, bulkSizeL),  #bulk size, first is highbulk
  replications = 1000, 
  intervals = c(95, 99) # confidence interval
)

# 输出QTLseqr结算结果
outtb <- df %>% select(CHROM, POS, REF, ALT, LowBulk.LPgeno.AD = AD_REF.LOW, LowBulk.HPgeno.AD = AD_ALT.LOW, 
                       SNPindex.LOW, HighBulk.LPgeno.AD = AD_REF.HIGH, HighBulk.HPgeno.AD = AD_ALT.HIGH, SNPindex.HIGH, 
                       deltaSNP, nSNPs, tricubeDeltaSNP, CI_95, CI_99)
write_tsv(outtb, paste(outPrefix, "QTLseqrdeltaSNPindex.txt", sep = "."))
write_csv(outtb, paste(outPrefix, "QTLseqrdeltaSNPindex.csv", sep = "."))

# 
df <- df %>% filter(nSNPs > minN)

#plotQTLStats(SNPset = df, var = "Gprime", plotThreshold = TRUE, q = 0.01)
#plotQTLStats(SNPset = df, var = "deltaSNP", plotIntervals = TRUE) + theme_half_open()
# add chrom color infomation
source("./plot.R")
#len <- chr %>% left_join(chr, by = "CHROM")
df$LABEL <- factor(df$LABEL, levels = chr$LABEL)
##
pdf(file = paste(outPrefix, "deltaSNPindex.95CI.pdf", sep = "."), width = width, height = height)
plotIndex(df = outtb, chr = chr, len = len, CI = "CI_95", nSNPs = "nSNPs", n = minN, 
          index = "tricubeDeltaSNP", band = 0.005, color = c("#177cb0","#f36838"), 
          ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()
#
png(filename = paste(outPrefix, "deltaSNPindex.95CI.png", sep = "."), width = width, height = height, units = "in", res = 500)
plotIndex(df = outtb, chr = chr, len = len, CI = "CI_95", nSNPs = "nSNPs", n = minN, 
          index = "tricubeDeltaSNP", band = 0.005, color = c("#177cb0","#f36838"), 
          ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()

##
pdf(file = paste(outPrefix, "deltaSNPindex.99CI.pdf", sep = "."), width = width, height = height)
plotIndex(df = outtb, chr = chr, len = len, CI = "CI_99", nSNPs = "nSNPs", n = minN, 
          index = "tricubeDeltaSNP", band = 0.005, color = c("#177cb0","#f36838"), 
          ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()

png(filename = paste(outPrefix, "deltaSNPindex.99CI.png", sep = "."), width = width, height = height, units = "in", res = 500)
plotIndex(df = outtb, chr = chr, len = len, CI = "CI_99", nSNPs = "nSNPs", n = minN, 
          index = "tricubeDeltaSNP", band = 0.005, color = c("#177cb0","#f36838"), 
          ylim = c(-1, 1), ylab = "Delta SNP index", size = 1.5, axis.size = 1, 
          axis.lab.size = 1, axis.title.size = 1.5, type = "l")
dev.off()


getQTLTable(SNPset = df, interval =95, export = TRUE, fileName = paste(outPrefix, "95CI.csv", sep = "."), method="QTLseq")
getQTLTable(SNPset = df, interval =99, export = TRUE, fileName = paste(outPrefix, "99CI.csv", sep = "."), method="QTLseq")


## plot target chrom
options(scipen=200)
plotTargetChrom(df = df, CI = 95, minN = minN, outPrefix = outPrefix, chr = chr, len = len)
plotTargetChrom(df = df, CI = 99, minN = minN, outPrefix = outPrefix, chr = chr, len = len)

