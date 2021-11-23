#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for BSA, TABEL file from GATK as input")

## Add command line arguments
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
p <- add_argument(p, "--winSize", help = "Window size for sliding window statistics", type = "numeric", default = 2000000)
#
p <- add_argument(p, "--minN", help = "Minimum SNP number in a window for plot", type = "numeric", default = 10)
#
p <- add_argument(p, "--width", help = "Delta SNP index plot width", type = "numeric", default = 15)
p <- add_argument(p, "--height", help = "Delta SNP index plot height", type = "numeric", default = 5)

# Parse the command line arguments
argv <- parse_args(p)
write.table(argv, paste(argv$out, "QTLseqrparameter.txt", sep = "."), quote = F, sep = "\t", row.names = F)

library(tidyverse)
library(QTLseqr)
library(cowplot)
library(ggsci)

filename <- argv$input
outPrefix <- argv$out
minQ <- argv$minQ
bulkSuieH <- argv$bulkSizeH
bulkSuieL <- argv$bulkSizeL
highP <- argv$highP
lowP <- argv$lowP
highB <- argv$highB
lowB <- argv$lowB
popType <- argv$popType
minHPdp <- argv$minHPdp
maxHPdp <- argv$maxHPdp
minLPdp <- argv$minLPdp
maxLPdp <- argv$maxLPdp
minHBdp <- argv$minHBdp
maxHBdp <- argv$maxHBdp
minLBdp <- argv$minLBdp
maxLBdp <- argv$maxLBdp
winSize <- argv$winSize
minN <- argv$minN
width <- argv$width
height <- argv$height

if (FALSE) {
  ## file name
  filename <- "./xq_bsa_M1F_N1F.filter.SNPs.table"
  outPrefix <- "xqBSA"
  ## filter SNP Quality
  minQ <- 300
  ##
  bulkSizeH <- 30
  bulkSizeL <- 30
  ## sample name
  highP <- "s8417"
  lowP <- "FM195BML"
  highB <- "N1F"
  lowB <- "M1F"
  
  ##
  popType <- "F2"  # F2 ro RIL
  
  ## filter depth parameter
  minHPdp <- 10
  maxHPdp <- 100
  
  minLPdp <- 10
  maxLPdp <- 100
  
  minHBdp <- 10
  maxHBdp <- 100
  
  minLBdp <- 10
  maxLBdp <- 100
  
  ## sliding window parameter
  winSize <- 2000000
  #
  minN <- 10
  ##
  width <- 15
  height <- 5
}


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
dd1 <- df %>% filter(highParent.GT == paste(REF, REF, sep = "/") | highParent.GT == paste(ALT, ALT, sep = "/")) %>%
  filter(lowParent.GT == paste(REF, REF, sep = "/") | lowParent.GT == paste(ALT, ALT, sep = "/")) %>%
  filter(highParent.GT != lowParent.GT) %>% 
  filter(highParent.DP > minHPdp, highParent.DP < maxHPdp,
         lowParent.DP > minLPdp, lowParent.DP < maxLPdp,
         highBulk.DP > minHBdp, highBulk.DP < maxHBdp,
         lowBulk.DP > minLBdp, lowBulk.DP < maxLBdp)

dd2 <- dd1 %>% separate(highBulk.AD, c("highBulk.AD_0", "highBulk.AD_1"), sep = ",") %>%
  separate(lowBulk.AD, c("lowBulk.AD_0", "lowBulk.AD_1"), sep = ",") %>%
  separate(highBulk.PL, c("highBulk.PL_00", "highBulk.PL_01", "highBulk.PL_11"), sep = ",") %>%
  separate(lowBulk.PL, c("lowBulk.PL_00", "lowBulk.PL_01", "lowBulk.PL_11"), sep = ",")


dd3  <- dd2 %>% mutate(highBulk.AD = if_else(highParent.GT == paste(ALT, ALT, sep = "/"), 
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

write_tsv(dd3, paste(outPrefix, "reform.txt", sep = "."))

df <- importFromGATK(file = paste(outPrefix, "reform.txt", sep = "."),
                     highBulk = "highBulk",
                     lowBulk = "lowBulk",
                     #chr=c('C1','C2',C3)
                     )
chromColor <- read_tsv("./chromColor.txt")
df <- chromColor %>% left_join(df, by = "CHROM")

df <- runGprimeAnalysis(
  SNPset = df,
  windowSize = winSize,    # window size
  outlierFilter = "deltaSNP")

df <- runQTLseqAnalysis(
  SNPset = df,
  windowSize = winSize, #window size
  popStruc = popType,  #F2 or RIL
  bulkSize = c(bulkSuieH, bulkSuieL),  #bulk size, first is highbulk
  replications = 1000, 
  intervals = c(95, 99) # confidence interval
)

# 输出QTLseqr结算结果
outtb <- df %>% select(CHROM, POS, REF, ALT, AD_REF.LOW, AD_ALT.LOW, 
              SNPindex.LOW, AD_REF.HIGH, AD_ALT.HIGH, SNPindex.HIGH, 
              deltaSNP, nSNPs, tricubeDeltaSNP, CI_95, CI_99)
write_tsv(outtb, paste(outPrefix, "QTLseqrdeltaSNPindex.txt", sep = "."))
write_csv(outtb, paste(outPrefix, "QTLseqrdeltaSNPindex.csv", sep = "."))

#plotQTLStats(SNPset = df, var = "Gprime", plotThreshold = TRUE, q = 0.01)
#plotQTLStats(SNPset = df, var = "deltaSNP", plotIntervals = TRUE) + theme_half_open()
# add chrom color infomation

p <- df %>% filter(nSNPs > minN) %>%
  ggplot() +
  geom_line(aes(x = POS, y = CI_95), color = "gray") +
  geom_line(aes(x = POS, y = -CI_95), color = "gray") +
  geom_line(aes(x = POS, y = tricubeDeltaSNP, color = COLOR), size = 1) +
  ylim(-1, 1) +
  labs(x = NULL, y = "delta SNP index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_half_open() +
  theme(legend.position = "none") +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
ggsave(filename = paste(outPrefix, "deltaSNPindex.95CI.pdf", sep = "."), width = width, height = height)
ggsave(filename = paste(outPrefix, "deltaSNPindex.95CI.png", sep = "."), width = width, height = height, dpi = 500)

p <- df %>% filter(nSNPs > minN) %>%
  ggplot() +
  geom_line(aes(x = POS, y = CI_99), color = "gray") +
  geom_line(aes(x = POS, y = -CI_99), color = "gray") +
  geom_line(aes(x = POS, y = tricubeDeltaSNP, color = COLOR), size = 1) +
  ylim(-1, 1) +
  labs(x = NULL, y = "delta SNP index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_half_open() +
  theme(legend.position = "none") +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
ggsave(filename = paste(outPrefix, "deltaSNPindex.99CI.pdf", sep = "."), width = width, height = height)
ggsave(filename = paste(outPrefix, "deltaSNPindex.99CI.png", sep = "."), width = width, height = height, dpi = 500)

#export summary CSV
#getQTLTable(SNPset = df, alpha = 0.01, export = TRUE, fileName = "Gprime_QTL.csv",method="Gprime")
getQTLTable(SNPset = df, interval =95, export = TRUE, fileName = paste(outPrefix, "95CI.csv", sep = "."), method="QTLseq")
getQTLTable(SNPset = df, interval =99, export = TRUE, fileName = paste(outPrefix, "99CI.csv", sep = "."), method="QTLseq")
