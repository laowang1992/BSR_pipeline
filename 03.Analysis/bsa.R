#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for BSA")

## Add command line arguments
#
p <- add_argument(p, "--input", help = "input file for BSA", type = "character")
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
  filename <- "./xq_bsa_M1F_N1F.filter.SNPs.txt"
  outPrefix <- "xq_bsa_M1F_N1F"
  ## filter SNP Quality
  minQ <- 300
  
  ## sample name
  # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s8417	FM195BML	M1F	N1F
  highP <- "s8417"
  lowP <- "FM195BML"
  highB <- "M1F"
  lowB <- "N1F"
  
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
  winSize <- 500000
  winStep <- 250000
  
  ## plot parameter
  minN <- 10
  width <-18
  height <- 5
}

library(tidyverse)
library(cowplot)
library(ggsci)
library(windowscanr)
library(RColorBrewer)

df <- read_tsv(file = filename)
df %>% select(CHROM, POS, REF, ALT, QUAL, all_of(c(highP, lowP, highB, lowB)))


chromColor <- read_tsv("./chromColor.txt")
dd <- df %>% filter(QUAL > minQ) #过滤QUAL小于等于50的行
l1<-nchar(dd$REF)	#求REF字符串长度
l2<-nchar(dd$ALT)	#求ALT字符串长度

#过滤indel行，保留SNP信息
lse <- l1 == 1 & l2 == 1
dd <- dd[lse, ]

dd1 <- dd %>% separate(highP, c(paste(highP, "geno", sep = "."), highP, "Y2"), sep=":") %>% 
  separate(lowP, c(paste(lowP, "geno", sep = "."), lowP, "Y4"), sep=":") %>%
  separate(highB, c(paste(highB, "geno", sep = "."), highB, "Y8"), sep=":") %>% 
  separate(lowB, c(paste(lowB, "geno", sep = "."), lowB, "Y10"), sep=":")
geno <- dd1 %>% select(paste(highP, "geno", sep = "."), 
                       paste(lowP, "geno", sep = "."), 
                       paste(highB, "geno", sep = "."), 
                       paste(lowB, "geno", sep = "."))	#sample基因型信息

#得到各sample每种基因型的覆盖度
dd2 <- dd1 %>% select(all_of(c(highP, lowP, highB, lowB)))
dd3 <- dd2 %>% 
  separate(highP,c(paste(highP, "ref.depth", sep = "."), paste(highP, "alt.depth", sep = ".")), sep=",", convert = T) %>%
  separate(lowP,c(paste(lowP, "ref.depth", sep = "."), paste(lowP, "alt.depth", sep = ".")), sep=",", convert = T) %>%
  separate(highB,c(paste(highB, "ref.depth", sep = "."), paste(highB, "alt.depth", sep = ".")), sep=",", convert = T)%>%
  separate(lowB,c(paste(lowB, "ref.depth", sep = "."), paste(lowB, "alt.depth", sep = ".")), sep=",", convert = T)

hh <- cbind(dd, geno, dd3) %>% as_tibble()
##########################筛选出两亲本均为纯合并且不同基因型的位点
hh <- hh %>% 
  filter((hh[, paste(highP, "geno", sep = ".")] == "0/0" |hh[ , paste(highP, "geno", sep = ".")] == "1/1") &
           (hh[, paste(lowP, "geno", sep = ".")] == "0/0" | hh[ , paste(lowP, "geno", sep = ".")] == "1/1") &
           (hh[, paste(highP, "geno", sep = ".")] != hh[, paste(lowP, "geno", sep = ".")]))
########去除混池样本中同一位点有第二种突变类型的行，即带有2的行，其实这一步有点多余
hh <- hh %>% 
  filter(!hh[, paste(highP, "geno", sep = ".")] =="0/2" | hh[, paste(highP, "geno", sep = ".")] =="1/2" | hh[, paste(highP, "geno", sep = ".")] =="2/2") %>%
  filter(!hh[, paste(lowP, "geno", sep = ".")] =="0/2" | hh[, paste(lowP, "geno", sep = ".")] =="1/2" | hh[, paste(lowP, "geno", sep = ".")] =="2/2")
############################计算每个sample每个位点总的覆盖度，个别位点并不等于DP值
hh[, paste(highP, "sum", sep = ".")] <- hh[, paste(highP, "ref.depth", sep = ".")] + hh[, paste(highP, "alt.depth", sep = ".")]
hh[, paste(lowP, "sum", sep = ".")] <- hh[, paste(lowP, "ref.depth", sep = ".")] + hh[, paste(lowP, "alt.depth", sep = ".")]
hh[, paste(highB, "sum", sep = ".")] <- hh[, paste(highB, "ref.depth", sep = ".")] + hh[, paste(highB, "alt.depth", sep = ".")]
hh[, paste(lowB, "sum", sep = ".")] <- hh[, paste(lowB, "ref.depth", sep = ".")] + hh[, paste(lowB, "alt.depth", sep = ".")]

hh<-na.omit(hh)	#删除含有NA的行


######过滤混池样本中总覆盖深度小于等于4的位点
pdf(paste(outPrefix, "depth_desity.pdf", sep = "."), height = 6, width = 6)
par(mfrow = c(2, 2))
plot(density(hh[, paste(highP, "sum", sep = ".")][[1]], width = 2), main = highP)
plot(density(hh[, paste(lowP, "sum", sep = ".")][[1]], width = 2), main = lowP)
plot(density(hh[, paste(highB, "sum", sep = ".")][[1]], width = 2), main = highB)
plot(density(hh[, paste(lowB, "sum", sep = ".")][[1]], width = 2), main = lowB)
dev.off()
png(paste(outPrefix, "depth_desity.png", sep = "."), units = "in", height = 6, width = 6, res = 300)
par(mfrow = c(2, 2))
plot(density(hh[, paste(highP, "sum", sep = ".")][[1]], width = 2), main = highP)
plot(density(hh[, paste(lowP, "sum", sep = ".")][[1]], width = 2), main = lowP)
plot(density(hh[, paste(highB, "sum", sep = ".")][[1]], width = 2), main = highB)
plot(density(hh[, paste(lowB, "sum", sep = ".")][[1]], width = 2), main = lowB)
dev.off()

hh <- hh %>% filter(hh[, paste(highP, "sum", sep = ".")] > minHPdp & hh[, paste(highP, "sum", sep = ".")] < maxHPdp)
hh <- hh %>% filter(hh[, paste(lowP, "sum", sep = ".")] > minLPdp & hh[, paste(lowP, "sum", sep = ".")] < maxLPdp)
hh <- hh %>% filter(hh[, paste(highB, "sum", sep = ".")] > minHBdp & hh[, paste(highB, "sum", sep = ".")] < maxHBdp)
hh <- hh %>% filter(hh[, paste(lowB, "sum", sep = ".")] > minLBdp & hh[, paste(lowB, "sum", sep = ".")] < maxLBdp)


##
SNPnumber <- hh %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, path = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
write_csv(x = SNPnumber, path = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))
##
options(scipen = 200)
colourCount = dim(chromColor)[[1]]
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
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


#########################
### get reads depth of P1 and P2 in pool
hh[, paste(highB, highP, "depth", sep = ".")] <- hh[, paste(highB, "ref.depth", sep = ".")]
hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(highB, highP, "depth", sep = ".")] <- hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(highB, "alt.depth", sep = ".")]

hh[, paste(highB, lowP, "depth", sep = ".")] <- hh[, paste(highB, "alt.depth", sep = ".")]
hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(highB, lowP, "depth", sep = ".")] <- hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(highB, "ref.depth", sep = ".")]

hh[, paste(lowB, highP, "depth", sep = ".")] <- hh[, paste(lowB, "ref.depth", sep = ".")]
hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(lowB, highP, "depth", sep = ".")] <- hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(lowB, "alt.depth", sep = ".")]

hh[, paste(lowB, lowP, "depth", sep = ".")] <- hh[, paste(lowB, "alt.depth", sep = ".")]
hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(lowB, lowP, "depth", sep = ".")] <- hh[hh[, paste(highP, "geno",sep = ".")] == "1/1", paste(lowB, "ref.depth", sep = ".")]

table(paste(hh[, paste(highP, "geno", sep = ".")][[1]], hh[, paste(lowP, "geno", sep = ".")][[1]], sep="_"))
#############################
hh$chr<-hh$CHROM

w <- hh %>% select(CHROM, POS, all_of(c(paste(highB, highP, "depth", sep = "."), paste(highB, lowP, "depth", sep = "."))))
b <- hh %>% select(CHROM, POS, all_of(c(paste(lowB, highP, "depth", sep = "."), paste(lowB, lowP, "depth", sep = "."))))
#############################
#ww$R.p1.depth<-ww$R.p1.depth+0.5
#bb$S.p1.depth<-bb$S.p1.depth+0.5
w1 <- winScan(x = w,
              groups = "CHROM",
              position = "POS",
              values = c(paste(highB, highP, "depth", sep = "."), paste(highB, lowP, "depth", sep = ".")),
              win_size = winSize,
              win_step = winStep,
              funs = c("sum"))
b1 <- winScan(x = b,
              groups = "CHROM",
              position = "POS",
              values = c(paste(lowB, highP, "depth", sep = "."), paste(lowB, lowP, "depth", sep = ".")),
              win_size = winSize,
              win_step = winStep,
              funs = c("sum"))
#w1$w<-w1$R.p2.depth_sum/(w1$R.p1.depth_sum)
#b1$b<-b1$S.p2.depth_sum/(b1$S.p1.depth_sum)
#w1$w<-w1$R.p2.depth_sum/(w1$R.p1.depth_sum+w1$R.p2.depth_sum)
#b1$b<-b1$S.p2.depth_sum/(b1$S.p1.depth_sum+b1$S.p2.depth_sum)
###
w1$w <- w1[, paste(highB, highP, "depth_sum", sep = ".")] / (w1[, paste(highB, highP, "depth_sum", sep = ".")] + w1[, paste(highB, lowP, "depth_sum", sep = ".")])
b1$b <- b1[, paste(lowB, highP, "depth_sum", sep = ".")] / (b1[, paste(lowB, highP, "depth_sum", sep = ".")] + b1[, paste(lowB, lowP, "depth_sum", sep = ".")])
#plot(w1$w - b1$b,col=factor(w1$CHROM),pch=20,cex=0.5,xaxt="n",ylab="SNP reads ratio")

d <- tibble(CHROM = w1$CHROM, w1$win_start, w1$win_mid, w1$win_end, 
            SNP_N = w1[, paste(highB, highP, "depth_n", sep = ".")], 
            INDEX_W = w1$w, INDEX_B = b1$b, delta_INDEX = w1$w - b1$b)
colnames(d) <- c("CHROM", "win_start", "win_mid", "win_end", "SNP_N", 
                 paste(highB, "SNP_INDEX", sep = "."), 
                 paste(lowB, "SNP_INDEX", sep = "."), "delta_SNP_INDEX")

write_tsv(x = d, path = paste(outPrefix, "SlidingWindow_SNPindex.txt", sep = "."))
write_csv(x = d, path = paste(outPrefix, "SlidingWindow_SNPindex.csv", sep = "."))

# 加上绘图颜色信息
d <- chromColor %>% left_join(d, by = "CHROM")

#mytheme<-theme(axis.title = element_text(face = "bold",size = 10),
#               axis.text = element_text(face = "bold",size = 9),
#               panel.background = element_rect(fill = "white",color = "darkblue"),
#               panel.grid.major.y = element_line(color = "grey",linetype = 1),
#               panel.grid.minor.y = element_line(color = "grey",linetype = 2),
#               panel.grid.minor.x = element_blank(),
#               legend.position = "none")

Pdelta <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = delta_SNP_INDEX)) +
  geom_point(aes(color = COLOR), size = 0.7) +
  ylim(-1, 1) +
  labs(x="", y="delta SNP index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_minimal_grid() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
ggsave(Pdelta, filename = paste(outPrefix, "delta_SNP_index.pdf", sep = "."), width = width, height = height)
ggsave(Pdelta, filename = paste(outPrefix, "delta_SNP_index.png", sep = "."), width = width, height = height, dpi = 500)
#plot(b1$b/w1$w,pch=20,xaxt="n",ylab="SNP reads ratio")

Phigh <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = d[d$SNP_N > minN, paste(highB, "SNP_INDEX", sep = ".")][[1]])) +
  geom_point(aes(color = COLOR), size = 0.7) +
  ylim(0, 1) +
  labs(x="", y="SNP index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_minimal_grid() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
ggsave(Phigh, filename = paste(outPrefix, highB, "SNP_index.pdf", sep = "."), width = width, height = height)
ggsave(Phigh, filename = paste(outPrefix, highB, "SNP_index.png", sep = "."), width = width, height = height, dpi = 500)

Plow <- ggplot(filter(d, SNP_N > minN), aes(x = win_mid, y = d[d$SNP_N > minN, paste(lowB, "SNP_INDEX", sep = ".")][[1]])) +
  geom_point(aes(color = COLOR), size = 0.7) +
  ylim(0, 1) +
  labs(x="", y="SNP index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_minimal_grid() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
ggsave(Plow, filename = paste(outPrefix, lowB, "SNP_index.pdf", sep = "."), width = width, height = height)
ggsave(Plow, filename = paste(outPrefix, lowB, "SNP_index.png", sep = "."), width = width, height = height, dpi = 500)

