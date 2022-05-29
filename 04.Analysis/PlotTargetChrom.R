library(tidyverse)
library(cowplot)
options(scipen=200)

inputPrefix <- "BSA"
CI <- "95"
qtlfile <- paste(inputPrefix, paste(CI, "CI", sep = ""), "csv", sep = ".")
indexfile <- paste(inputPrefix, "QTLseqrdeltaSNPindex", "csv", sep = ".")
chr <- read_tsv("./chromColor.txt")
len <- read_tsv("./ref.len", col_names = c("CHROM", "Len"))
minN <- 5

#
qtl <- read_csv(qtlfile) %>% 
  dplyr::filter(length > 0, nSNPs > minN)
index <- read_csv(indexfile)

for (i in seq_along(qtl$CHROM)) {
  chr <- qtl[i, "CHROM"][[1]]
  start <- qtl[i, "start"][[1]]
  end <- qtl[i, "end"][[1]]
  df <- index %>% filter(CHROM == chr)
  p <- ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300", alpha = .005) +
    geom_line(aes(x = POS, y = tricubeDeltaSNP), color = "blue") +
    geom_line(aes(x = POS, y = CI_99), color = "gray60") + 
    geom_line(aes(x = POS, y = -CI_99), color = "gray60") +
    geom_hline(aes(yintercept=0)) + 
    labs(x = chr, y = "delta SNP index") +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_half_open()
  ggsave(p, filename = paste(paste(chr, start, end, sep = "_"), paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
  print(paste(chr, "has been done...", sep = " "))
}

options(scipen=200)
plotTargetChrom <- function(df, CI = 95, minN = 0, outPrefix, chr, len){
  interval <- read_csv(file = paste(outPrefix, paste(CI, "CI", sep = ""), "csv", sep = "."))
  if (nrow(interval) > 0) {
    for (i in unique(interval$CHROM)) {
      d <- df %>% filter(CHROM == i, nSNPs > minN) %>% left_join(chr, by = "CHROM") %>% 
        dplyr::select(CHROM, LABEL, POS, CI = all_of(paste("CI", CI, sep = "_")), tricubeDeltaSNP)
      inter <- interval %>% filter(CHROM == i) %>% left_join(chr, by = "CHROM")
      chrLen <- len$Len[len$CHROM == i]
      p <- ggplot() +
        geom_rect(data = inter, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300") +
        geom_line(data = d, aes(x = POS, y = CI), color = "gray60") + 
        geom_line(data = d, aes(x = POS, y = -CI), color = "gray60") +
        geom_line(data = d, aes(x = POS, y = tricubeDeltaSNP), color = "blue") +
        coord_cartesian(xlim = c(0, chrLen), ylim = c(-1, 1)) +
        geom_hline(aes(yintercept=0)) + 
        labs(x = chr$LABEL[chr$CHROM == i], y = "Delta SNP index") +
        theme_half_open()
      ggsave(p, filename = paste(outPrefix, chr$LABEL[chr$CHROM == i], paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      print(paste(chr$LABEL[chr$CHROM == i], "has been done...", sep = " "))
    }
  }
}

plotTargetChrom(df = index, CI = CI, minN = minN, outPrefix = inputPrefix, chr = chr, len = len)
