library(tidyverse)
library(cowplot)

# 定义参数
input <- "./RY.QTLseqrdeltaSNPindex.txt"
output <- "RY_chrA07.CI99"
chr <- "scaffoldA07"
start <- 21963525
end <- 30953727
minN <- 10

# 读取数据
df <- read_tsv(input) %>%
  filter(CHROM == chr, nSNP >= minN)

# 绘图并输出
p <- df %>% filter(CHROM == "scaffoldA07", nSNPs >= 10) %>%
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),fill = '#FF3300', alpha = .005) +
  geom_line(aes(x = POS, y = tricubeDeltaSNP), color = "blue") +
  geom_line(aes(x = POS, y = CI_99), color = "gray60") + 
  geom_line(aes(x = POS, y = -CI_99), color = "gray60") +
  geom_hline(aes(yintercept=0)) + 
  scale_x_continuous(breaks = c(0, 10000000, 20000000, 30000000), 
                     labels = c("0M", "10M", "20M", "30M")) +
  labs(x = "chrA07", y = "delta SNP index") +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_half_open()
#ggsave(pRY, filename = "RY_chrA07.CI99.pdf", height = 3.5, width = 4.5)
ggsave(pRY, filename = paste(output, "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)

