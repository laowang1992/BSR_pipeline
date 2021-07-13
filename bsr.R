library(tidyverse)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(windowscanr)

input <- "./BSR.filter.SNPs.txt"
(df <- read_tsv(input))
(chromColor <- read_tsv("./chromColor.txt"))
dd <- df %>% filter(QUAL > 300)

l1<-nchar(dd$REF)	#求REF字符串长度
l2<-nchar(dd$ALT)	#求ALT字符串长度

#过滤indel行，保留SNP信息
lse <- l1 == 1 & l2 == 1
dd <- dd[lse, ]

dd1 <- dd %>% separate(LL, c("LL.geno", "LL", "Y2"), sep = ":") %>%
  separate(RL, c("RL.geno", "RL", "Y2"), sep = ":")
geno <- dd1 %>% select(LL.geno, RL.geno)

dd2 <- dd1 %>% select(LL, RL)
dd3 <- dd2 %>% 
  separate(LL, c("LL.ref.depth", "LL.alt.depth"), sep = ",", convert = T) %>%
  separate(RL, c("RL.ref.depth", "RL.alt.depth"), sep = ",", convert = T) 

hh <- cbind(dd, geno, dd3) %>% as_tibble()
hh %>% group_by(LL.geno, RL.geno) %>% count()

hh <- hh %>% mutate(LL.sum = LL.ref.depth + LL.alt.depth,
                    RL.sum = RL.ref.depth + RL.alt.depth)

pdf("depth_desity.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
plot(density(hh$LL.sum, width = 1), main = "LL", xlim = c(0, 100))
plot(density(hh$RL.sum, width = 1), main = "RL", xlim = c(0, 100))
dev.off()
png("depth_desity.png", height = 3, width = 6, units = "in", res = 500)
par(mfrow = c(1, 2))
plot(density(hh$LL.sum, width = 1), main = "LL", xlim = c(0, 100))
plot(density(hh$RL.sum, width = 1), main = "RL", xlim = c(0, 100))
dev.off()
hh <- hh %>% filter(LL.sum > 10, RL.sum > 10)

SNPnumber <- hh %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, path = "SNP_number_per_chr.txt")
write_csv(x = SNPnumber, path = "SNP_number_per_chr.csv")


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

ggsave(Phist, filename = "SNP_distribution_histogram.pdf", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5)
ggsave(Phist, filename = "SNP_distribution_histogram.png", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5, dpi = 500)

hh <- hh %>% mutate(LL.ref.rate = LL.ref.depth / LL.sum,
                    LL.alt.rate = LL.alt.depth / LL.sum,
                    RL.ref.rate = RL.ref.depth / RL.sum,
                    RL.alt.rate = RL.alt.depth / RL.sum,
                    ED = sqrt((LL.ref.rate - RL.ref.rate)^2 + (LL.alt.rate - RL.alt.rate)^2),
                    ED4 = ED^4)

w <- winScan(x = hh,
             groups = "CHROM",
             position = "POS",
             values = c("ED", "ED4"),
             win_size = 500000,
             win_step = 250000,
             funs = c("mean"))

w <- w %>% select(CHROM, win_start, win_end, win_mid, ED = ED_mean, ED4 = ED4_mean, SNP_number = ED_n)
write_tsv(file = "BSR_ED.txt", x = w)
write_csv(file = "BSR_ED.csv", x = w)

d <- chromColor %>% left_join(w, by = "CHROM")
d
P1 <- ggplot(filter(d, SNP_number > 10), aes(x = win_mid, y = ED)) +
  geom_point(aes(color = COLOR), size = 0.7) +
  #ylim(-1, 1) +
  labs(x="", y="ED") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
P1
ggsave(P1, filename = "BSR_ED.pdf", width = 15, height = 4)
ggsave(P1, filename = "BSR_ED.png", width = 15, height = 4, dpi = 500)

P2 <- ggplot(filter(d, SNP_number > 10), aes(x = win_mid, y = ED4)) +
  geom_point(aes(color = COLOR), size = 0.7) +
  #ylim(-1, 1) +
  labs(x="", y="ED4") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
P2
ggsave(P2, filename = "BSR_ED4.pdf", width = 15, height = 4)
ggsave(P2, filename = "BSR_ED4.png", width = 15, height = 4, dpi = 500)

###
# chrA10
P_A10_1 <- ggplot(filter(d, SNP_number > 10, CHROM == "scaffoldA10"), aes(x = win_mid, y = ED)) +
  geom_point(color = "orange", size = 1) +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0M", "5M", "10M", "15M", "20M", "25M")) +
  labs(x = "Position", y = "ED") +
  theme_half_open()
P_A10_1
ggsave(P_A10_1, filename = "A10_ED.pdf", height = 3, width = 4)
ggsave(P_A10_1, filename = "A10_ED.png", height = 3, width = 4, dpi = 500)

P_A10_2 <- ggplot(filter(d, SNP_number > 10, CHROM == "scaffoldA10"), aes(x = win_mid, y = ED4)) +
  geom_point(color = "orange", size = 1) +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0M", "5M", "10M", "15M", "20M", "25M")) +
  labs(x = "Position", y = "ED4") +
  theme_half_open()
P_A10_2
ggsave(P_A10_2, filename = "A10_ED4.pdf", height = 3, width = 4)
ggsave(P_A10_2, filename = "A10_ED4.png", height = 3, width = 4, dpi = 500)



x <- hh %>% filter((LL.geno == "0/0" & RL.geno == "1/1") | (LL.geno == "1/1" & RL.geno == "0/0"))
x <- chromColor %>% left_join(x, by = "CHROM")
ggplot(x, aes(x = POS)) + 
  geom_histogram() + 
  facet_grid(LABEL ~ ., scales = "free_x", space = "free_x")

