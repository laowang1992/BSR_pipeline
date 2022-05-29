library(argparser, quietly=TRUE)

p <- arg_parser("Filter primer for 2 parent lines")

# Add command line arguments
p <- add_argument(p, "--primer", short = "-p", help = "A primer file containing several sample", type = "character")
p <- add_argument(p, "--sample1", short = "-v", help = "sample1", type = "character")
p <- add_argument(p, "--sample2", short = "-w", help = "sample2", type = "character")
p <- add_argument(p, "--chr", short = "-c", help = "chromosome", type = "character", default = NA)
p <- add_argument(p, "--start", short = "-s", help="start", type="numeric", default = -Inf)
p <- add_argument(p, "--end", short = "-e", help="end", type="numeric", default = Inf)
p <- add_argument(p, "--list", short = "-l", help = "list samples", flag = TRUE)
# Parse the command line arguments
argv <- parse_args(p)

primer <- argv$primer
s1 <- argv$sample1
s2 <- argv$sample2
chr <- argv$chr
s <- argv$start
e <- argv$end

test <- FALSE
if (test) {
  primer <- "./public.primer.csv.gz"
  s1 <- "s8800"
  s2 <- "s8804"
  chr <- "scaffoldA01"
  s <- 10000
  e <- 30000
}


if (argv$list) {
  df <- read.csv(primer, nrows = 1)
  col_name <- colnames(df)
  print(col_name)
  quit(save = "no")
}


library(tidyverse)
df <- read_csv(primer)

d1 <- df %>% select("CHROM", "POS", "ID", "REF", "ALT", "REF_len", "ALT_len", "DIFF", "QUAL", s1 = all_of(s1), s2 = all_of(s2), 
              "PrimerL1", "PrimerR1", "Hit1", "TmL1", "TmR1", "GCL1", "GCR1", "Length1",
              "PrimerL2", "PrimerR2", "Hit2", "TmL2", "TmR2", "GCL2", "GCR2", "Length2",
              "PrimerL3", "PrimerR3", "Hit3", "TmL3", "TmR3", "GCL3", "GCR3", "Length3") %>% 
  filter((s1 == "0/0" & s2 == "1/1") | (s1 == "1/1" & s2 == "0/0")) %>% filter(!is.na(Hit1))

colnames(d1) = c("CHROM", "POS", "ID", "REF", "ALT", "REF_len", "ALT_len", "DIFF", "QUAL", s1 , s2, 
                 "PrimerL1", "PrimerR1", "Hit1", "TmL1", "TmR1", "GCL1", "GCR1", "Length1",
                 "PrimerL2", "PrimerR2", "Hit2", "TmL2", "TmR2", "GCL2", "GCR2", "Length2",
                 "PrimerL3", "PrimerR3", "Hit3", "TmL3", "TmR3", "GCL3", "GCR3", "Length3")
if ((!is.na(chr)) && s == -Inf && e == Inf) {
  d2 <- d1 %>% filter(CHROM == chr) 
  write_csv(x = d2, file = paste(paste(s1, s2, sep = "_"), chr, "indel_primer", "csv", sep = "."))
}else if ((!is.na(chr)) && s != -Inf && e == Inf) {
  d2 <- d1 %>% filter(CHROM == chr, POS > s, POS < e)
  write_csv(x = d2, file = paste(paste(s1, s2, sep = "_"), paste(chr, s, "end", sep = "_"), "indel_primer", "csv", sep = "."))
}else if ((!is.na(chr)) && s == -Inf && e != Inf) {
  d2 <- d1 %>% filter(CHROM == chr, POS > s, POS < e)
  write_csv(x = d2, file = paste(paste(s1, s2, sep = "_"), paste(chr, "0", e, sep = "_"), "indel_primer", "csv", sep = "."))
}else if ((!is.na(chr)) && s != -Inf && e != Inf) {
  d2 <- d1 %>% filter(CHROM == chr, POS > s, POS < e)
  write_csv(x = d2, file = paste(paste(s1, s2, sep = "_"), paste(chr, s, e, sep = "_"), "indel_primer", "csv", sep = "."))
}else if (is.na(chr)) {
  d2 <- d1
  write_csv(x = d2, file = paste(paste(s1, s2, sep = "_"), "indel_primer", "csv", sep = "."))
}
