library(tidyverse)
#library(xlsx)

argv <- commandArgs(trailingOnly = T)

test <- FALSE
if (test) {
  argv[1] <- "public.indels.filter.txt"
  argv[2] <- "public.primer.txt"
}

gt <- read_tsv(file = argv[1]) %>%
  mutate(REF_len = str_length(REF), ALT_len = str_length(ALT), 
         DIFF = abs(REF_len - ALT_len))
primer <- read_tsv(file = argv[2])

p0 <- primer %>% filter(NUM == 0) %>% select(ID, PrimerL1 = PrimerL, PrimerR1 = PrimerR, Hit1 = Hit, TmL1 = TmL, TmR1 = TmR, GCL1 = GCL, GCR1 = GCR, Length1 = Length)
p1 <- primer %>% filter(NUM == 1) %>% select(ID, PrimerL2 = PrimerL, PrimerR2 = PrimerR, Hit2 = Hit, TmL2 = TmL, TmR2 = TmR, GCL2 = GCL, GCR2 = GCR, Length2 = Length)
p2 <- primer %>% filter(NUM == 2) %>% select(ID, PrimerL3 = PrimerL, PrimerR3 = PrimerR, Hit3 = Hit, TmL3 = TmL, TmR3 = TmR, GCL3 = GCL, GCR3 = GCR, Length3 = Length)

df <- gt %>% left_join(p0, by = "ID") %>% 
  left_join(p1, by = "ID") %>%
  left_join(p2, by = "ID")

#df <- gt %>% right_join(primer, by = "ID")

write_csv(x = df, file = str_replace(argv[2], "primer.txt", "primer.csv"))
#write.xlsx(x = df, file = str_replace(argv[2], "primer.txt", "primer.xlsx"))
