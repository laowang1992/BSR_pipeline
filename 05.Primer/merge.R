library(tidyverse)
#library(xlsx)

argv <- commandArgs(trailingOnly = T)

#argv[1] <- "public.indels.filter.txt"
#argv[2] <- "public.primer.txt"
gt <- read_tsv(file = argv[1]) %>%
  mutate(REF_len = str_length(REF), ALT_len = str_length(ALT), 
         DIFF = abs(REF_len - ALT_len))
primer <- read_tsv(file = argv[2])

df <- gt %>% right_join(primer, by = "ID")

write_csv(x = df, file = str_replace(argv[2], "primer.txt", "primer.csv"))
#write.xlsx(x = df, file = str_replace(argv[2], "primer.txt", "primer.xlsx"))
