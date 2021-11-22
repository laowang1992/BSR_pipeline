library(jsonlite)
library(tidyverse)

all <- tibble()

#sample <- "DiParent"
for (sample in list.files(pattern = "json") %>% basename() %>% str_remove(".json")) {
  df <- jsonlite::fromJSON(paste(sample, "json", sep = "."))
  raw <- as_tibble(df$summary$before_filtering) %>% mutate(sample = sample)
  clean <- as_tibble(df$summary$after_filtering) %>% mutate(sample = sample)
  
  stat <- tibble(sample = sample) %>% 
    left_join(raw, by = "sample") %>%
    left_join(clean, by = "sample", suffix = c(".raw", ".clean"))
  
  all <- rbind(all, stat)
}

all <- all %>% mutate(effective.rate = total_bases.clean / total_bases.raw)
write_csv(x = all, path = "data_stat.csv")
write_tsv(x = all, path = "data_stat.txt")
