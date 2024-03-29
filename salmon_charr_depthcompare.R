library(tidyverse)
library(data.table)
library(matrixStats)

setwd("~/Desktop/Working/Salmon_eDNA_nuc/")

#import charr depth data

Charr.5kdepths <- fread("Charr.5kdepths") %>%  mutate(Species = "Charr")
Charr_rmed <- rowMedians(as.matrix(Charr.5kdepths[,4:54]))
Charr_rsum <- rowSums(as.matrix(Charr.5kdepths[,4:54]))

Charr.5kdepths <- bind_cols(Charr.5kdepths, Charr_rsum)  %>%  
  mutate(sumdepth  = ...56) %>%  
  select(-...56)
colnames(Charr.5kdepths)[1:3] <- c("Chrom", "Winstart", "Winend")


Charr.5kdepths <- Charr.5kdepths %>%  group_by(Chrom) %>%  summarise(median_depth = median(sumdepth)) %>%  inner_join(Charr.5kdepths)
Charr.5kdepths <-  Charr.5kdepths %>%  mutate(chromdepth_sum_std = sumdepth/median_depth)


#import salmon depth data

Salmon.5kdepths <- fread("Salmo.5kdepths") %>%  mutate(Species = "Salmon")
Salmon_rmed <- rowMedians(as.matrix(Salmon.5kdepths[,4:63]))
Salmon_rsum <- rowSums(as.matrix(Salmon.5kdepths[,4:63]))

Salmon.5kdepths <- bind_cols(Salmon.5kdepths, Salmon_rsum)  %>%  
  mutate(sumdepth  = ...65) %>%  
  select(-...65)
colnames(Salmon.5kdepths)[1:3] <- c("Chrom", "Winstart", "Winend")


Salmon.5kdepths <- Salmon.5kdepths %>%  group_by(Chrom) %>%  summarise(median_depth = median(sumdepth)) %>%  inner_join(Salmon.5kdepths)
Salmon.5kdepths <-  Salmon.5kdepths %>%  mutate(chromdepth_sum_std = sumdepth/median_depth)

#get regions of 0 coverage in charr, write out bed format file

Charr_0depth <- Charr.5kdepths %>% filter(sumdepth == 0) %>%   select(Chrom, Winstart, Winend, sumdepth, median_depth)
Charr_salm <- Salmon.5kdepths %>%  select(Chrom, Winstart, Winend, sumdepth, median_depth) %>% 
  inner_join(Charr_0depth, by = c("Chrom", "Winstart", "Winend"))

Charr_salm  <- Charr_salm %>%  mutate(deltadepth = abs(sumdepth.x - sumdepth.y))
Charr_salm_top10_deltacov.bed <- Charr_salm  %>%  slice_max(deltadepth, n = 10) %>%  select(Chrom, Winstart, Winend)

write.table(Charr_salm_top10_deltacov.bed, "Charr_salm_top10_deltacov.bed", col.names = F, row.names = F, sep = "\t", quote = F)




