library(DESeq2)
library(ggplot2)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())

load("../Data/TTdds_inter.Rda")
load("../Data/TTdds_pool.Rda")
load("../Data/TTdds_condition.Rda")


gene_map_table <- read_lines(file = "../../../LearnMemTolQTL/ProcessedData/gene_map_table_fb_2018_04.tsv",
                             skip = 6) %>% 
  str_split(pattern = "\t", simplify = TRUE) %>% 
  as_tibble() %>% 
  filter(V1 == "Dmel")

colnames(gene_map_table) <- c('spp','gname','FBgn','v3','cyt','pos')


gene_map_table$chr <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,1]
temp.s1 <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,2] %>% 
  str_split(fixed(".."),simplify=TRUE)
gene_map_table$startp <- as.numeric(temp.s1[,1])

gene_map_table$stopp <- as.numeric(str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1])


plotCounts(TTlrt_pool_deseq, gene="FBgn0011672", intgroup=c("pool","condition"))

plotCounts(TTlrt_pool_deseq, gene="FBgn0038851", intgroup=c("pool","condition"))

plotCounts(TTlrt_pool_deseq, gene="FBgn0025865", intgroup=c("pool","condition"))

