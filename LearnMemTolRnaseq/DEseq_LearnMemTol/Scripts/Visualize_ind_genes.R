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

load(file ="../../../LearnMemTolQTL/ProcessedData/Peaks_wCIs.rda")

lowc <- ci.peak$Tolsqrtvariable[which.max(ci.peak$Tolsqrtvariable$LL), 'lpR6']
highc <- ci.peak$Tolsqrtvariable[which.max(ci.peak$Tolsqrtvariable$LL), 'upR6']


gglist <- subset(gene_map_table, chr=="3R" & startp <= highc & stopp >=lowc)



#loop plots
#output all


pdf(file="../Plots/Ind_genes_TT_peak.pdf", width=16, height = 24)
par(mfrow=c(6,4))

for(ii in 1:nrow(gglist))
{

  
plotCounts(TTlrt_pool_deseq, gene=gglist$FBgn[ii], intgroup=c("pool","condition"), main=gglist$gname[ii])

}

dev.off()


