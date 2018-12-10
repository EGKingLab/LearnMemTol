library(DSPRqtl)
library(tidyverse)


#load in all the datasets

#significant qtl peaks
load(file ="../ProcessedData/Peaks_wCIs.rda")
str(ci.peak)

#DE genes for Learning
load(file = "../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/LearnresSVOrder.Rda")
str(LearnresSVorder)
LearnresSVorder$FBgn <- rownames(LearnresSVorder)
colnames(LearnresSVorder)

Learn_sig_genes <- as.data.frame(LearnresSVorder)

#colnames(LearnresSVorder)[which(names(LearnresSVorder) == "rownames")] <- "FBgn"

#DE genes for Memory
load(file = "../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/MemresSVOrder.Rda")
str(MemresSVorder)
MemresSVorder$FBgn <- rownames(MemresSVorder)
colnames(MemresSVorder)

Mem_sig_genes <- as.data.frame(MemresSVorder)

#DE genes thermal Tolerance 

#interaction
load(file="../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/TTlrt_inter.Rda")
str(TTlrt_inter)
TTlrt_inter$FBgn <- rownames(TTlrt_inter)
colnames(TTlrt_inter)

ThermTol_inter_sig_genes <- as.data.frame(TTlrt_inter)

#condition
load(file="../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/TTlrt_condition.Rda")
str(TTlrt_condition)
TTlrt_condition$FBgn <- rownames(TTlrt_condition)
colnames(TTlrt_condition)

ThermTol_condition_sig_genes <- as.data.frame(TTlrt_condition)


#pool
load(file="../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/TTlrt_pool.Rda")
str(TTlrt_pool)
TTlrt_pool$FBgn <- rownames(TTlrt_pool)
colnames(TTlrt_pool)

ThermTol_pool_sig_genes <- as.data.frame(TTlrt_pool)


#gene list from Fly Base

gene_map_table <- read_lines(file = "../ProcessedData/gene_map_table_fb_2018_04.tsv",
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

str(gene_map_table)
colnames(gene_map_table)


#gene_map_table <- as.numeric(gene_map_table$startp)

#merge DE genes with gene list from Fly base (merge by FBgn num)
Learn_gene_list_merged <- merge(Learn_sig_genes, gene_map_table, by="FBgn")
Mem_gene_list_merged <- merge(Mem_sig_genes, gene_map_table, by="FBgn")
ThermTol_inter_gene_list_merged <- merge(ThermTol_inter_sig_genes, gene_map_table, by="FBgn")
ThermTol_pool_gene_list_merged <- merge(ThermTol_pool_sig_genes, gene_map_table, by="FBgn")
ThermTol_condition_gene_list_merged <- merge(ThermTol_condition_sig_genes, gene_map_table, by="FBgn")



str(Learn_gene_list_merged)
colnames(Learn_gene_list_merged)


#select specific columns

Learn_gene_list_sub <- Learn_gene_list_merged[c(1,3,6,7,9,13:15)]
colnames(Learn_gene_list_sub)
Learn_gene_list_sub <- subset(Learn_gene_list_sub, padj <=0.05)

Mem_gene_list_sub <- Mem_gene_list_merged[c(1,3,6,7,9,13:15)]
colnames(Mem_gene_list_sub)
Mem_gene_list_sub <- subset(Mem_gene_list_sub, padj <=0.05)

ThermTol_inter_gene_list_sub <- ThermTol_inter_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_inter_gene_list_sub)
ThermTol_inter_gene_list_sub <- subset(ThermTol_inter_gene_list_sub, padj <=0.05)

ThermTol_pool_gene_list_sub <- ThermTol_pool_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_pool_gene_list_sub)
ThermTol_pool_gene_list_sub <- subset(ThermTol_pool_gene_list_sub, padj <=0.05)

ThermTol_condition_gene_list_sub <- ThermTol_condition_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_condition_gene_list_sub)
ThermTol_conditiion_gene_list_sub <- subset(ThermTol_condition_gene_list_sub, padj <=0.05)


#filter dataset to narrow down genes under the qtl peaks
#USE R6

#Learning
foc.peak <- ci.peak[[1]][1,]

learn.list <-vector(mode='list', length=nrow(ci.peak$LearnPowerTrans))

Learn_totdf<-data.frame(NULL)
for(row in 1:nrow(ci.peak$LearnPowerTrans)) 
  
  {
  
  foc.peak<-ci.peak[['LearnPowerTrans']][row,]
  lw <- (which(Learn_gene_list_sub$chr==foc.peak$chr & 
                 ((Learn_gene_list_sub$startp <= foc.peak$upR6 & Learn_gene_list_sub$stopp >= foc.peak$lpR6))))  
  
  Learn_genes_under_peak<-Learn_gene_list_sub[lw,]
  learn.list[[row]] <-Learn_genes_under_peak
  
 }

save(learn.list,file="../ProcessedData/Learn_genes_under_peak.rda")


#Memory 

Mem_foc_peak <- ci.peak[["Memory_Mean"]][1,]

Mem.list <-vector(mode='list', length=nrow(ci.peak$Memory_Mean))

Mem_totdf<-data.frame(NULL)
for(row in 1:nrow(ci.peak$Memory_Mean)) 
  
  {
  
  foc.peak<-ci.peak[['Memory_Mean']][row,]
  mw <- (which(Mem_gene_list_sub$chr==foc.peak$chr & 
                 ((Mem_gene_list_sub$startp <= foc.peak$upR6 & Mem_gene_list_sub$stopp >= foc.peak$lpR6))))  
  
 
  Mem_genes_under_peak <- Mem_gene_list_sub[mw,]
  
  Mem.list[[row]] <-Mem_genes_under_peak
  
}


Mem_genes_under_peak


save(Mem.list,file="../ProcessedData/Mem_genes_under_peak.rda")


