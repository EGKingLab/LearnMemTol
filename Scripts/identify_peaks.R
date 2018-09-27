library(DSPRqtl)
library(tidyverse)
data("positionlist_wgenetic")
load(file = "/home/kingeg/Projects/DSPRgeneral/Convert5_6/DSPR_5_6_parsed.rda")

source("mappingfunctions.R")

load(file="../Data/sig_ths.rda")#fdr.out
load(file="../Data/Lodscores_3traits.rda")#obs.LL

#th.l <- fdr.out$fwer
th.l <- fdr.out$fdr$th[min(which(fdr.out$fdr$fdr<=0.05))]

#find peaks
peak.i <- apply(obs.LL[,1:3],2,function(x) peakInfo(x,cM=poslist[,c('chr','Gpos')] ,th=th.l, tol.dist=3))



ci.peak <- vector(mode="list", length=3)
names(ci.peak)<-names(peak.i)

for(kk in 1:3){
qq <- poslist[,1:3]
qq$LOD <- obs.LL[,kk]
pp.set <- peak.i[[kk]]
pp.set$up <- NA
pp.set$lp <- NA
pp.set$ug <- NA
pp.set$lg <- NA
pp.set$ulod <- NA
pp.set$llod <- NA


for(jj in 1:nrow(pp.set))
{
ff<-findCI(pp.set$chr[jj],pp.set$Ppos[jj], qtldat=qq, method='BCI')
pp.set$up[jj] <- ff$Ppos[2]
pp.set$lp[jj] <- ff$Ppos[1]
pp.set$ug[jj] <- ff$Gpos[2]
pp.set$lg[jj] <- ff$Gpos[1]
pp.set$ulod[jj] <- ff$LOD[2]
pp.set$llod[jj] <- ff$LOD[1] 
} 

pp.set[,c('chr','Ppos')] <- merge(pp.set, coord.table, by.x=c('chr','Ppos'), by.y=c('R5chr','R5pos'),sort=FALSE)[,c('R6chr','R6pos')]
pp.set[,c('up')] <- merge(pp.set, coord.table, by.x=c('chr','up'), by.y=c('R5chr','R5pos'),sort=FALSE)[,c('R6pos')]
pp.set[,c('lp')] <- merge(pp.set[,c('chr','lp')], coord.table, by.x=c('chr','lp'), by.y=c('R5chr','R5pos'),sort=FALSE)[,c('R6pos')]
pp.set$Ppos<- as.numeric(pp.set$Ppos)
pp.set$up<- as.numeric(pp.set$up)
pp.set$lp<- as.numeric(pp.set$lp)

pp.set<-pp.set[order(pp.set$chr,pp.set$Ppos),]

ci.peak[[kk]]<-pp.set
}



save(ci.peak,file="../Data/Peaks_wCIs.rda")


###########


#load in all the datasets

#significant qtl peaks
load(file ="../Data/Peaks_wCIs.rda")
str(ci.peak)

#DE genes for Learning
load(file = "../Data/LearnresSVOrder.Rda")
str(LearnresSVorder)
LearnresSVorder$FBgn <- rownames(LearnresSVorder)
colnames(LearnresSVorder)

Learn_sig_genes <- as.data.frame(LearnresSVorder)

#colnames(LearnresSVorder)[which(names(LearnresSVorder) == "rownames")] <- "FBgn"

#DE genes for Memory
load(file = "../Data/MemresSVOrder.Rda")
str(MemresSVorder)
MemresSVorder$FBgn <- rownames(MemresSVorder)
colnames(MemresSVorder)

Mem_sig_genes <- as.data.frame(MemresSVorder)

#DE genes thermal Tolerance 
load(file="../Data/TTlrt_inter.Rda")
str(TTlrt_inter)
TTlrt_inter$FBgn <- rownames(TTlrt_inter)
colnames(TTlrt_inter)

ThermTol_sig_genes <- as.data.frame(TTlrt_inter)


#gene list form Fly Base
gene_map_table <-read.table(file = "../Data/gene_map_table_fb_2015_03.tsv.gz", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
str(gene_map_table)
colnames(gene_map_table) <- c('gname','FBgn','v3','cyt','pos')

###
####
#is this truncating early? file has 244185 lines total

gene_map_table$chr <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,1]
temp.s1 <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,2] %>% 
  str_split(fixed(".."),simplify=TRUE)
gene_map_table$startp <- temp.s1[,1]

gene_map_table$stopp <- str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1]

str(gene_map_table)


#gene_map_table <- as.numeric(gene_map_table$startp)

#merge DE genes with gene list from Fly base (merge by FBgn num)
Learn_gene_list_merged <- merge(Learn_sig_genes, gene_map_table, by="FBgn")
Mem_gene_list_merged <- merge(Mem_sig_genes, gene_map_table, by="FBgn")
ThermTol_gene_list_merged <- merge(ThermTol_sig_genes, gene_map_table, by="FBgn")

str(Learn_gene_list_merged)
colnames(Learn_gene_list_merged)


#select specific columns

Learn_gene_list_sub <- Learn_gene_list_merged[c(1,3,6,7,12:14)]
colnames(Learn_gene_list_sub)

Mem_gene_list_sub <- Mem_gene_list_merged[c(1,3,6,7,12:14)]
colnames(Mem_gene_list_sub)

ThermTol_gene_list_sub <- ThermTol_gene_list_merged[c(1,3,6,7,12:14)]
colnames(ThermTol_gene_list_sub)

#filter dataset to narrow down options
#by chromosome number 
TT_genes_3Rpeak <- ThermTol_gene_list_sub %>% filter(chr %in% c("3R"))
str(TT_genes_3Rpeak)

#by startp number
TT_genes_startp <-TT_genes_3Rpeak %>% filter(startp %in% 20000000:21100000)
str(TT_genes_startp)


#TT_genes_FBgn <- TT_genes_startp %>% filter(FBgn %in% c("FBgn0001217"))
#str(TT_genes_FBgn)
