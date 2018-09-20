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
pp.set[,c('chr','Ppos')] <- merge(pp.set[,c('chr','Ppos')], coord.table, by.x=c('chr','Ppos'), by.y=c('R5chr','R5pos'))[,c('R6chr','R6pos')]
pp.set[,c('up')] <- merge(pp.set[,c('chr','up')], coord.table, by.x=c('chr','up'), by.y=c('R5chr','R5pos'))[,c('R6pos')]
pp.set[,c('lp')] <- merge(pp.set[,c('chr','lp')], coord.table, by.x=c('chr','lp'), by.y=c('R5chr','R5pos'))[,c('R6pos')]
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
#these paths are incorrect
load(file = "../Data/LearnresSVOrder.Rda")
str(LearnresSVorder)

#DE genes for Memory
load(file = "../Data/MemresSVOrder.Rda")
str(MemresSVorder)


#gene list form Fly Base
gene_map_table <-read.table(file = "/home/pwilliams/DSPR/RawData/gene_map_table_fb_2015_03.tsv", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
str(gene_map_table)
colnames(gene_map_table) <- c('gname','FBgn','v3','cyt','pos')
#is this truncating early? file has 244185 lines total

gene_map_table$chr <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,1]
temp.s1 <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,2] %>% 
  str_split(fixed(".."),simplify=TRUE)
gene_map_table$startp <- temp.s1[,1]

gene_map_table$stopp <- str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1]



#sort gene_map_table dataset


gene_map_table_sort <- data.frame('current_symbol ')character(length=length(gene_map_table)),
                                 'recombination_loc'=numeric(length=length(gene_map_table)), 
                                 'cytogenetic_loc'=numeric(length=length(gene_map_table)),
                                'sequence_loc'=numeric(length=length(Incappeaks)), stringsAsFactors = FALSE)

#hack
gene_map_table <- rbind(gene_map_table_sort,c('2R','21410000','NA','NA'))

for(kk in 1:length(gene_map_table))
  
  
#change colnames to match up with the other dataset


#merge DE genes with gene list from Fly base
Learn_gene_list_merged <- merge(LearnresSVorder, gene_map_table, by="gene_name")
Mem_gene_list_merged <- merge(MemresSVorder, gene_map_table, by="gene_name")








