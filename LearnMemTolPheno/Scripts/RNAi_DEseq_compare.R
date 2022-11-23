library(tidyverse)
library(readxl)
library(viridis)
library(cowplot)
library(DESeq2)
library(ggrepel)
library(paletteer)

theme_set(theme_cowplot())
source(file="../../Functions/ggplot_theme.R")

ThermDat_all <- read_csv(file="../ProcessedData/Combined_RNAi_tracked.csv")

ThermDat_all$Gene <- str_split(ThermDat_all$Gname, "_", simplify=TRUE)[,1]

ch.g <- unique(ThermDat_all[,c("genotype","set_id","Gname")])
ch.g <- ch.g[order(ch.g$Gname),]

ThermDat_all$backg_id[ThermDat_all$genotype=="36303"] <- "attP2"
ThermDat_all$backg_id[ThermDat_all$genotype=="36304"] <- "attP40"

ThermDat_all$Type <- ThermDat_all$cross
ThermDat_all$Type[ThermDat_all$cross==3954] <- "all cells"
ThermDat_all$Type[ThermDat_all$cross==25750] <- "neuro_1"
ThermDat_all$Type[ThermDat_all$cross==51941] <- "neuro_2"
ThermDat_all$Type[is.na(ThermDat_all$cross)] <- "inbred"

ThermDat_all$Type <- as.factor(ThermDat_all$Type)

TTmeans <- ThermDat_all %>% 
  group_by(set_id,genotype,Type,backg_id,Gene) %>%
  summarise(MeanI = mean(incapacitation))

att2_mu_all <- as.numeric(TTmeans[TTmeans$genotype==36303 & TTmeans$Type=="all cells", "MeanI", drop=TRUE])
att2_mu_n1 <- as.numeric(TTmeans[TTmeans$genotype==36303 & TTmeans$Type=="neuro_1", "MeanI", drop=TRUE])

att4_mu_all <- as.numeric(TTmeans[TTmeans$genotype==36304 & TTmeans$Type=="all cells", "MeanI", drop=TRUE])
att4_mu_n1 <- as.numeric(TTmeans[TTmeans$genotype==36304 & TTmeans$Type=="neuro_1", "MeanI", drop=TRUE])

TTmeans <- TTmeans %>% filter(!genotype %in% c(36303,36304) & !Type %in% c("inbred","neuro_2"))

TTmeans$Difference <- NA
TTmeans$Difference[TTmeans$backg_id=="attP2" & TTmeans$Type=="all cells"]<- TTmeans$MeanI[TTmeans$backg_id=="attP2" & TTmeans$Type=="all cells"] - att2_mu_all
TTmeans$Difference[TTmeans$backg_id=="attP2" & TTmeans$Type=="neuro_1"]<- TTmeans$MeanI[TTmeans$backg_id=="attP2" & TTmeans$Type=="neuro_1"] - att2_mu_n1
TTmeans$Difference[TTmeans$backg_id=="attP40" & TTmeans$Type=="all cells"]<- TTmeans$MeanI[TTmeans$backg_id=="attP40" & TTmeans$Type=="all cells"] - att4_mu_all
TTmeans$Difference[TTmeans$backg_id=="attP40" & TTmeans$Type=="neuro_1"]<- TTmeans$MeanI[TTmeans$backg_id=="attP40" & TTmeans$Type=="neuro_1"] - att4_mu_n1




load("../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/TTlrt_pool.Rda")
TTlrt_pool$FBgn <- rownames(TTlrt_pool)
TTpool_sig_genes <- as.data.frame(TTlrt_pool)

gene_map_table <- read_lines(file = "../../LearnMemTolQTL/ProcessedData/gene_map_table_fb_2018_04.tsv",
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

load(file ="../../LearnMemTolQTL/ProcessedData/Peaks_wCIs.rda")

lowc <- ci.peak$Tolsqrtvariable[which.max(ci.peak$Tolsqrtvariable$LL), 'lpR6']
highc <- ci.peak$Tolsqrtvariable[which.max(ci.peak$Tolsqrtvariable$LL), 'upR6']

gglist <- subset(gene_map_table, chr=="3R" & startp <= highc & stopp >=lowc)

sort(c(which(TTmeans$Gene %in% gglist$gname), which(TTmeans$Gene %in% gglist$FBgn)))

sm.tab <- inner_join(TTmeans, gglist[,c("gname","FBgn")], by=c("Gene"="FBgn"))
sm.tab$FBgn <- sm.tab$Gene
sm.tab2 <- inner_join(TTmeans, gglist[,c("gname","FBgn")], by=c("Gene"="gname"))

TTmeans <- rbind(sm.tab[,colnames(sm.tab2)], sm.tab2)

alldat <- left_join(TTmeans, TTpool_sig_genes, by="FBgn")
alldat$labelpos <- 75

#make color palette for genes - diverge for each
#order by log fold change
#add fold change value to top with color by sig - make this
#add label to fold change
alldat$GeneID <- factor(alldat$Gene) 
alldat$GeneID <- reorder(alldat$GeneID, alldat$log2FoldChange)

alldat$desig <- "black"
alldat$desig[alldat$padj < 0.05] <- "darkred"
alldat$desig[alldat$padj < 0.01] <- "red"



ggplot(alldat, aes(GeneID,Difference, shape=Type)) +
  geom_point() +
  ylim(c(-130, 75)) +
  geom_vline(xintercept=9.5, color="grey50")+ z
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=seq(1,14), lty=3, color="grey50")+
  geom_label(aes(GeneID, labelpos, label=round(log2FoldChange,2)), color=alldat$desig) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  ylab("Difference")
  
