library(tidyverse)
library(readxl)
library(viridis)
library(cowplot)

theme_set(theme_cowplot())
source(file="../../Functions/ggplot_theme.R")

ThermDat_all <- read_csv(file="../ProcessedData/Combined_RNAi_tracked.csv")
ThermDat_all$Gene <- str_split(ThermDat_all$Gname, "_", simplify=TRUE)[,1]

ThermDat_all$Type <- ThermDat_all$cross
ThermDat_all$Type[ThermDat_all$cross==3954] <- "all_cells"
ThermDat_all$Type[ThermDat_all$cross==25750] <- "neuro_1"
ThermDat_all$Type[ThermDat_all$cross==51941] <- "neuro_2"
ThermDat_all$Type[is.na(ThermDat_all$cross)] <- "inbred"

ThermDat_all <- ThermDat_all %>% filter(Type %in% c("all_cells", "neuro_1"))


gene_map_table <- read_lines(file = "../../LearnMemTolQTL/ProcessedData/gene_map_table_fb_2018_04.tsv",
                             skip = 6) %>% 
  str_split(pattern = "\t", simplify = TRUE) %>% 
  as_tibble() %>% 
  filter(V1 == "Dmel")

colnames(gene_map_table) <- c('spp','gname','FBgn','v3','cyt','pos')

sm.tab <- inner_join(ThermDat_all, gene_map_table[,c("gname","FBgn")], by=c("Gene"="FBgn"))
sm.tab$FBgn <- sm.tab$Gene
sm.tab2 <- inner_join(ThermDat_all,gene_map_table[,c("gname","FBgn")], by=c("Gene"="gname"))

#grab background data
nogene <- ThermDat_all[ThermDat_all$Gene %in% c("attP2","attP40"),]
nogene$FBgn <- NA
nogene$gname <- NA

ThermDat_all <- rbind(sm.tab[,colnames(sm.tab2)], sm.tab2)
ThermDat_all <- left_join(ThermDat_all, gene_map_table[,c("gname","FBgn")], by="FBgn")

#assign background id to background crosses

nogene$backg_id[nogene$genotype=="36303"] <- "attP2"
nogene$backg_id[nogene$genotype=="36304"] <- "attP40"
nogene$rnaicol <- "#000000" 

maxn <- ThermDat_all %>% group_by(Gene) %>% summarise("mn"=length(unique(Gname)))
col_r <- viridis(max(maxn[,2])/2)
ThermDat_all$rnaicol <- col_r[1]

maxp_i <- max(ThermDat_all$incapacitation) + 50

plist <- vector(mode="list", length= length(unique(ThermDat_all$Gene)))
counter <- 1

foc.gene.set <- sort(unique(ThermDat_all$Gene))

for(foc.gene in foc.gene.set)
{
  GG <- ThermDat_all %>%
    filter(Gene==foc.gene)
  
  #assign colors to lines
  ggs <- unique(GG$genotype)
  for(ii in 1:length(ggs))
  {
    GG$rnaicol[GG$Gene==foc.gene & GG$genotype==ggs[ii]] <- col_r[ii]
  }
  
  addon <- nogene %>% filter(backg_id %in% unique(GG$backg_id) & Type %in% unique(GG$Type))
  addon$Gene <- foc.gene
  
  gorder <- character(length=0)
  laborder <- character(length=0)
  for(ff in unique(addon$set_id[order(addon$genotype,addon$cross)]))
  {
    aa <- subset(addon, set_id==ff)
    gorder <- c(gorder,ff)
    gorder <- c(gorder,unique(GG$set_id[GG$backg_id==aa$backg_id[1] & GG$cross==aa$cross[1]]))
    laborder <- c(laborder, paste0("BK(",aa$backg_id[1],")"))
    laborder <- c(laborder,unique(GG$genotype[GG$backg_id==aa$backg_id[1] & GG$cross==aa$cross[1]]))
    
  }
  
  GG <- rbind(GG, addon)
  GG$plot_x <- factor(GG$set_id, levels=gorder)  
  
  plist[[counter]] <- ggplot(GG, aes(plot_x, incapacitation, shape=Type)) +
    geom_point(position = position_jitter(width = 0.3), alpha=1/3, color=GG$rnaicol, size=2) +
    stat_summary(fun = mean, geom = "point", size = 2.5, color = "red") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
                 color = "black", linewidth = 0.6) +
    annotate("label", x= length(unique(GG$plot_x))*0.25, 
             y= maxp_i - 10, label=GG$gname[1], size=3)+
    ylim(c(0,maxp_i)) +
    scale_x_discrete(breaks= gorder, labels=laborder) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    xlab("id")+
    #scale_shape_discrete(labels=c("all cells", "all neurons")) +
    theme(legend.position = "none") + 
    my_theme
  
  counter <- counter + 1
}

pp1 <- plot_grid(plotlist=plist[1:8 ], ncol=2, 
                labels=paste0(LETTERS[1:8],"."), 
                label_size = 10)
pp2 <- plot_grid(plotlist=plist[9:14], ncol=2, 
                 labels=paste0(LETTERS[9:14],"."), 
                 label_size = 10)

ggsave(pp1, filename = "../Plots/FigS8_1.pdf", height= 10, width=6.5)
ggsave(pp2, filename = "../Plots/FigS8_2.pdf", height= 7.5, width=6.5)
