library(tidyverse)
library(readxl)
library(viridis)
library(cowplot)

theme_set(theme_cowplot())
source(file="../../Functions/ggplot_theme.R")

ThermDat_all <- read_csv(file="../ProcessedData/Combined_RNAi_tracked.csv")

ThermDat_all$Gene <- str_split(ThermDat_all$Gname, "_", simplify=TRUE)[,1]

ch.g <- unique(ThermDat_all[,c("genotype","set_id","Gname")])
ch.g <- ch.g[order(ch.g$Gname),]

TTmeans <- ThermDat_all %>% 
  group_by(set_id) %>%
  summarise(MeanI = mean(incapacitation))

#assign background id to some inbred lines

ThermDat_all$backg_id[ThermDat_all$genotype=="36303"] <- "attP2"
ThermDat_all$backg_id[ThermDat_all$genotype=="36304"] <- "attP40"

ThermDat_all$Type <- ThermDat_all$cross
ThermDat_all$Type[ThermDat_all$cross==3954] <- "all cells"
ThermDat_all$Type[ThermDat_all$cross==25750] <- "neuro_1"
ThermDat_all$Type[ThermDat_all$cross==51941] <- "neuro_2"
ThermDat_all$Type[is.na(ThermDat_all$cross)] <- "inbred"

ThermDat_all$Type <- as.factor(ThermDat_all$Type)

ThermDat_all <- ThermDat_all %>% 
  mutate(Type = fct_relevel(Type, "inbred", after=3))

unique(ThermDat_all$Gname)

for(foc.gene in unique(ThermDat_all$Gene))
{
#foc.gene <- ThermDat_all$Gene[1]

GG <- ThermDat_all %>%
  filter(Gene==foc.gene)

bb <- c(paste0(GG$backg_id[1],"_pan"),paste0(GG$backg_id[1],"_neuro_1"))

BB <- ThermDat_all %>%
  filter(Gname %in% bb)

GG <- rbind(GG,BB)

GG <- GG %>% 
  filter(cross != 51941)

GG <- GG %>% 
  mutate(Gname = fct_relevel(Gname, bb[1],bb[2]))

mu_pan <- mean(GG[GG$Gname==bb[1], "incapacitation", drop=TRUE])
mu_neuro <- mean(GG[GG$Gname==bb[2], "incapacitation", drop=TRUE])

p1 <- ggplot(GG, aes(Gname, incapacitation, color=Type)) +
  geom_point(position = position_jitter(width = 0.1), alpha=1/2) +
  stat_summary(fun = mean, geom = "point",shape=16, size = 2, color = "coral") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.05,
               color = "coral", size = 0.7) +
  geom_hline(yintercept=mu_pan, color=viridis(2)[1],lty=2) +
  geom_hline(yintercept=mu_neuro, color=viridis(2)[2], lty=2) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("id")+
  scale_color_viridis(discrete = TRUE) +
  my_theme

ggsave(p1,filename = paste0("../Plots/Ind_RNAi_",foc.gene,".pdf"),width=6, height = 5)

cat(foc.gene, "\n")
}


pch <- ggplot(GG, aes(Gname, incapacitation, color=date)) +
  geom_point(position = position_jitter(width = 0.1), alpha=1/2)
