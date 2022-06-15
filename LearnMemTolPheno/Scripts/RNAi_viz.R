library(tidyverse)
library(readxl)
library(viridis)
library(cowplot)

theme_set(theme_cowplot())
source(file="../../Functions/ggplot_theme.R")


ThermDat1 <- read.csv(file="../ProcessedData/RNAi1_tracked_data.csv")
str(ThermDat1)
ThermDat1[1:5,]
ThermDat1$genotype  <- as.character(ThermDat1$genotype)


ThermDat2 <- read.csv(file="../ProcessedData/RNAi2_tracked_data.csv")
str(ThermDat2)
ThermDat2[1:5,]
ThermDat2$genotype  <- as.character(ThermDat2$genotype)

#add a col set id 
ThermDat1$set_id<- paste(ThermDat1$genotype, ThermDat1$cross, sep="_")
ThermDat2$set_id<- paste(ThermDat2$genotype, ThermDat2$cross, sep="_")

#fix errors
ThermDat1[which(ThermDat1$set_id=="60071_25752"), 'cross'] <-25750
ThermDat1[which(ThermDat1$set_id=="33545_3954"), 'genotype'] <-33646
ThermDat1$set_id<- paste(ThermDat1$genotype, ThermDat1$cross, sep="_")


#load background data
set_ids <- read_csv(file="../ProcessedData/SetID_names.csv")

RNAi_backg <- read_excel("../ProcessedData/HeatPlate_autotrack/RNAi_Bckg.xlsx")
colnames(RNAi_backg) <-c("genotype", "backg_id", "Gname")
RNAi_backg <- RNAi_backg[which(is.na(RNAi_backg$genotype)==FALSE),]

#combine data sets
ThermDat1 <- left_join(ThermDat1, set_ids, by="set_id")
ThermDat1 <- left_join(ThermDat1, RNAi_backg[,1:2], by="genotype")

ThermDat2 <- left_join(ThermDat2, set_ids, by="set_id")
ThermDat2 <- left_join(ThermDat2, RNAi_backg[,1:2], by="genotype")


ThermDat_all <- rbind(ThermDat1[,c("file", "Chamber","incapacitation","date","group","rep","genotype","cross","set_id","backg_id","Gname")], 
                      ThermDat2[,c("file", "Chamber","incapacitation","date","group","rep","genotype","cross","set_id","backg_id","Gname")])


ThermDat_all <- ThermDat_all[which(is.na(ThermDat_all$incapacitation)==FALSE),]
ThermDat_all <- ThermDat_all %>% filter(incapacitation > 5)

write_csv(ThermDat_all, file="../ProcessedData/Combined_RNAi_tracked.csv")

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



attP2 <- ThermDat_all %>%
  filter(backg_id=="attP2")

attP2 <- attP2 %>% 
  mutate(Gname = fct_relevel(Gname, "inbred_attP2","attP2_neuro_1","attP2_pan"))

mu_pan_2 <- mean(attP2[attP2$Gname=="attP2_pan", "incapacitation"])
mu_neuro_2 <- mean(attP2[attP2$Gname=="attP2_neuro_1", "incapacitation"])


p1 <- ggplot(attP2, aes(Gname, incapacitation, color=Type)) +
  geom_point(position = position_jitter(width = 0.4), alpha=1/2) +
  stat_summary(fun = mean, geom = "point",shape=16, size = 2, color = "coral") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.05,
               color = "coral", size = 0.7) +
  geom_hline(yintercept=mu_pan_2, color=viridis(4)[1],lty=2) +
  geom_hline(yintercept=mu_neuro_2, color=viridis(4)[2], lty=2) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("id")+
  scale_color_viridis(discrete = TRUE) +
  my_theme
p1


attP40 <- ThermDat_all %>%
  filter(backg_id=="attP40")

attP40 <- attP40 %>% 
  mutate(Gname = fct_relevel(Gname, "inbred_attP40","attP40_neuro_1","attP40_pan"))


mu_pan_40 <- mean(attP40[attP40$Gname=="attP40_pan", "incapacitation"])
mu_neuro_40 <- mean(attP40[attP40$Gname=="attP40_neuro_1", "incapacitation"])


p2 <- ggplot(attP40, aes(Gname, incapacitation, color=Type)) +
  geom_point(position = position_jitter(width = 0.4), alpha=1/2) +
  stat_summary(fun = mean, geom = "point",shape=16, size = 2, color = "coral") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.05,
               color = "coral", size = 0.7) +
  geom_hline(yintercept=mu_pan_40, color=viridis(4)[1],lty=2) +
  geom_hline(yintercept=mu_neuro_40, color=viridis(4)[2], lty=2) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("id")+
  scale_color_viridis(discrete = TRUE) +
  my_theme
p2

pall <- plot_grid(p1,p2, nrow=2,  labels = c("a.","b."))

ggsave(pall, filename = "../Plots/RNAi_phenos_all.pdf", width=6.5,height=9)