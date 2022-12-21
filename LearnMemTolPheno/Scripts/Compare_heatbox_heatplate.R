library(nlme)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

tt_plate <- readRDS("../ProcessedData/Incap_processed_VAL.Rds")
# looks like 2021-01-28_VAL-11228-2_group1_segmentation did not get to temp - plate elevated??
#will drop 11228 because data only from one day then
#will drop outlier values too

tt_plate <- tt_plate[which(!(tt_plate$genotype %in% "11228")),]
tt_plate <- subset(tt_plate, incapacitation >5 & incapacitation < 400)

genos <- unique(tt_plate$genotype)

tt_box <- read.table("../ProcessedData/ThermalTol_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)

hap_code <- read.table(file="../../LearnMemTolQTL/ProcessedData/HL_thermQTL_lines.txt", sep="\t", header=TRUE)
hap_code$genotype <- as.character(hap_code$patRIL)

plM <- tt_plate %>% group_by(genotype) %>%
  summarise("Mtol" = median(incapacitation))

blM <- tt_box %>% group_by(patRIL) %>%
  summarise("Mtol" = median(incapacitation))
blM$genotype <- as.character(blM$patRIL)

cts <- tt_box %>% group_by(patRIL) %>% tally() 

plM$type <- "plate"
blM$type <- "box"
allM <- rbind(plM, blM[,c("genotype","Mtol","type")])
allM <- left_join(hap_code, allM, by="genotype")

p1 <- ggplot(allM, aes(hard,Mtol, color=type)) +
  geom_point() 
p1

wideAll <- inner_join(plM, blM, by="genotype")
wideAll <- inner_join(hap_code, wideAll, by="genotype")
ww_exclude <- c("11316", "11360","11465","12097","12316")
wideAll <- wideAll[-which(wideAll$genotype %in% ww_exclude),]
cor(wideAll$Mtol.x, wideAll$Mtol.y)

mms <- data.frame("Mtol.x"=c(mean(wideAll[wideAll$hard=="A4","Mtol.x"]),
                             mean(wideAll[wideAll$hard=="A5","Mtol.x"])),
                  "Mtol.y"=c(mean(wideAll[wideAll$hard=="A4","Mtol.y"]),
                             mean(wideAll[wideAll$hard=="A5","Mtol.y"])),
                  "Haplotype"=c("A4","A5","A4","A5"))
sds <- data.frame("Mtol.x"=c(sd(wideAll[wideAll$hard=="A4","Mtol.x"]),
                             sd(wideAll[wideAll$hard=="A5","Mtol.x"])),
                  "Mtol.y"=c(sd(wideAll[wideAll$hard=="A4","Mtol.y"]),
                             sd(wideAll[wideAll$hard=="A5","Mtol.y"])),
                  "Haplotype"=c("A4","A5","A4","A5"))

wideAll$Haplotype <- wideAll$hard
p2 <- ggplot(wideAll, aes(Mtol.x, Mtol.y, color=Haplotype)) +
  geom_point(size = 2) +
  #geom_point(data=mms, size=5, pch=1) +
  #geom_abline(slope=1, intercept=0) +
  geom_label(label=wideAll$genotype) + 
  xlab("Heat Plate Score") +
  ylab("Heat Box Score")
p2


p2 <- ggplot(wideAll, aes(Mtol.y, Mtol.x, color=hard)) +
  geom_point()
p2

ggeno <- unique(tt_plate$genotype)
hist(as.numeric(tt_plate[tt_plate$genotype=="12169","incapacitation",drop=TRUE]),20)
hist(tt_plate$incapacitation, 100)


tester <- tt_plate[tt_plate$genotype==ggeno[23],]
tester %>% group_by(filename) %>% summarise("M"= mean(incapacitation))

p3 <- ggplot(tt_plate, aes(genotype, incapacitation, color=date)) +
  geom_point(alpha=1/2) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3

#anova with date? difference 

dds <- tt_plate %>% group_by(date) %>% summarise("M"=mean(incapacitation))

aa <- aov(incapacitation ~ genotype + date, data=tt_plate)
summary(aa)


tt_plate_hap <- inner_join(hap_code[,c("genotype","hard")],tt_plate, by="genotype")
aa <- lme(incapacitation ~ hard, random = ~ 1|genotype, data=tt_plate_hap)
summary(aa)

date_sum <- tt_plate %>%
  group_by(date) %>%
  summarise(date_mean = mean(incapacitation))

date_sum %>%
  ggplot(aes(date, date_mean)) +
  geom_point()

p4 <- ggplot(tt_plate, aes(date, incapacitation, color=genotype)) +
  geom_point(alpha=1/2) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p4

