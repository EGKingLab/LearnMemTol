library(nlme)
library(tidyverse)
library(cowplot)
library(lubridate)
theme_set(theme_cowplot())
source(file="../../Functions/ggplot_theme.R")

tt_plate <- readRDS("../ProcessedData/Incap_processed_VAL.Rds")
# looks like 2021-01-28_VAL-11228-2_group1_segmentation did not get to temp - plate elevated??
#will drop 11228 because data only from one day then
#will drop outlier values too

tt_plate <- tt_plate[which(!(tt_plate$genotype %in% "11228")),]

min(table(tt_plate$genotype))

genos <- unique(tt_plate$genotype)

tt_box <- read.table("../ProcessedData/ThermalTol_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)

hap_code <- read.table(file="../../LearnMemTolQTL/ProcessedData/HL_thermQTL_lines.txt", sep="\t", header=TRUE)
hap_code$genotype <- as.character(hap_code$patRIL)

plM <- tt_plate %>% group_by(genotype) %>%
  summarise("Mtol" = mean(incapacitation))

ctsp <- tt_plate %>% group_by(genotype) %>% tally() 

blM <- tt_box %>% group_by(patRIL) %>%
  summarise("Mtol" = mean(incapacitation))
blM$genotype <- as.character(blM$patRIL)

ctsb <- tt_box %>% group_by(patRIL) %>% tally() 

plM$type <- "plate"
blM$type <- "box"
allM <- rbind(plM, blM[,c("genotype","Mtol","type")])
allM <- left_join(hap_code, allM, by="genotype")

p0 <- ggplot(allM, aes(hard,Mtol, color=type)) +
  geom_point() 
p0

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
p1 <- ggplot(wideAll, aes(Mtol.x, Mtol.y)) +
  geom_point(size = 2) +
  geom_abline(slope=1, intercept=0) +
  geom_label(label=wideAll$genotype, size=2) +
  xlim(c(60,150)) +
  xlab("Heat Plate Score") +
  ylab("Heat Box Score") +
  my_theme
p1


## Check hand scores

hh <- read_csv(file="../ProcessedData/TT-VAL_groundtruthing_CO_blanksNA.csv")[,1:5]
#fix inconsistant date coding
unique(hh$Genotype)
hh$Genotype[hh$Genotype=="2020-12-4_VAL-12145_group1"] <- "2020-12-04_VAL-12145_group1"
hh$Genotype[hh$Genotype=="2021-1-27_VAL-11374-2_group1"] <- "2021-01-27_VAL-11374-2_group1"
hh$Genotype[hh$Genotype=="2021-2-8_VAL-12191-2_group1"] <- "2021-02-08_VAL-12191-2_group1"
hh$Genotype[hh$Genotype=="2021-2-3_VAL-12075-2_group1"] <- "2021-02-03_VAL-12075-2_group1"
#this last group (12075 had a temp problem- excluded)

tt_hh <- str_split(hh$Incap_time_CO, ":", simplify=TRUE)
hh$handtime <- ms(paste0(tt_hh[,1],":", tt_hh[,2]))

hh$file <- paste0(hh$Genotype, "_", hh$Chamber)
hh$handIncap <- as.period(hh$handtime, unit="seconds")

jj <- left_join(hh, tt_plate, by="file")

#12145 17:57
#11377 18:22
#11374 17:50
#12191 17:12

cor(as.period(jj$incapacitation, unit="seconds"), jj$handIncap, use="complete.obs")

p2 <- jj |>
  drop_na() |>
  ggplot(aes(as.period(incapacitation, unit="seconds"), handIncap, color=genotype)) +
  geom_point() +
  ylab("Hand Score") +
  xlab("Automated Score") +
  my_theme

pp_g <- plot_grid(p1,p2, ncol=2, 
                 labels=c("A.","B."), 
                 label_size = 10, rel_widths = c(1,1.3))

ggsave(pp_g, filename = "../Plots/FigS3.pdf", height= 3, width=6.5)

