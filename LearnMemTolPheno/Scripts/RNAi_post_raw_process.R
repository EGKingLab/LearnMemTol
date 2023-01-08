library(tidyverse)
library(readxl)

ThermDat1 <- read.csv(file="../ProcessedData/RNAi1_tracked_data.csv")
str(ThermDat1)
ThermDat1[1:5,]
ThermDat1$genotype  <- as.character(ThermDat1$genotype)


ThermDat2 <- read.csv(file="../ProcessedData/RNAi2_tracked_data.csv")
str(ThermDat2)
ThermDat2[1:5,]
ThermDat2$genotype  <- as.character(ThermDat2$genotype)

#remove duplicate row
ThermDat2 <- ThermDat2[-which(duplicated(ThermDat2, MARGIN=1)),]

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
