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

#add a col set id 
ThermDat1$set_id<- paste(ThermDat1$genotype, ThermDat1$cross, sep="_")
ThermDat2$set_id<- paste(ThermDat2$genotype, ThermDat2$cross, sep="_")

#load background data
RNAi_backg <- read_excel("../ProcessedData/HeatPlate_autotrack/RNAi_Bckg.xlsx")
str(RNAi_backg)
RNAi_backg[1:5,]

#assign column names
colnames(RNAi_backg) <-c("genotype", "backg_id", "Gname")

#combine data sets
ThermDat1 <- left_join(ThermDat1, RNAi_backg, by="genotype")
ThermDat2 <- left_join(ThermDat2, RNAi_backg, by="genotype")

ThermDat_all <- rbind(ThermDat1[,c("file", "Chamber","incapacitation","date","group","rep","genotype","cross","set_id","backg_id","Gname")], 
                      ThermDat2[,c("file", "Chamber","incapacitation","date","group","rep","genotype","cross","set_id","backg_id","Gname")])


#write_csv(ThermDat_all, file="../ProcessedData/Combined_RNAi_tracked.csv")

TTmeans <- ThermDat_all %>% 
  group_by(set_id) %>%
  summarise(MeanI = mean(incapacitation))

p1 <- ggplot(ThermDat_all, aes(genotype, incapacitation, color=as.character(cross))) +
  geom_point(position = position_jitter(width = 0.4)) 
p1
