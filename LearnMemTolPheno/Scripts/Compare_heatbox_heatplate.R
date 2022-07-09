library(tidyverse)

tt_plate <- readRDS("../ProcessedData/Incap_processed_VAL.Rds")

tt_box <- read.table("../ProcessedData/ThermalTol_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)


plM <- tt_plate %>% group_by(genotype) %>%
  summarise("Mtol" = mean(incapacitation))

blM <- tt_box %>% group_by(patRIL) %>%
  summarise("Mtol" = mean(incapacitation))
blM$genotype <- as.character(blM$patRIL)

cts <- tt_box %>% group_by(patRIL) %>% tally() 

allM <- inner_join(plM, blM, by="genotype")

cor(allM$Mtol.x, allM$Mtol.y)
plot(allM$Mtol.x[-9], allM$Mtol.y[-9])
