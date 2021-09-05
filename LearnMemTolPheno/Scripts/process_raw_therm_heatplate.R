
slideThermo <- function(xx, tt)
{
  st <- 1
  steps <- 1
  width <- 60
  tol <- 15
  set <- 20
  # disp <- 20
  #dtol <-10
  starts <- seq(from=1, to=(max(tt)-(width-1)),by=1)
  ends <- seq(from = width, to = max(tt), by = 1)
  counter <- 1
  while((set > tol) & (counter <= length(starts)))
  {
    if(length(xx[tt >= starts[counter] & tt < ends[counter]]) > 1)
    {set <- var(xx[tt >= starts[counter] & tt < ends[counter]])}
    st <- st+steps
    counter <- counter+1
  }
  
  return(counter)  
}

slideRec <- function(st, xx, tt)
{
  
  steps <- 1 
  width <- 60
  
  iis<-seq(from = tt[st], to = (max(tt)-width), by=steps)
  vvs <- numeric(length(iis))
  for(ii in 1:length(iis))
  {
    vvs[ii] <- var(xx[tt >= iis[ii] & tt < (iis[ii] + width)])
  }
  
  return(max(vvs))  
}



library(tidyverse)
library(data.table)
library(readxl)



#set as data frame ###change 
ThermTol_HeatPlate<-data.frame('file'=character(length=0),
                               'incapacitation'= numeric(length=0), 
                               'Rvar'= numeric(length=0),
                               stringsAsFactors = FALSE)

proj <- "RNAi_2"

TTdat <- readRDS(file=paste0("../ProcessedData/Combined_tracks_",proj,".Rds"))
TTdat <- as.data.table(TTdat)

ffs <- unique(TTdat$id)
cc <- 1

for(ff in ffs) {
  tt1 <- TTdat[TTdat$id==ff,]
  
  if(mean(tt1$likelihood) > 0.65) 
  {
    
      Th.set <- data.frame( 'file'=character(length=1),
                            'incapacitation'= numeric(length=1), 
                            'Rvar'= numeric(length=1),
                            stringsAsFactors = FALSE)
      
      #remove low likelihood positions
      tt1 <- tt1[tt1$likelihood >= 0.6,]
      
      incap.i <- slideThermo(xx=tt1$x,tt=tt1$Second)
      
      Th.set$Rvar <- slideRec(st <- incap.i, xx=tt1$x,tt=tt1$Second)
      Th.set$incapacitation <- tt1$Second[incap.i]
      Th.set$file <- tt1$id[1]
      
      ThermTol_HeatPlate <- rbind(ThermTol_HeatPlate, Th.set)
    
    
  }
  cat(cc,"\t")
  cc<-cc+1
}

#Add 07-2 rerun
TT72 <- readRDS(file="../ProcessedData/2021-07-02b_Fly_tracks.Rds")
TT72 <- bind_rows(TT72)

ffs <- unique(TT72$id)

for(ff in ffs) {
  tt1 <- TT72[TT72$id==ff,]
  
if(mean(tt1$likelihood) > 0.65) 
{
  
  Th.set <- data.frame( 'file'=character(length=1),
                        'incapacitation'= numeric(length=1), 
                        'Rvar'= numeric(length=1),
                        stringsAsFactors = FALSE)
  
  #remove low likelihood positions
  tt1 <- tt1[tt1$likelihood >= 0.6,]
  
  incap.i <- slideThermo(xx=tt1$x,tt=tt1$Second)
  
  #Th.set$Rvar <- slideRec(st <- incap.i, xx=tt1$x,tt=tt1$Second)
  Th.set$incapacitation <- tt1$Second[incap.i]
  Th.set$file <- tt1$id[1]
  
  ThermTol_HeatPlate <- rbind(ThermTol_HeatPlate, Th.set)
  
}
}
##07-2

ss <- ThermTol_HeatPlate$file %>%
  str_split("_",simplify=TRUE)

ThermTol_HeatPlate$date <- ss[,1]
ThermTol_HeatPlate$group <- ss[,3]
ThermTol_HeatPlate$chamber <- ss[,4]


ss2 <- ss[,2] %>%
  str_split("-",simplify=TRUE)

ThermTol_HeatPlate$rep <- ss2[,3]

ss3 <- ss2[,2] %>%
  str_split(fixed("."),simplify=TRUE)

ThermTol_HeatPlate$genotype <- ss3[,1]
ThermTol_HeatPlate$cross <- ss3[,2]


#saveRDS(ThermTol_HeatPlate, file=paste0("../ProcessedData/Incap_processed_",proj,".Rds"))


chambers <- read_csv(file="../ProcessedData/HeatPlate_autotrack/Chambers_RNAi_round1.csv")
chambers$FileName <- str_split(chambers$FileName,"_s", simplify=TRUE)[,1]

chambers$file <- paste0(chambers$FileName,"_",chambers$Chamber)

RNAi1 <- left_join(chambers,ThermTol_HeatPlate, by="file")
RNAi1 <- subset(RNAi1, Chamber != 25)
RNAi1 <- subset(RNAi1, FileName != "2021-04-06_RNAi-32871.3954-2_group1")

filt <- RNAi1[which(is.na(RNAi1$incapacitation)),]


chambers <- read_excel(path="../ProcessedData/HeatPlate_autotrack/RNAi_Val2_Chamber_positions.xlsx")

chambers$file <- paste0(chambers$FileName,"_",chambers$Chamber)

RNAi2 <- left_join(chambers,ThermTol_HeatPlate, by="file")
RNAi2 <- subset(RNAi2, Chamber != 25)

filt <- RNAi2[which(is.na(RNAi2$incapacitation)),]

write_csv(RNAi1, file="../ProcessedData/RNAi1_tracked_data.csv")
write_csv(RNAi2, file="../ProcessedData/RNAi2_tracked_data.csv")
