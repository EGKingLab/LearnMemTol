
slideThermo <- function(xx, tt)
{
  st <- 1
  steps <- 1
  width <- 60
  tol <- 20
  set <- 25
  # disp <- 20
  #dtol <-10
  starts <- seq(from=1, to=(max(tt)-(width-1)),by=1)
  ends <- seq(from = width, to = max(tt), by = 1)
  counter <- 0
  while((set > tol) & (counter <= length(starts)))
  {
    counter <- counter+1
    if(length(xx[tt >= starts[counter] & tt < ends[counter]]) > 1)
    {set <- var(xx[tt >= starts[counter] & tt < ends[counter]])}
    st <- st+steps
    
  }
  
  if(counter <= length(starts))
  {
  return(min(tt[tt >= starts[counter] & tt < ends[counter]],na.rm=TRUE))  
  }else{
    return(max(tt,na.rm=TRUE))    
  }
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
  
  if(mean(tt1$likelihood) > 0.4) 
  {
    
      Th.set <- data.frame( 'file'=character(length=1),
                            'incapacitation'= numeric(length=1), 
                            'Rvar'= numeric(length=1),
                            stringsAsFactors = FALSE)
      
      lastpos <- tt1[(nrow(tt1)-200):nrow(tt1),]
      lastpos <- lastpos[lastpos$likelihood >= 0.6,]

      Th.set$Rvar <- max(lastpos$x) - min(lastpos$x)
      
      #remove low likelihood positions
      tt1 <- tt1[tt1$likelihood >= 0.6,]
      
      incap.i <- slideThermo(xx=tt1$x,tt=tt1$Second)
      
      
      #Th.set$Rvar <- slideRec(st <- incap.i, xx=tt1$x,tt=tt1$Second)
      Th.set$incapacitation <- incap.i
      Th.set$file <- tt1$id[1]
      
      ThermTol_HeatPlate <- rbind(ThermTol_HeatPlate, Th.set)
    
    
  }
  cat(cc,"\t")
  cc<-cc+1
}

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

proj <- "RNAi_1"
ThermTol_HeatPlate <- readRDS(file=paste0("../ProcessedData/Incap_processed_",proj,".Rds"))
chambers <- read_csv(file="../ProcessedData/HeatPlate_autotrack/Chambers_RNAi_round1.csv")
chambers$FileName <- str_split(chambers$FileName,"_s", simplify=TRUE)[,1]

chambers$file <- paste0(chambers$FileName,"_",chambers$Chamber)

RNAi1 <- left_join(chambers,ThermTol_HeatPlate, by="file")
RNAi1 <- subset(RNAi1, Chamber != 25)
RNAi1 <- subset(RNAi1, FileName != "2021-04-06_RNAi-32871.3954-2_group1")

filt <- RNAi1[which(is.na(RNAi1$incapacitation)),]

proj <- "RNAi_2"
ThermTol_HeatPlate <- readRDS(file=paste0("../ProcessedData/Incap_processed_",proj,".Rds"))

chambers <- read_excel(path="../ProcessedData/HeatPlate_autotrack/RNAi_Val2_Chmb_positions.xlsx")

chambers$file <- paste0(chambers$FileName,"_",chambers$Chamber)

RNAi2 <- left_join(chambers,ThermTol_HeatPlate, by="file")
RNAi2 <- subset(RNAi2, Chamber != 25)

filt <- RNAi2[which(is.na(RNAi2$incapacitation)),]

write_csv(RNAi1, file="../ProcessedData/RNAi1_tracked_data.csv")
write_csv(RNAi2, file="../ProcessedData/RNAi2_tracked_data.csv")



#check for problem tracking

plot(RNAi1$Rvar, RNAi1$incapacitation)

plot(RNAi2$Rvar, RNAi2$incapacitation)


R1_c <- subset(RNAi1, Rvar >75)
R2_c <- subset(RNAi2, Rvar >75)

proj <- "RNAi_1"
TTdat <- readRDS(file=paste0("../ProcessedData/Combined_tracks_",proj,".Rds"))
TTdat <- as.data.table(TTdat)
ThermTol_HeatPlate <- readRDS(file=paste0("../ProcessedData/Incap_processed_",proj,".Rds"))

pdf(file="../Plots/RNA1_score.pdf",width=8, height=12)
par(mfrow=c(3,2))

for(ff in R1_c$file)
{
  tt1 <- TTdat[TTdat$id==ff,]
  tt1 <- tt1[tt1$likelihood >= 0.6,]
  plot(tt1$Second, tt1$x, main=ff)
  abline(v=ThermTol_HeatPlate[ThermTol_HeatPlate$file==ff,'incapacitation'])
}

dev.off()

proj <- "RNAi_2"
TTdat <- readRDS(file=paste0("../ProcessedData/Combined_tracks_",proj,".Rds"))
TTdat <- as.data.table(TTdat)

pdf(file="../Plots/RNA2_score.pdf",width=8, height=12)
par(mfrow=c(3,2))

for(ff in R2_c$file)
{
tt1 <- TTdat[TTdat$id==ff,]
tt1 <- tt1[tt1$likelihood >= 0.6,]
plot(tt1$Second, tt1$x, main=ff)
abline(v=ThermTol_HeatPlate[ThermTol_HeatPlate$file==ff,'incapacitation'])
}

dev.off()


#check random set
pdf(file="../Plots/rand_score.pdf",width=8, height=12)
par(mfrow=c(3,2))

for(ff in sample(c(RNAi1$file,RNAi2$file),6*4))
{
  tt1 <- TTdat[TTdat$id==ff,]
  tt1 <- tt1[tt1$likelihood >= 0.6,]
  plot(tt1$Second, tt1$x, main=ff)
  abline(v=ThermTol_HeatPlate[ThermTol_HeatPlate$file==ff,'incapacitation'])
}

dev.off()


