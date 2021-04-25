

slideThermo <- function(xx, tt)
{
  st <- 1
  steps <- 1
  width <- 60
  tol <- 1
  set <- 20
  # disp <- 20
  #dtol <-10
  starts <- seq(from=1, to=(max(tt)-(width-1)),by=1)
  ends <- seq(from = width, to = max(tt), by = 1)
  counter <- 1
  while((set > tol) & (counter <= length(starts)))
  {
    set <- var(xx[tt >= starts[counter] & tt < ends[counter]])
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


#set as data frame ###change 
ThermTol_HeatPlate<-data.frame('file'=character(length=0),
                               'incapacitation'= numeric(length=0), 
                               'Rvar'= numeric(length=0),
                               stringsAsFactors = FALSE)


TTdat <- readRDS(file="Combined_tracks.Rds")
TTdat <- as.data.table(TTdat)

ffs <- unique(TTdat$id)
cc <- 1

for(ff in ffs) {
  tt1 <- TTdat[TTdat$id==ff,]
  if(mean(tt1$likelihood) > 0.5 & length(grep("RNAi",tt1$id[1])) > 0) 
  {
    
    Th.set <- data.frame( 'file'=character(length=1),
                          'incapacitation'= numeric(length=1), 
                          'Rvar'= numeric(length=1),
                          stringsAsFactors = FALSE)
    
    incap.i <- slideThermo(xx=tt1$x,tt=tt1$Second)
    
    Th.set$Rvar <- slideRec(st <- incap.i, xx=tt1$x,tt=tt1$Second)
    Th.set$incapacitation <- tt1$Second[incap.i]
    Th.set$file <- tt1$id[1]
    
    ThermTol_HeatPlate <- rbind(ThermTol_HeatPlate, Th.set)
  }
  cat(cc,"\n")
  cc<-cc+1
}

saveRDS(ThermTol_HeatPlate, file="Incap_processed.Rds")

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


saveRDS(ThermTol_HeatPlate, file="Incap_processed.Rds")

