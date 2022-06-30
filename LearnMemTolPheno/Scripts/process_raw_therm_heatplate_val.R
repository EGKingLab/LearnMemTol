
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




proj <- "VAL"

TTdat <- readRDS(file=paste0("../ProcessedData/HeatPlate_autotrack/Combined_tracks_",proj,".Rds"))
extdat <- readRDS(file="../ProcessedData/2021-02-15_VAl-12243-2_group1_Fly_tracks.Rds")
TTdat <- rbind(TTdat,extdat)
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

#Date is fine
ThermTol_HeatPlate$date <- ss[,1]
ThermTol_HeatPlate$genotype <- NA
ThermTol_HeatPlate$group <- NA
ThermTol_HeatPlate$chamber <- NA



#multiple separators used :(
#split up 
ww_a <- which(ss[,2]=="VAL")
ww_b <- which(ss[,2]!="VAL")

ss_extra <- str_split(ss[ww_b,2],"-", simplify=TRUE)
unique(ss_extra[,1])
unique(ss_extra[,2])
unique(ss_extra[,3])
unique(ss_extra[,4])

length(ww_b)
length(ss_extra[,2])

ThermTol_HeatPlate$genotype[ww_b] <- ss_extra[,2]
ThermTol_HeatPlate$group[ww_b] <- ss_extra[,3]

ss_a_geno <- str_split(ss[ww_a,3],"-",simplify=TRUE)
ThermTol_HeatPlate$genotype[ww_a] <- ss_a_geno[,1]
ThermTol_HeatPlate$group[ww_a] <- ss_a_geno[,2]

ss_chamb <- ThermTol_HeatPlate$file %>% 
  str_sub(-3) %>%
  str_split("_", simplify=TRUE)

ThermTol_HeatPlate$chamber <- ss_chamb[,2]


chambers <- read_excel(path="../ProcessedData/HeatPlate_autotrack/VAL_alldata.xlsx")

chambers_f <- str_split(chambers$filename,"_segmentation", simplify = TRUE)
chambers$file <- paste0(chambers_f[,1],"_",chambers$chamber)

alldat <- left_join(chambers,ThermTol_HeatPlate, by="file")
which(alldat$chamber.x==25)
all.equal(alldat$chamber.x, as.numeric(alldat$chamber.y))

alldat <- alldat[is.na(alldat$incapacitation)==FALSE,]

saveRDS(alldat, file=paste0("../ProcessedData/Incap_processed_",proj,".Rds"))


filt <- alldat[which(is.na(alldat$incapacitation)),]

write_csv(filt, file="../ProcessedData/Val_nomatch.csv")


alldat_src <- full_join(chambers,ThermTol_HeatPlate, by="file")

inlist <- alldat_src[which(is.na(alldat_src$chamber.y)),]
indat <- alldat_src[which(is.na(alldat_src$chamber.x)),]
indat <- subset(indat, chamber.y != 25)

fname_inlist <- table(inlist$filename)

inlist[inlist$filename==names(fname_inlist[15]),"file"]

#2021-01-21_VAL-11465-2_group1 - massive highlight, need to drop
#2021-02-03_VAL-11239-2_group1 & 2021-02-18_VAL-A3_group1-2 - very bright, need to drop
#2021-02-11_VAL-11037-2_group1 & 2021-02-11_VAL-12043-2_group1  & 2021-02-15_VAL-11305-2_group1- camera off

#2021-02-09_VAL-12097-2_group2 - low probs for some reason - bright?
#2021-02-15_VAl-12243-2 - dropped because of project name - will read in separately
