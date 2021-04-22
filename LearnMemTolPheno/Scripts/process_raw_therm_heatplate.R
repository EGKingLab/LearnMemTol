library(stringr)


slideThermo <- function(xx, tt)
{
  st <- 0
  steps <- 1
  width <- 50
  tol <- 1
  set <- 20
 # disp <- 20
  #dtol <-10
  
  while((set > tol)) #& (disp > dtol))
  {
    if(st<=(length(xx)-width)){
      set <- var(xx[st:(st+width)])
      #disp <- max(abs(diff(xx[st:(st+width)],lag=100)))
      #cat(subR$time[st],"\t",subR$time[st+300],"\n")
      st <- st+steps
    }else{
      set<-0
      st<-length(xx)+steps
    }
  }
  
  return(tt[(st)])  
}

slideRec <- function(st, xx, tt)
{
  
  steps <- 1 
  width <- 50
  
  iis<-seq(from = st, to = (length(xx) - width), by=steps)
  vvs <- numeric(length(iis))
  for(ii in 1:length(iis))
  {
    st <- iis[ii]
    vvs[ii] <- var(xx[st:(st+width)])
  }
  
  return(max(vvs))  
}

#set as data frame ###change 
ThermTol_HeatPlate<-data.frame('patRIL'= numeric(length=1) , 'chamber'= numeric(length=1), 
                      'incapacitation'= numeric(length=1), 
                      'Rvar'= numeric(length=1),
                      'group'=character(length=1),
                      'date' =character(length=1),
                      'file'=character(length=1), stringsAsFactors = FALSE)

#basef <- "/home/pwilliams/DSPR/RawData/Founders_Incapacitation_/"
#folds <-list.files(basef)

#fcheck <- data.frame("file"=character(length=0),"rows"=numeric(length=0), stringsAsFactors=FALSE)


#for(kk in folds)
#{
#  fiset<-list.files(file.path(basef,kk,fsep=""),pattern=".asc$")
 # RR <- substr(kk,1,5)
 # for(gg in fiset)
 # {
  
 therm.set<-read.csv(file="../ProcessedData/HeatPlate_autotrack/2021-01-21_VAL-11465-2_group1_9DLC_resnet50_fly_trackerNov28shuffle1_1030000.csv", skip = 2)

 #####might want to check for low likelihood
  #min(therm.set$likelihood)  
 
 
 #########################################################################
 
 RNAi.therm.set <- read.csv(file="../ProcessedData/HeatPlate_autotrack/2021-01-21_VAL-11465-2_group1_9DLC_resnet50_fly_trackerNov28shuffle1_1030000.csv",
                            skip = 6) %>% 
   str_split(pattern = "_", simplify = TRUE) %>% 
   as_tibble() %>% 
   filter(V5 == ".csv")
 
 colnames(fname) <- c('date','genotype','group','v1','v2','v3','v4' ,'v5') #we need to add rep column, split on -2
 
 RNAi.therm.set$genotype <- str_split(RNAi.therm.set$genotype, ":",simplify=TRUE)
 RNAi.therm.set$genotype <- as.numeric(temp.s1[,1])
 
 RNAi.therm.set$v5 <- as.numeric(str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1])
 
 
 
 #splitting file name into different columns 

 fname <- "2021-01-21_RNAi-32871.51941-2_group1_9DLC_resnet50_fly_trackerNov28shuffle1_1030000.csv"
 
 fname<- str_split(fname, pattern = "_", simplify = TRUE)
 colnames(fname) <- c('date','genotype','group','v1','v2','v3','v4' ,'v5')
 
 
 
 fname<- str_split_fixed(fname, pattern =  "-2_",n =1)
 
 fnamet$genotype <- as.numeric(temp.s1[,1])
 
 RNAi.therm.set$v5 <- as.numeric(str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1])
 #########
 ##############################################################
 
 
 incap <- slideThermo(xx=therm.set$x, tt = therm.set$coords)
 

       if(incap >= length(therm.set$x) ){
        res <- NA 
      }else{
        res<-slideRec(incap, xx=therm.set$x, tt = therm.set$coords )
      }
  
  ThermTol_HeatPlate[,'incapacitation'] <-incap  
  