slideThermo <- function(xx, tt)
{
  st <- 301
  steps <- 50
  width <- 600
  #tol <- 1
  tol <- 5
  set <- 20
  disp <- 20
  dtol <-10
  
  while((set > tol) & (disp > dtol))
  {
    if(st<=(length(xx)-width)){
      set <- var(xx[st:(st+width)])
      disp <- max(abs(diff(xx[st:(st+width)],lag=100)))
      #cat(subR$time[st],"\t",subR$time[st+300],"\n")
      st <- st+50
    }else{
      set<-0
      st<-length(xx)+50
    }
  }
  
  return(tt[(st-50)])  
}

slideRec <- function(st, xx, tt)
{
  steps <- 50
  width <- 600
  
  iis<-seq(st,5701,by=50)
  vvs <- numeric(length(iis))
  for(ii in 1:length(iis))
  {
    st <- iis[ii]
    vvs[ii] <- var(xx[st:(st+width)])
  }
  
  return(max(vvs))  
}

#set as data frame
ThermTol <-data.frame('patRIL'= numeric(length=0) , 'chamber'= numeric(length=0), 
                      'incapacitation'= numeric(length=0), 
                      'target.Temp.test.A'= numeric(length=0), 
                      'actual.Temp.test.A'= numeric(length=0),
                      'target.Temp.test.B'= numeric(length=0), 
                      'actual.Temp.test.B'= numeric(length=0),
                      'max.actual.Temp.test.A'= numeric(length=0),
                      'max.actual.Temp.test.B'= numeric(length=0),
                      'target.Temp.train.A'= numeric(length=0), 
                      'actual.Temp.train.A'= numeric(length=0),
                      'target.Temp.train.B'= numeric(length=0), 
                      'actual.Temp.train.B'= numeric(length=0),
                      'max.actual.Temp.train.A'= numeric(length=0),
                      'max.actual.Temp.train.B'= numeric(length=0),
                      'Rvar'= numeric(length=0),
                      'ThermAct'= numeric(length=0),
                      'group'=character(length=0),
                      'file'=character(length=0), stringsAsFactors = FALSE)

basef <- "/home/pwilliams/DSPR/RawData/Main_Incapacitation/"
folds <-list.files(basef)

fcheck <- data.frame("file"=character(length=0),"rows"=numeric(length=0), stringsAsFactors=FALSE)


for(kk in folds)
{
  fiset<-list.files(file.path(basef,kk,fsep=""),pattern=".asc$")
  RR <- substr(kk,1,5)
  for(gg in fiset)
  {
    therm.set<-read.table(file.path(basef,kk,gg),sep="\t", header=TRUE, skip=29)
    fcheck <- rbind(fcheck,data.frame("file"=file.path(basef,kk,gg),"rows"=nrow(therm.set), stringsAsFactors=FALSE))
    all.dat <-data.frame('patRIL'= rep(RR,16) , 'chamber'= numeric(length=16), 
                         'incapacitation'=numeric(length=16), 
                         'target.Temp.test.A'=numeric(length=16), 
                         'actual.Temp.test.A'=numeric(length=16),
                         'target.Temp.test.B'=numeric(length=16), 
                         'actual.Temp.test.B'=numeric(length=16),
                         'max.actual.Temp.test.A'= numeric(length=16),
                         'max.actual.Temp.test.B'= numeric(length=16),
                         'target.Temp.train.A'=numeric(length=16), 
                         'actual.Temp.train.A'=numeric(length=16),
                         'target.Temp.train.B'=numeric(length=16), 
                         'actual.Temp.train.B'=numeric(length=16),
                         'max.actual.Temp.train.A'= numeric(length=16),
                         'max.actual.Temp.train.B'= numeric(length=16),
                         'Rvar'=numeric(length=16),
                         'ThermAct'=numeric(length=16),
                         'group'=character(length=16),
                         'file'=character(length=16),stringsAsFactors = FALSE)
    chambs<-seq(2,16)
    allR<-therm.set[,1:10]
    all.dat[1,'chamber'] <- allR$chamber[1]
    incap <- slideThermo(xx=allR$pos, tt = allR$time)
    all.dat[1,'incapacitation'] <- incap
    all.dat[1, "ThermAct"] <- var(allR$pos) 
    if(incap > 629999){
      all.dat[1,'Rvar'] <- NA 
    }else{
      all.dat[1,'Rvar'] <- slideRec(which(allR$time==incap),xx=allR$pos, tt = allR$time)
    }
    gg.1 <- substring(gg, 6)
    
    all.dat[1,'group'] <- as.numeric(gsub("[^\\d]+", "", gg.1,perl=TRUE))
    
    all.dat[1,'file'] <- gg
    #first set incap done above
    #target temp for accl - time 0 -300 target
    all.dat[1,'target.Temp.test.A'] <-mean(c(therm.set[1:301, 7]))
    all.dat[1,'actual.Temp.test.A'] <-mean(c(therm.set[1:301, 8]))
    all.dat[1,'target.Temp.test.B'] <-mean(c(therm.set[1:301, 9]))
    all.dat[1,'actual.Temp.test.B'] <-mean(c(therm.set[1:301, 10]))
    all.dat[1,'target.Temp.train.A'] <-mean(c(therm.set[302:6301, 7]))
    all.dat[1,'actual.Temp.train.A'] <-mean(c(therm.set[302:6301, 8]))
    all.dat[1,'target.Temp.train.B'] <-mean(c(therm.set[302:6301, 9]))
    all.dat[1,'actual.Temp.train.B'] <-mean(c(therm.set[302:6301, 10]))
    all.dat[1,'max.actual.Temp.test.A'] <-max(c(therm.set[1:301, 8]))
    all.dat[1,'max.actual.Temp.train.A'] <-max(c(therm.set[302:6301, 8]))
    all.dat[1,'max.actual.Temp.test.B'] <-max(c(therm.set[1:301, 10]))
    all.dat[1,'max.actual.Temp.train.B'] <-max(c(therm.set[302:6301, 10]))
    
    
    counter<-2
    for(ii in chambs)
    {
      ss<-counter+9
      ee<-ss+8
      subR<-therm.set[,c(1,ss:ee)]
      colnames(subR)<-colnames(allR)
      all.dat[ii,'chamber'] <- subR$chamber[1]
      all.dat[ii,'incapacitation'] <- slideThermo(xx=subR$pos, tt = subR$time)
      incap <- slideThermo(xx=subR$pos, tt = subR$time)
      all.dat[ii,'incapacitation'] <- incap
      all.dat[ii, "ThermAct"] <- var(subR$pos) 
      if(incap > 629999){
        all.dat[ii,'Rvar'] <- NA 
      }else{
        all.dat[ii,'Rvar'] <- slideRec(which(subR$time==incap),xx=subR$pos, tt = subR$time)
      }
      gg.1 <- substring(gg, 6)
      all.dat[ii,'group'] <- as.numeric(gsub("[^\\d]+", "", gg.1,perl=TRUE))
      all.dat[ii,'file'] <- gg
      all.dat[ii,'target.Temp.test.A'] <-mean(c(subR[1:301, 7]))
      all.dat[ii,'actual.Temp.test.A'] <-mean(c(subR[1:301, 8]))
      all.dat[ii,'target.Temp.test.B'] <-mean(c(subR[1:301, 9]))
      all.dat[ii,'actual.Temp.test.B'] <-mean(c(subR[1:301, 10]))
      all.dat[ii,'target.Temp.train.A'] <-mean(c(subR[302:6301, 7]))
      all.dat[ii,'actual.Temp.train.A'] <-mean(c(subR[302:6301, 8]))
      all.dat[ii,'target.Temp.train.B'] <-mean(c(subR[302:6301, 9]))
      all.dat[ii,'actual.Temp.train.B'] <-mean(c(subR[302:6301, 10]))
      all.dat[ii,'max.actual.Temp.test.A'] <-max(c(subR[1:301, 8]))
      all.dat[ii,'max.actual.Temp.train.A'] <-max(c(subR[302:6301, 8]))
      all.dat[ii,'max.actual.Temp.test.B'] <-max(c(subR[1:301, 10]))
      all.dat[ii,'max.actual.Temp.train.B'] <-max(c(subR[302:6301, 10]))
      
      
      counter<-counter+9
      
    } #ii
    ThermTol <- rbind(ThermTol, all.dat)
    #cat(min(allR$pos),"\t", max(allR$pos),"\n")
  } #gg
}#kk

ThermTol$incapacitation <- ThermTol$incapacitation/1000
ThermTol$incapacitation <- ThermTol$incapacitation - 30

tester <- ThermTol[which(ThermTol$incapacitation==0),]

#ThermTol <- ThermTol[-which(ThermTol$incapacitation==0),]

#save the R object
save(ThermTol, file="../Data/Incapacitation.rda")

