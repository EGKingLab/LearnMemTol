options(scipen=999)


#load genome scan data
#load(file="../ProcessedData/Tolpeaks.rda")
load(file="../ProcessedData/Memorypeaks.rda")
load(file="../ProcessedData/Learningpeaks.rda")


#read in signigicant gene list
L.DEG <- read.csv("../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/LearnDEseqSVA_resOrdered_padj.csv",header=TRUE,stringsAsFactors = FALSE)
M.DEG <- read.csv("../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/MemDEseqSVA_resOrdered_padj.csv",header=TRUE,stringsAsFactors = FALSE)

#I.DEG <- read.csv("../RawData/Learn_Mem_Tol_rna-seq/Significant_Tolgenes_Results.csv",header=TRUE,stringsAsFactors = FALSE)


load(file = "/home/kingeg/Projects/DSPRgeneral/Convert5_6/DSPR_5_6_parsed.rda")

str(coord.table)

Conv_5_6<-function(ccc,ppp,coord.table)

  {
  return(coord.table[which(coord.table$R5chr==ccc & coord.table$R5pos==ppp),c('R6chr','R6pos')])
}


Conv_5_6('X',1200000,coord.table)
str(Learningpeaks[[1]])

Learn_p <- data.frame('CHROM'=character(length=length(Learningpeaks)), 
                      'PeakP'=numeric(length=length(Learningpeaks)), 
                      'Gpos'=numeric(length=length(Learningpeaks)), 
                      'LOD'=numeric(length=length(Learningpeaks)),
                      'UL'=numeric(length=length(Learningpeaks)),
                      'LL'=numeric(length=length(Learningpeaks)),
                      'UB'=numeric(length=length(Learningpeaks)),
                      'LB'=numeric(length=length(Learningpeaks)),stringsAsFactors = FALSE)


for(kk in 1:length(Learningpeaks))
{
  cc <- Learningpeaks[[kk]]$peak$chr
  Learn_p[kk,'CHROM']<- cc
  Learn_p[kk,'PeakP']<- as.numeric(Conv_5_6(cc, Learningpeaks[[kk]]$peak$Ppos, coord.table)$R6pos)
  Learn_p[kk,'Gpos']<- Learningpeaks[[kk]]$peak$Gpos
  Learn_p[kk,'LOD']<- Learningpeaks[[kk]]$peak$LOD
  Learn_p[kk,'UL']<- as.numeric(Conv_5_6(Learningpeaks[[kk]]$CI$LODdrop[2,1], Learningpeaks[[kk]]$CI$LODdrop[2,2], coord.table)$R6pos)
  Learn_p[kk,'LL']<- as.numeric(Conv_5_6(Learningpeaks[[kk]]$CI$LODdrop[1,1], Learningpeaks[[kk]]$CI$LODdrop[1,2], coord.table)$R6pos)
  Learn_p[kk,'UB']<- as.numeric(Conv_5_6(Learningpeaks[[kk]]$CI$LODdrop[2,1], Learningpeaks[[kk]]$CI$BCI[2,2], coord.table)$R6pos)
  Learn_p[kk,'LB']<- as.numeric(Conv_5_6(Learningpeaks[[kk]]$CI$LODdrop[1,1], Learningpeaks[[kk]]$CI$BCI[1,2], coord.table)$R6pos)
 cat(kk,"\n") 
  
}




#converting 5-6 memory

Mem_p <- data.frame('CHROM'=character(length=length(Memorypeaks)), 
                      'PeakP'=numeric(length=length(Memorypeaks)), 
                      'Gpos'=numeric(length=length(Memorypeaks)), 
                      'LOD'=numeric(length=length(Memorypeaks)),
                      'UL'=numeric(length=length(Memorypeaks)),
                      'LL'=numeric(length=length(Memorypeaks)),
                      'UB'=numeric(length=length(Memorypeaks)),
                      'LB'=numeric(length=length(Memorypeaks)),stringsAsFactors = FALSE)


for(kk in 1:length(Memorypeaks))
{
  cc <- Memorypeaks[[kk]]$peak$chr
  Mem_p[kk,'CHROM']<- cc
  Mem_p[kk,'PeakP']<- as.numeric(Conv_5_6(cc, Memorypeaks[[kk]]$peak$Ppos, coord.table)$R6pos)
  Mem_p[kk,'Gpos']<- Memorypeaks[[kk]]$peak$Gpos
  Mem_p[kk,'LOD']<- Memorypeaks[[kk]]$peak$LOD
  Mem_p[kk,'UL']<- as.numeric(Conv_5_6(Memorypeaks[[kk]]$CI$LODdrop[2,1], Memorypeaks[[kk]]$CI$LODdrop[2,2], coord.table)$R6pos)
  Mem_p[kk,'LL']<- as.numeric(Conv_5_6(Memorypeaks[[kk]]$CI$LODdrop[1,1], Memorypeaks[[kk]]$CI$LODdrop[1,2], coord.table)$R6pos)
  Mem_p[kk,'UB']<- as.numeric(Conv_5_6(Memorypeaks[[kk]]$CI$LODdrop[2,1], Memorypeaks[[kk]]$CI$BCI[2,2], coord.table)$R6pos)
  Mem_p[kk,'LB']<- as.numeric(Conv_5_6(Memorypeaks[[kk]]$CI$LODdrop[1,1], Memorypeaks[[kk]]$CI$BCI[1,2], coord.table)$R6pos)
  cat(kk,"\n") 
  
}


#Thermal Tolerance
Incap_p <- data.frame('CHROM'=character(length=length(Incappeaks)), 
                      'PeakP'=numeric(length=length(Incappeaks)), 
                      'Gpos'=numeric(length=length(Incappeaks)), 
                      'LOD'=numeric(length=length(Incappeaks)),
                      'UL'=numeric(length=length(Incappeaks)),
                      'LL'=numeric(length=length(Incappeaks)),
                      'UB'=numeric(length=length(Incappeaks)),
                      'LB'=numeric(length=length(Incappeaks)),stringsAsFactors = FALSE)

#hack
coord.table <- rbind(coord.table,c('2R','21410000','NA','NA'))

for(kk in 1:length(Incappeaks))
{
  cc <- Incappeaks[[kk]]$peak$chr
  Incap_p[kk,'CHROM']<- cc
  Incap_p[kk,'PeakP']<- as.numeric(Conv_5_6(cc, Incappeaks[[kk]]$peak$Ppos, coord.table)$R6pos)
  Incap_p[kk,'Gpos']<- Incappeaks[[kk]]$peak$Gpos
  Incap_p[kk,'LOD']<- Incappeaks[[kk]]$peak$LOD
  Incap_p[kk,'UL']<- as.numeric(Conv_5_6(Incappeaks[[kk]]$CI$LODdrop[2,1], Incappeaks[[kk]]$CI$LODdrop[2,2], coord.table)$R6pos)
  Incap_p[kk,'LL']<- as.numeric(Conv_5_6(Incappeaks[[kk]]$CI$LODdrop[1,1], Incappeaks[[kk]]$CI$LODdrop[1,2], coord.table)$R6pos)
  Incap_p[kk,'UB']<- as.numeric(Conv_5_6(Incappeaks[[kk]]$CI$LODdrop[2,1], Incappeaks[[kk]]$CI$BCI[2,2], coord.table)$R6pos)
  Incap_p[kk,'LB']<- as.numeric(Conv_5_6(Incappeaks[[kk]]$CI$LODdrop[1,1], Incappeaks[[kk]]$CI$BCI[1,2], coord.table)$R6pos)
  cat(kk,"\n") 
  
}

#hist(L.DEG[,'LOCATION_MAX']-L.DEG[,'LOCATION_MIN'])
I.dd.L<-numeric(0)
I.dd.B<-numeric(0)
for(ii in 1:nrow(L.DEG))
{
  lw <- (which(Learn_p$CHROM==L.DEG$LOCATION_ARM[ii] & 
          ((Learn_p$UL>L.DEG$LOCATION_MIN[ii] & Learn_p$LL<L.DEG$LOCATION_MIN[ii])  
           |(Learn_p$LL<L.DEG$LOCATION_MIN[ii] & Learn_p$UL>L.DEG$LOCATION_MIN[ii]))))
  
  bw <-which(Learn_p$CHROM==L.DEG$LOCATION_ARM[ii] & 
              ((Learn_p$UB>L.DEG$LOCATION_MIN[ii] & Learn_p$LB<L.DEG$LOCATION_MIN[ii])  
               |(Learn_p$LB<L.DEG$LOCATION_MIN[ii] & Learn_p$UB>L.DEG$LOCATION_MIN[ii])))
  
  #cat(ii, lw, bw,"\n")
  if(length(lw)>0){
    I.dd.L <-c(I.dd.L,ii)
  }
  if(length(bw)>0){
    I.dd.B <-c(I.dd.B,ii)
  }
  }
L.DEG[I.dd.L,'SYMBOL']
L.DEG[I.dd.B,'SYMBOL']

I.dd.L<-numeric(0)
I.dd.B<-numeric(0)
for(ii in 1:nrow(M.DEG))
{
  lw <- (which(Mem_p$CHROM==M.DEG$LOCATION_ARM[ii] & 
                 ((Mem_p$UL>M.DEG$LOCATION_MIN[ii] & Mem_p$LL<M.DEG$LOCATION_MIN[ii])  
                  |(Mem_p$LL<M.DEG$LOCATION_MIN[ii] & Mem_p$UL>M.DEG$LOCATION_MIN[ii]))))
  
  bw <-which(Mem_p$CHROM==M.DEG$LOCATION_ARM[ii] & 
               ((Mem_p$UB>M.DEG$LOCATION_MIN[ii] & Mem_p$LB<M.DEG$LOCATION_MIN[ii])  
                |(Mem_p$LB<M.DEG$LOCATION_MIN[ii] & Mem_p$UB>M.DEG$LOCATION_MIN[ii])))
  
  #cat(ii, lw, bw,"\n")
  if(length(lw)>0){
    I.dd.L <-c(I.dd.L,ii)
  }
  if(length(bw)>0){
    I.dd.B <-c(I.dd.B,ii)
  }
}
M.DEG[I.dd.L,'SYMBOL']
M.DEG[I.dd.B,'SYMBOL']


I.dd.L<-numeric(0)
I.dd.B<-numeric(0)
for(ii in 1:nrow(I.DEG))
{
  lw <- (which(Incap_p$CHROM==I.DEG$LOCATION_ARM[ii] & 
                 ((Incap_p$UL>I.DEG$LOCATION_MIN[ii] & Incap_p$LL<I.DEG$LOCATION_MIN[ii])  
                  |(Incap_p$LL<I.DEG$LOCATION_MIN[ii] & Incap_p$UL>I.DEG$LOCATION_MIN[ii]))))
  
  bw <-which(Incap_p$CHROM==I.DEG$LOCATION_ARM[ii] & 
               ((Incap_p$UB>I.DEG$LOCATION_MIN[ii] & Incap_p$LB<I.DEG$LOCATION_MIN[ii])  
                |(Incap_p$LB<I.DEG$LOCATION_MIN[ii] & Incap_p$UB>I.DEG$LOCATION_MIN[ii])))
  
  #cat(ii, lw, bw,"\n")
  if(length(lw)>0){
    I.dd.L <-c(I.dd.L,ii)
  }
  if(length(bw)>0){
    I.dd.B <-c(I.dd.B,ii)
  }
}

I.DEG[I.dd.L,'SYMBOL']
I.DEG[I.dd.B,'SYMBOL']


