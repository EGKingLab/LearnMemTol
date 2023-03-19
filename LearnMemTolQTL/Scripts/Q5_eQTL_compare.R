library(DSPRqtlDataA)
library(DSPRqtl)
library(tidyverse)
library(DESeq2)

fbgns <-read_csv(file = "../../LearnMemTolPheno/ProcessedData/Q5_geneids.csv")


load(file="/home/kingeg/Archived/wfitchNOV2015/Fheads_eQTLs/eQTLs/AbyB/Analysis/NewRMA/QTLmapping/Final_EXPR_MAT.rda")

load(file="../../LearnMemTolPheno/ProcessedData/T_TDATA.rda")
str(T_TDATA)
tww<- which(T_TDATA$patRIL %in% rownames(qn.mat))
twwe<- which(rownames(qn.mat) %in% T_TDATA$patRIL)
TT.sub <- T_TDATA[tww,]
TTqn.mat <- qn.mat[twwe,]

all.equal(as.numeric(row.names(TTqn.mat)), TT.sub$patRIL)

ccs.t <- cor(TT.sub$Incapacitation_Mean, TTqn.mat)

cors.tt <- data.frame('ID' = colnames(TTqn.mat), 'rThermTol' = as.numeric(ccs.t), stringsAsFactors=FALSE)

#load eQTL data from king et  al 2014
cgn<-read.table(file="/home/kingeg/Archived/wfitchNOV2015/Fheads_eQTLs/eQTLs/AbyB/Analysis/NewRMA/QTLmapping/fbgn_fbtr_fbpp_fb_2012_02.tsv",header=FALSE,stringsAsFactors=FALSE,sep="\t",fill=TRUE)

load(file='/home/kingeg/Archived/wfitchNOV2015/Fheads_eQTLs/eQTLs/probes/CDS/completeCDS.rda')

ss<-merge(cgn,bl, by.x="V2",by.y="FBtr")

ww.tt <- which(ss$CGid %in% colnames(TTqn.mat))
ww.tg <- which(ss$CGgn %in% colnames(TTqn.mat))
decode.tn <- ss[,c('V1','CGid','CGgn')]
decode.tn$ID <- NA
decode.tn$ID[ww.tt]<- ss$CGid[ww.tt]
decode.tn$ID[ww.tg]<- ss$CGgn[ww.tg]
decode.tn <- decode.tn[,c('V1','ID')]
decode.tn <- unique(decode.tn, MARGIN=1)
colnames(decode.tn) <- c("FBgn","ID")

cors.tt <- merge(cors.tt, decode.tn,by='ID',sort=FALSE, all.x=TRUE)

cors.tnas <- cors.tt[is.na(cors.tt$FBgn),]
cors.tt <- cors.tt[is.na(cors.tt$FBgn)==FALSE,]

cors.tnas[which.max(abs(cors.tnas$rThermTol)),]
cors.tnas[which(abs(cors.tnas$rThermTol) >0.1),]


q5c <- subset(cors.tt, FBgn %in% fbgns$FBgn)

data("positionlist_wgenetic")


load(file ="../ProcessedData/Peaks_wCIs.rda")

rrs <- read.table(file="../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/RILs_in_pools.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)

rrcode <- c('LH','LL',
            'MH','ML',
            'TH','TL')

cc <- 1

for(kk in 1:3)
{
  
  list.peak <- ci.peak[[kk]]
  list.peak[,43:58]<-NA
  colnames(list.peak)[43:58] <- c("L1","L2","L3","L4","L5","L6","L7","L8",
                                  "H1","H2","H3","H4","H5","H6","H7","H8")
  for(jj in 1:nrow(list.peak))
  {
    ppeak <- list.peak[jj,]
    
    objname<-paste("A_",ppeak$chr,"_",ppeak$Ppos,sep="")
    data(list=objname)
    geno <- get(objname)
    geno <- geno[,c("ril", "AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8")]
    lowset <- geno %>% filter(ril %in% rrs[rrs$Pool==rrcode[(cc+1)],'patRIL'])
    highset <- geno %>% filter(ril %in% rrs[rrs$Pool==rrcode[cc],'patRIL'])
    ppeak[,43:50]<-round(colSums(lowset[,2:9]))
    ppeak[,51:58]<-round(colSums(highset[,2:9]))
    list.peak[jj,]<-ppeak
  }
  cc<-cc+1
  ci.peak[[kk]]<-list.peak
}

load(file = "../../LearnMemTolRnaseq/DEseq_LearnMemTol/Data/FoldChangeAll.Rda")
FC_All$FBgn <- rownames(FC_All)

gene_map_table <- read_lines(file = "../ProcessedData/gene_map_table_fb_2018_04.tsv",
                             skip = 6) %>% 
  str_split(pattern = "\t", simplify = TRUE) %>% 
  as_tibble() %>% 
  filter(V1 == "Dmel")

colnames(gene_map_table) <- c('spp','gname','FBgn','v3','cyt','pos')


gene_map_table$chr <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,1]
temp.s1 <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,2] %>% 
  str_split(fixed(".."),simplify=TRUE)
gene_map_table$startp <- as.numeric(temp.s1[,1])

gene_map_table$stopp <- as.numeric(str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1])

TT_gene_list_merged <- merge(FC_All, gene_map_table, by="FBgn")


ThermTol.list <-vector(mode='list', length=nrow(ci.peak$Tolsqrtvariable))

for(rowN in 1:nrow(ci.peak$Tolsqrtvariable)) 
{
  
  foc.peak<-ci.peak[['Tolsqrtvariable']][rowN,]
  
  #Q1 spans centromere so needs a different logical
  if(rowN ==1)
  {
    t_all <- c(which(TT_gene_list_merged$chr==foc.peak$lpchr & 
                       TT_gene_list_merged$startp >= foc.peak$lpR6), 
               which(TT_gene_list_merged$chr==foc.peak$upchr & 
                       TT_gene_list_merged$startp <= foc.peak$upR6)
    )
  }else{
    t_all <- (which(TT_gene_list_merged$chr==foc.peak$chr & 
                      ((TT_gene_list_merged$startp <= foc.peak$upR6 & TT_gene_list_merged$stopp >= foc.peak$lpR6))))  
    
  }
  TT_genes_under_peak<-TT_gene_list_merged[t_all,]
  
  ThermTol.list[[rowN]] <-TT_genes_under_peak
  
}
#files from old proj
load(file="/home/kingeg/Archived/wfitchNOV2015/Fheads_eQTLs/eQTLs/AbyB/Analysis/NewRMA/QTLmapping/goodpeaksinfo_matrix.rda")
load(file="/home/kingeg/Archived/wfitchNOV2015/Fheads_eQTLs/eQTLs/AbyB/Analysis/NewRMA/QTLmapping/means_wcorrected.rda")

effects <- effects[pinfo$cis==TRUE,]

pinfo <- pinfo[pinfo$cis==TRUE,]

effects <- merge(effects, decode.tn, by='ID',sort=FALSE, all.x=TRUE)

for(jj in 1:length(ThermTol.list))
{
  
  ll.foc <- ThermTol.list[[jj]]
  pp.foc <- ci.peak[[3]][jj,]
  
  tw <- (which(pinfo$chr==pp.foc$chr & 
                 ((pinfo$peaklpB <= pp.foc$up & pinfo$peakupB >= pp.foc$lp))))  
  
  ee <- effects[effects$ID %in% pinfo$ID[tw],]
  
  f.high <- which(apply(rbind(as.numeric(pp.foc[,35:42]), as.matrix(ee[,34:41])),2, min) > 4)
  
  qtl.e <- pp.foc[,19:26]
  eqtl.e <- ee[,2:9]
  
  qtl.e <- qtl.e[,f.high]
  eqtl.e <- eqtl.e[,f.high]
  
  cors.effs <- data.frame('FBgn'= ee$FBgn, 'corEff'= apply(eqtl.e, 1, function(x) cor(x,as.numeric(qtl.e))), stringsAsFactors = FALSE)
  
  ll.all <- merge(ll.foc, cors.effs, by='FBgn', all=TRUE)
  
  ll.all <- merge(ll.all, cors.tt, by='FBgn', all.x=TRUE)
  ThermTol.list[[jj]] <- ll.all
}

