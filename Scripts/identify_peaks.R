library(DSPRqtl)
data("positionlist_wgenetic")
load(file = "/home/kingeg/Projects/DSPRgeneral/Convert5_6/DSPR_5_6_parsed.rda")

source("mappingfunctions.R")

load(file="../Data/sig_ths.rda")#fdr.out
load(file="../Data/Lodscores_3traits.rda")#obs.LL

#th.l <- fdr.out$fwer
th.l <- fdr.out$fdr$th[min(which(fdr.out$fdr$fdr<=0.05))]

#find peaks
peak.i <- apply(obs.LL[,1:3],2,function(x) peakInfo(x,cM=poslist[,c('chr','Gpos')] ,th=th.l, tol.dist=3))



ci.peak <- vector(mode="list", length=3)
names(ci.peak)<-names(peak.i)

for(kk in 1:3){
qq <- poslist[,1:3]
qq$LOD <- obs.LL[,kk]
pp.set <- peak.i[[kk]]
pp.set$up <- NA
pp.set$lp <- NA
pp.set$ug <- NA
pp.set$lg <- NA
pp.set$ulod <- NA
pp.set$llod <- NA


for(jj in 1:nrow(pp.set))
{
ff<-findCI(pp.set$chr[jj],pp.set$Ppos[jj], qtldat=qq, method='BCI')
pp.set$up[jj] <- ff$Ppos[2]
pp.set$lp[jj] <- ff$Ppos[1]
pp.set$ug[jj] <- ff$Gpos[2]
pp.set$lg[jj] <- ff$Gpos[1]
pp.set$ulod[jj] <- ff$LOD[2]
pp.set$llod[jj] <- ff$LOD[1] 
} 
pp.set[,c('chr','Ppos')] <- merge(pp.set[,c('chr','Ppos')], coord.table, by.x=c('chr','Ppos'), by.y=c('R5chr','R5pos'))[,c('R6chr','R6pos')]
pp.set[,c('up')] <- merge(pp.set[,c('chr','up')], coord.table, by.x=c('chr','up'), by.y=c('R5chr','R5pos'))[,c('R6pos')]
pp.set[,c('lp')] <- merge(pp.set[,c('chr','lp')], coord.table, by.x=c('chr','lp'), by.y=c('R5chr','R5pos'))[,c('R6pos')]
pp.set<-pp.set[,order(c('chr','Ppos'))]

ci.peak[[kk]]<-pp.set
}

save(ci.peak,file="../Data/Peaks_wCIs.rda")


#sort by position



