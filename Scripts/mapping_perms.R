library(DSPRqtl)
source("mappingfunctions.R")
set.seed(87523771)



##Load datasets
load(file="../Data/L_MDATA.rda")
load(file="../Data/T_TDATA.rda")

#View dataset
str(L_MDATA)
str(T_TDATA)

#merge by patRIL
All_3_Traits_Data <- merge(L_MDATA, T_TDATA, by='patRIL')


#Genome scan of all three traits concurrently
myGenos <- DSPRgenos(design='inbredA', All_3_Traits_Data, id.col='patRIL')

#Save for later, so you can just load instead of running DSPRgenos again
save(myGenos, file="../Data/myGenos.rda")

#Look at the possitions. These are in 5.x coordinates
positions <- myGenos$positions

#Put column names of phenotypes here
pp<-myGenos[['phenotype']][,c('LearnPowerTrans','Memory_Mean','Tolsqrtvariable')]


Null.Mods<-logLik.multi(lm(as.matrix(pp)~1))

Null.Ls<-unlist(Null.Mods)/log(10)

gg <- myGenos$genolist
Phen.Mods<-lapply(gg, function(mat) logLik.multi(lm(as.matrix(pp) ~ mat[,1] + mat[,2]+ mat[,3]+ mat[,4]+ mat[,5]+ mat[,6]+ mat[,7])))

Phen.Ls<-lapply(Phen.Mods, function(xx) (xx/log(10)) - Null.Ls)

obs.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) obs.LL[zz,] <- Phen.Ls[[zz]]
colnames(obs.LL)<-colnames(pp)

#############################################################
######### PERMUTATIONS ########################
#############################################################
pp<-as.matrix(pp)
#RANDOMIZE PHENOS 1000 TIMES
pp.all<-matrix(NA, nrow(pp), 3000)
cc<-1
for(i in 1:1000)
{
  pp.all[,cc:(cc+2)]<-pp[sample(seq(1,nrow(pp))),]
  cc<-cc+3
}


Null.Mods<-logLik.multi(lm(pp.all~1))

Null.Ls<-unlist(Null.Mods)/log(10)

Phen.Mods<-lapply(gg, function(mat) logLik.multi(lm(pp.all ~ mat[,1] + mat[,2]+ mat[,3]+ mat[,4]+ mat[,5]+ mat[,6]+ mat[,7])))

Phen.Ls<-lapply(Phen.Mods, function(xx) (xx/log(10)) - Null.Ls)

all.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) all.LL[zz,] <- Phen.Ls[[zz]]

apply(all.LL, 2,max)
colnames(all.LL)<-colnames(pp)

save(all.LL, file="Data/Perm_LODS_FDR.rda")

#### Get FDR for different thresholds

th.set<-seq(5,9, by=0.25)

N.positives<-numeric(length(th.set))

for(tt in 1:length(th.set))
{
  N.positives[tt]<-sum(apply(obs.LL,2,function(x) pfind(x,cM=positions[,c('chr','Gpos')] ,th=th.set[tt], tol.dist=5)))
}


N.pos.mat<-matrix(NA,1000,length(th.set))
Max.L<-numeric(1000)

cc<-1
for(ii in 1:1000)
{
  for(tt in 1:length(th.set))
  {
    N.pos.mat[ii,tt]<-sum(apply(all.LL[,cc:(cc+2)],2,function(x) pfind(x,cM=positions[,c('chr','Gpos')] ,th=th.set[tt], tol.dist=5)))
  }
  Max.L[ii]<-max(all.LL[,cc:(cc+2)])
  cc<-cc+3
  
}

fp<-colMeans(N.pos.mat)
fp/N.positives

quantile(Max.L,0.95)

rbind(th.set, fp, N.positives)

tt<-seq(1, 3000,by=3)
tt.m<-apply(all.LL[,tt], 2,max)
ppp<-getP(tt.m,80,7,72)
quantile(ppp, 0.95)
