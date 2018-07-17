library(DSPRqtl)

source("logLikmulti.R")

#phenotype.dat should be all three phenotypes with patRIL 
#use merge and all = TRUE, missing vals will be NA

myGenos <- DSPRgenos(design='inbredA', phenotype.dat, id.col='patRIL')

#probably a good idea to save for later so you can just load instead of running DSPRgenos again
save(myGenos, file="../Data/myGenos.rda")

#these are in 5.x coordinates
positions <- myGenos$positions

#put column names of phenotypes here
pp<-myGenos[['phenotype']][,c('Learn','Mem','Therm')]

Null.Mods<-logLik.multi(lm(as.matrix(pp)~1))

Null.Ls<-unlist(Null.Mods)/log(10)

gg <- myGenos$genolist
Phen.Mods<-lapply(gg, function(mat) logLik.multi(lm(as.matrix(pp) ~ mat[,1] + mat[,2]+ mat[,3]+ mat[,4]+ mat[,5]+ mat[,6]+ mat[,7])))

Phen.Ls<-lapply(Phen.Mods, function(xx) (xx/log(10)) - Null.Ls)

all.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) all.LL[zz,] <- Phen.Ls[[zz]]
colnames(all.LL)<-colnames(pp)




