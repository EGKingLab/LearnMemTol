library(DSPRqtl)
source("../../Functions/mappingfunctions.R")


##Load datasets

load(file="../../LearnMemTolPheno/ProcessedData/L_MDATA.rda")
load(file="../../LearnMemTolPheno/ProcessedData/T_TDATA.rda")

#View dataset
str(L_MDATA)
str(T_TDATA)

#merge by patRIL
All_3_Traits_Data <- merge(L_MDATA, T_TDATA, by='patRIL')


load(file="../ProcessedData/myGenos.rda")

#Look at the possitions. These are in 5.x coordinates
positions <- myGenos$positions

#Put column names of phenotypes here
pp<-myGenos[['phenotype']][,c('LearnPowerTrans','Memory_Mean','Tolsqrtvariable')]

#load previous mapping results
load(file="../ProcessedData/Lodscores_3traits.rda")

tqtl <- which.max(obs.LL[,3])

gg_cor <- myGenos$genolist[[tqtl]]

LL0 <- logLik(lm(pp[,3] ~ gg_cor[,1] + 
                   gg_cor[,2] +
                   gg_cor[,3] +
                   gg_cor[,4] +
                   gg_cor[,5] +
                   gg_cor[,6] +
                   gg_cor[,7]))/log(10)


gg <- myGenos$genolist
Phen.Mods<-lapply(gg, function(mat) logLik(lm(pp[,3] ~ gg_cor[,1] + 
                                                gg_cor[,2] +
                                                gg_cor[,3] +
                                                gg_cor[,4] +
                                                gg_cor[,5] +
                                                gg_cor[,6] +
                                                gg_cor[,7] +
                                                mat[,1] + 
                                                mat[,2] + 
                                                mat[,3] + 
                                                mat[,4] + 
                                                mat[,5] + 
                                                mat[,6] + 
                                                mat[,7]))/log(10))

Phen.Ls <- unlist(lapply(Phen.Mods, function(xx) (xx - LL0)))

plot(Phen.Ls, type='l')
lines(obs.LL[,3], col='blue')


pp<-myGenos[['phenotype']]

pp$subpop <- 'A'
pp$subpop[pp$patRIL>12000]<-'B'

LL0 <- logLik(lm(pp[,"Tolsqrtvariable"] ~ pp[,'subpop']))/log(10)


gg <- myGenos$genolist
Phen.Mods<-lapply(gg, function(mat) logLik(lm(pp[,"Tolsqrtvariable"] ~ pp[,'subpop'] +
                                                mat[,1] + 
                                                mat[,2] + 
                                                mat[,3] + 
                                                mat[,4] + 
                                                mat[,5] + 
                                                mat[,6] + 
                                                mat[,7]))/log(10))

Phen.Ls <- unlist(lapply(Phen.Mods, function(xx) (xx - LL0)))


plot(obs.LL[,3], type='l')
lines(Phen.Ls, col='blue')

# Cox hazard??
