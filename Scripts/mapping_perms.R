library(DSPRqtl)

source("logLikmulti.R")


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

all.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) all.LL[zz,] <- Phen.Ls[[zz]]
colnames(all.LL)<-colnames(pp)

#########
#
#
#
#

##general information of genom scan












#make qtl maps

#set dataset as a date frame

myGenos <- as.data.frame(myGenos$phenotype)


#Create peaks for qtl
myGenospeaks <- DSPRpeaks(myGenos, threshold = 6.8, LODdrop = 2)


save(Learningpeaks,file= "../Data/myGenospeaks.rda")


#Create qtl map

myGenos_qtlmap <- DSPRplot(list(myGenos), threshold=6.8)

myGenos_qtlmap

ggsave(myGenos_qtlmap, file="../Plots/myGenos_qtlmap.pdf", width=10, height=4)

dev.off()

