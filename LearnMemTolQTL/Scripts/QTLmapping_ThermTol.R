library(DSPRqtl)
library(DSPRqtlDataA)


#Load in processed data, named IncapData or IncapDatatransformed (if you transformed the data)
load(file="../LearnMemTolPheno/ProcessedData/T_TDATA.rda")

cat("data in","\n")

str(T_TDATA)

#Turn into data frame
T_TDATA <- as.data.frame(T_TDATA)


#genome scans 
scanresults_ThermTol <- DSPRscan (Tolsqrtvariable ~ 1, 
                                    design = 'inbredA', 
                                    phenotype.dat = T_TDATA, 
                                    id.col = "patRIL")

cat("Incap scan done","\n")

#to view file
scanresults_ThermTol

#to save genome scans
save(scanresults_ThermTol, file="../ProcessedData/scanresults_TherrmTol.rda")




#running permutation test for significant thereshold
permresults_ThermTol <- DSPRperm(Tolsqrtvariable ~ 1,
                                    design = "inbredA",
                                    phenotype.dat = T_TDATA,
                                    id.col='patRIL',
                                    batch = 1000, niter = 1000,
                                    alpha = 0.05)
cat("incap perms done","\n")


permresults_ThermTol

save(permresults_TheermTol, file="../ProcessedData/permresults_ThermTol.rda")



#
#
#


load(file="../ProcessedData/scanresults_ThermTol.rda")

load(file="../ProcessedData/permresults_ThermTol.rda")

str(permresults_ThermTol)

scanresults_ThermTol$phenotype <- as.data.frame(scanresults_ThermTol$phenotype)


#to create peaks for qtl
ThermTolpeaks <- DSPRpeaks(scanresults_ThermTol, threshold = 6.7, LODdrop = 2)

save(ThermTolpeaks,file= "../ProcessedData/ThermTolpeaks.rda")

#to view confidence intervals of significant peaks
ThermTolpeaks[[1]]
ThermTolpeaks[[11]]


#to zoom in on the qtl map
scanresults_ThermTol$LODscores[14500000:14700000,]



plot(scanresults_ThermTol$LODscores$LOD,type='l')


#to view data
glimpse(ThermTolpeaks)


#general information on qtl
DSPRpeaks(ThermTolpeaks, method, threshold, LODdrop, BCIprob)


#to see the tallest peak
max(scanresults_ThermTol$LOD)


#shows the location of the tallest peaks in genome
which.max(scanresults_ThermTol$LOD)


#to look at data for significant peaks
str(ThermTolpeaks)
ThermTolpeaks[[3]]

#to make a QTL map
ThermTolqtl <- DSPRplot(list(scanresults_ThermTol), threshold=6.7)

#to view qtl map
ThermTolqtl

#to save qtl map
ggsave(ThermTolqtl, file="../Plots/ThermTolqtlmap.pdf", width = 10, height = 4)
       