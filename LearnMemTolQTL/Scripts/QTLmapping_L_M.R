library(DSPRqtl)
library(DSPRqtlDataA)


#Load in processed data, named DSPRDATA
load(file="../../LearnMemTolPheno/ProcessedData/L_MDATA.rda")

cat("data in","\n")

#Turn into data frame
L_MDATA <- as.data.frame(L_MDATA)


#genome scans 

#stop!!!!
#use code below if data is norally distributed, but if transformed use the second one.
#be sure to change code to fit the type of transforation done.

#use this code if the data is normally distributed 
scanresults_Learning <- DSPRscan (Learning_Mean ~ 1, 
                                    design = 'inbredA', 
                                    phenotype.dat = L_MDATA, 
                                    id.col = "patRIL")

cat("learning scan done","\n")

#use this code if the data is transformed 
scanresults_Learning <- DSPRscan (LearnPowerTrans ~ 1, 
                                  design = 'inbredA', 
                                  phenotype.dat = L_MDATA, 
                                  id.col = "patRIL")

cat("learning scan done","\n")

save(scanresults_Learning, file="../ProcessedData/scanresults_Learning.rda")


scanresults_Memory <- DSPRscan (Memory_Mean ~ 1, 
                              design = 'inbredA', 
                              phenotype.dat = L_MDATA, 
                              id.col = "patRIL")
cat("memory scan done","\n")


save(scanresults_Memory, file="../ProcessedData/scanresults_Memory.rda")

#to view scan results 
scanresults_Learning
scanresults_Memory


#running permutation test for significant thereshold
permresults_Learning <- DSPRperm(Learning_Mean ~ 1,
                                    design = "inbredA",
                                    phenotype.dat = L_MDATA,
                                    id.col='patRIL',
                                    batch = 1000, niter = 1000,
                                    alpha = 0.05)
cat("learning perms done","\n")


#run this permutation test if data was transformed.

permresults_Learning <- DSPRperm(LearnPowerTrans ~ 1,
                                 design = "inbredA",
                                 phenotype.dat = L_MDATA,
                                 id.col='patRIL',
                                 batch = 1000, niter = 1000,
                                 alpha = 0.05)
cat("learning perms done","\n")

#to view results 
permresults_Learning

#to save results
save(permresults_Learning, file="../ProcessedData/permresults_Learning.rda")



permresults_Memory <- DSPRperm(Memory_Mean ~ 1,
                                  design = "inbredA",
                                 phenotype.dat = L_MDATA,
                                 id.col='patRIL',
                                 batch = 1000, niter = 1000,
                                 alpha = 0.05)
cat("memory perms done","\n")

save(permresults_Memory, file="../ProcessedData/permresults_Memory.rda")



#
#
#

load(file="../ProcessedData/scanresults_Learning.rda")
load(file="../ProcessedData/scanresults_Memory.rda")

scanresults_Learning$phenotype <- as.data.frame(scanresults_Learning$phenotype)
scanresults_Memory$phenotype <- as.data.frame(scanresults_Memory$phenotype)

#to create peaks for qtl
Learningpeaks <- DSPRpeaks(scanresults_Learning, threshold = 6.73, LODdrop = 2)
Memorypeaks <- DSPRpeaks(scanresults_Memory, threshold = 6.65, LODdrop = 2)

save(Learningpeaks,file= "../ProcessedData/Learningpeaks.rda")
save(Memorypeaks,file= "../ProcessedData/Memorypeaks.rda")



#to view confidence intervals of significant peaks
Learningpeaks[[1]]
Learningpeaks[[6]]
Memorypeaks[[3]]
Memorypeaks[[7]]


#to view specific location on chromosome
scanresults_Learning$LODscores[10800:10810,]

#to zoom in on the qtl map
plot(scanresults_Learning$LODscores$LOD,type='l')



glimpse(Memorypeaks)

#general information on qtl

DSPRpeaks(Learningpeaks, method, threshold, LODdrop, BCIprob)



#to see the tallest peak
max(scanresults_Learning$LOD)
max(scanresults_Memory$LOD)


#shows the location of the peaks in genome
which.max(scanresults_Learning$LOD)


#to find the main peaks
str(Memorypeaks)
Memorypeaks[[1]]

#to make a QTL map for learning 
learning_qtlmap <- DSPRplot(list(scanresults_Learning), threshold=6.73)

learning_qtlmap

ggsave(learning_qtlmap, file="../Plots/Learning_qtl.pdf", width=10, height=4)

dev.off()



#to make a QTL map for memory
memory_qtlmap <- DSPRplot(list(scanresults_Memory), threshold=6.65)

memory_qtlmap

ggsave(memory_qtlmap, file="../Plots/Memory_qtl.pdf", width=10, height=4)

dev.off()
