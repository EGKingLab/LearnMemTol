library(DSPRqtlDataA)

rawdLM <- read.table("../../LearnMemTolPheno/ProcessedData/LearnMem_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)
load(file="../../LearnMemTolPheno/ProcessedData/L_MDATA.rda")

load(file ="../ProcessedData/Peaks_wCIs.rda")

#Q4,Q11,Q14

data(A_2L_3190000)
hq4 <- apply(A_2L_3190000[,2:9],1,which.max)
hq4[apply(A_2L_3190000[,2:9],1,max)<0.95] <- NA

data(A_3L_8560000)
hq11 <- apply(A_3L_8560000[,2:9],1,which.max)
hq11[apply(A_3L_8560000[,2:9],1,max)<0.95] <- NA

data(A_3R_18280000)
hq14 <- apply(A_3R_18280000[,2:9],1,which.max)
hq14[apply(A_3R_18280000[,2:9],1,max)<0.95] <- NA

hinfo <- data.frame("patRIL"=as.integer(A_2L_3190000[,1]),
                    "HQ4" =hq4,
                    "HQ11" =hq11,
                    "HQ14" =hq14)

dd <- merge(L_MDATA,hinfo, by="patRIL")

z <- prcomp(L_MDATA[,c(2,4)], center = TRUE, scale. = TRUE)
dd$score <- z$x[,1]

dd <- dd[order(dd$score),]


low <- dd[1:100,]
low <- low[-which(is.na(low$HQ4)|is.na(low$HQ11)|is.na(low$HQ14)),]
low <- low[-which(low$N < 45),]
table(low$HQ4[1:30])
table(low$HQ11[1:30])
table(low$HQ14[1:30])

pp <- low[sort(unique(c(which(low$HQ14==5), which(low$HQ11[1:30] != 8), which(low$HQ4[1:30] %in% c(6,8))))),]

pp <- pp[pp$HQ11 != 8,]
pp<- pp[pp$HQ4 %in% c(3,4,6,7,8),]
pp <- rbind(pp[pp$HQ14 == 5,], pp[pp$HQ14 != 5,])

table(pp$HQ4[1:20])
table(pp$HQ11[1:20])
table(pp$HQ14[1:20])


high <- dd[642:741,]
high <- high[order(-high$score),]
high <- high[-which(is.na(high$HQ4)|is.na(high$HQ11)|is.na(high$HQ14)),]
high <- high[-which(high$N < 45),]
table(high$HQ4[1:30])
table(high$HQ11[1:30])
table(high$HQ14[1:30])

pph <- high[sort(unique(c(which(high$HQ4==5), which(high$HQ11 == 8),seq(1:10)))),]
table(pph$HQ4[1:20])
table(pph$HQ11[1:20])
table(pph$HQ14[1:20])

write.table(sort(c(pp$patRIL,pph$patRIL)), file="../ProcessedData/RTpcr.txt",sep="\t", row.names=FALSE, col.names = FALSE)
