
library(DSPRqtl)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

load(file="../../LearnMemTolPheno/ProcessedData/T_TDATA.rda")

rr <- readRDS("../ProcessedData/Q5_RIL_SNP_probs.rds")

check_s <- rowSums(rr)
max(check_s)

ww <- which(check_s < (741-4) & check_s > 4)

rr <- rr[ww,]
rr_pos <- as.numeric(str_split(rownames(rr), "_", simplify=TRUE)[,2])

all.equal(T_TDATA$patRIL, as.numeric(colnames(rr)))

output <- data.frame("pos" = rep(NA,nrow(rr)),
                     "pval" = rep(NA,nrow(rr)),
                     "effect" = rep(NA,nrow(rr)))

for(ii in 1:nrow(rr))
{
  pp <- data.frame("pheno" = T_TDATA$Tolsqrtvariable,
                   "geno" = as.numeric(rr[ii,]))
  mod <- lm(pheno ~ geno, data = pp)
  mod_sum <- summary(mod)
  output[ii, "pos"] <- rr_pos[ii]
  output[ii,"pval"] <- -pf(mod_sum$fstatistic[1], mod_sum$fstatistic[2], mod_sum$fstatistic[3], lower.tail = FALSE, log.p=TRUE)/log(10)
  output[ii,"effect"] <- mod_sum$coefficients[2,1]
}

ggplot(output, aes(pos/1e6, pval)) +
  geom_point(alpha = 1/2) 
  


