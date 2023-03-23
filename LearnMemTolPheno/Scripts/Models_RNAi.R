library(tidyverse)
library(multcomp)


TT <- read_csv(file="../ProcessedData/Combined_RNAi_tracked.csv")

TT$compare <- paste0(TT$genotype,"c",TT$cross)

alphacod <- c(paste0("a",letters), paste0("b",letters), paste0("c",letters))

TTcode <- tibble("compare" = unique(TT$compare),  "code"= alphacod[1:length(unique(TT$compare))])

TT <- left_join(TT, TTcode, by="compare")

TT$code <- as.factor(TT$code)

mod <- aov(incapacitation ~ code, data= TT)

#set up comparisons
g1 <- unique(TT[TT$backg_id=="attP2","genotype"]) %>% drop_na

g2 <- unique(TT[TT$backg_id=="attP40","genotype"]) %>% drop_na

pan_g <- "3954"
neuro_g <- "25750"

b1 <- "36303"
b2 <- "36304"

pan1_compare <- TTcode[TTcode$compare == paste0(b1,"c",pan_g),"code"]
pan2_compare <- TTcode[TTcode$compare == paste0(b2,"c",pan_g),"code"]

neuro1_compare <- TTcode[TTcode$compare == paste0(b1,"c",neuro_g),"code"]
neuro2_compare <- TTcode[TTcode$compare == paste0(b2,"c",neuro_g),"code"]

g1_pan <- TTcode[TTcode$compare %in% paste0(g1$genotype,"c",pan_g),"code"]
g2_pan <- TTcode[TTcode$compare %in% paste0(g2$genotype,"c",pan_g),"code"]

g1_neuro <- TTcode[TTcode$compare %in% paste0(g1$genotype,"c",neuro_g),"code"]
g2_neuro <- TTcode[TTcode$compare %in% paste0(g2$genotype,"c",neuro_g),"code"]

compare_genes <- c(paste0(g1_pan$code, " - ", pan1_compare, "= 0"),
                   paste0(g2_pan$code, " - ", pan2_compare, "= 0"),
                   paste0(g1_neuro$code, " - ", neuro1_compare, "= 0"),
                   paste0(g2_neuro$code, " - ", neuro2_compare, "= 0"))

#set up ones we want
SS <- summary(glht(mod, linfct = mcp(code = compare_genes)))

out_ps <- tibble("test"= names(SS$test$coefficients),
                 "Estimate" = SS$test$coefficients,
                 "SE" = SS$test$sigma,
                 "tval" = SS$test$tstat,
                 "Pval" = SS$test$pvalues
                 )

spl <- str_split(out_ps$test, "-", simplify=TRUE)
spl[,1] <- str_trim(spl[,1])
spl[,2] <- str_trim(spl[,2])

spl <- as.data.frame(spl)
colnames(spl) <- c("code", "g2")

spl_genotype <- left_join(spl, TTcode, by=c("code"))

colnames(spl) <- c("g1", "code")
spl_cross <- left_join(spl, TTcode, by=c("code"))

out_ps$cross <- spl_genotype$compare
out_ps$base <- spl_cross$compare

gid <- read_csv(file="../ProcessedData/RNAi_gene_number.csv")
gid$genotype <- as.character(gid$genotype)
ss <- str_split(out_ps$cross, "c", simplify=TRUE)
out_ps$genotype <- ss[,1]
out_ps$type <- ss[,2]
out_ps <- left_join(out_ps, gid, by="genotype")
out_ps$type[out_ps$type=="3954"] <- "all"
out_ps$type[out_ps$type=="25750"] <- "neuron"



write_csv(out_ps, file="../ProcessedData/RNAi_model_output.csv")

