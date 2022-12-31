library(DSPRqtl)
library(tidyverse)
library(cowplot)
library(viridis)
source("../../Functions/mappingfunctions.R")
source(file="../../Functions/ggplot_theme.R")
theme_set(theme_cowplot())




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

compD <- tibble("Position" = c(positions$Gaxis,positions$Gaxis),
                "LOD" = c(obs.LL[,3],Phen.Ls),
                "Scan"=c(rep("Raw",nrow(obs.LL)), rep("Q6 Corrected",length(Phen.Ls))))

top <- max(obs.LL$Tolsqrtvariable)

p1 <- ggplot(compD, aes(x=Position, y=LOD, color=Scan)) +
  geom_rect(xmin=66.3,ymin=-10,xmax=120.3,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_rect(xmin=174,ymin=-10,xmax=221,ymax=top+10,fill='aliceblue',color='aliceblue')+
  scale_x_continuous(breaks=c(0,66.3,120,174,221,277,33.15,98.3,145,205,249), 
                     labels=c(0,"66  0",54,"108  0",47,103,'\nX','\n2L','\n2R','\n3L','\n3R'),limits = c(-0.01,max(compD$Position)+25))+
  geom_line(alpha=0.6) +
  theme(axis.ticks.x=element_blank())+
  xlab("Position (cM)") +
  scale_color_viridis(discrete = TRUE)+
  my_theme

p1




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



compD <- tibble("Position" = c(positions$Gaxis,positions$Gaxis),
                "LOD" = c(obs.LL[,3],Phen.Ls),
                "Scan"=c(rep("Raw",nrow(obs.LL)), rep("Subpop",length(Phen.Ls))))

compD$Scan <- factor(compD$Scan, levels = c("Subpop","Raw"))

top <- max(obs.LL$Tolsqrtvariable)

p2 <- ggplot(compD, aes(x=Position, y=LOD, color=Scan)) +
  geom_rect(xmin=66.3,ymin=-10,xmax=120.3,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_rect(xmin=174,ymin=-10,xmax=221,ymax=top+10,fill='aliceblue',color='aliceblue')+
  scale_x_continuous(breaks=c(0,66.3,120,174,221,277,33.15,98.3,145,205,249), 
                     labels=c(0,"66  0",54,"108  0",47,103,'\nX','\n2L','\n2R','\n3L','\n3R'),limits = c(-0.01,max(compD$Position)+25))+
  geom_line(alpha=0.6) +
  theme(axis.ticks.x=element_blank())+
  xlab("Position (cM)") +
  scale_color_viridis(discrete = TRUE)+
  my_theme

p2

pall <- plot_grid(p2,p1, labels=c("A.","B."), nrow=2)

ggsave(pall, width = 6.5, height = 5, filename = "../Plots/TT_additional.pdf")

