library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
source("../../Functions/ggplot_theme.R")

therm.dat <- read.csv(file ="../ProcessedData/Thermotolerance_test.csv", header=TRUE, stringsAsFactors = FALSE)

#drop NAs
therm.dat <- therm.dat[-which(is.na(therm.dat$Chamber)),]

therm.dat<-therm.dat[-which(is.na(therm.dat$Thermotolerance)),]

#change chamber number
unique(therm.dat$Chamber)

therm.dat$Chamber <- therm.dat$Chamber+1

#group number added

therm.dat$Group <- NA

therm.dat.new <- therm.dat[0,]


rrs <- unique(therm.dat$PatRIL)

for(ril in rrs)
{
tester <- subset(therm.dat, PatRIL == ril)

dd <- which(diff(tester$Chamber) < 0)

for(ii in 1:length(dd))
{
  if(ii == 1)
  {
  tester$Group[1:dd[ii]] <- 1
  }else{
    tester$Group[(dd[(ii-1)]+1):dd[ii]]<-ii
  }
  if(ii==length(dd))
    {
    tester$Group[(dd[ii]+1):nrow(tester)]<-(ii+1)
  }
}
therm.dat.new <- rbind(therm.dat.new, tester)
}



ThermTol <- read.table(file= "../ProcessedData/ThermalTol_processed.txt", sep="\t", header=TRUE)

ThermTol$group <- as.numeric(ThermTol$group)

ThermTol$uid <- paste(ThermTol$patRIL,".",ThermTol$chamber, ".", ThermTol$group,sep="")

therm.dat.new$uid <- paste(therm.dat.new$PatRIL,".",therm.dat.new$Chamber, ".", therm.dat.new$Group,sep="")

#drop problem values
length(which(therm.dat.new$Thermotolerance>700))
#only 3

therm.dat.new <- subset(therm.dat.new, Thermotolerance < 700)

all <- merge(ThermTol[,c('uid','incapacitation')],therm.dat.new[,c('uid','Thermotolerance')])

cor(all$incapacitation, all$Thermotolerance)

p1 <- ggplot(all, aes(x=incapacitation, y = Thermotolerance)) +
  geom_point(alpha=1/5, size = 0.5) +
  xlab('Code Assigned') +
  ylab('Human Assigned') + 
  my_theme

p1

ggsave(p1, file="../Plots/HeatBox_Human_Val.pdf", width=3, height=2)

colnames(all)<-c('uid','code','human')
all$Diff<-abs(all$code - all$human)

all <- all[order(-all$Diff),]
write.table(all, "../ProcessedData/Validation.txt", sep="\t", row.names=FALSE)


