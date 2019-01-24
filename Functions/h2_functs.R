h2.lme<-function(geno,pheno)
{
  
  #make sure geno is factor
  geno<-as.factor(geno)
  
  #get genotype IDs
  geno.ids<-unique(geno)
  
  #how many genotypes?
  ngeno<-length(geno.ids)
  
  #set up output vector
  H2.p<-numeric(length=ngeno)
  
  #get estimate
  afit<-lme(pheno ~ 1 , random= ~1 | geno)
  vars<-VarCorr(afit)
  
  va<-as.numeric(vars[1])
  vw<-as.numeric(vars[2])
  
  #Heritability
  H2.0 <- va/(va+vw)
  return(H2.0)
}


h2.aov<-function(geno,pheno)
{
  
  #make sure geno is factor
  geno<-as.factor(geno)
  
  #get genotype IDs
  geno.ids<-unique(geno)
  
  #how many genotypes?
  ngeno<-length(geno.ids)
  
  #set up output vector
  H2.p<-numeric(length=ngeno)
  
  #get estimate
  afit<-aov(pheno~geno)
  suma<-unlist(summary(afit))
  #get variance within groups
  sw<-suma['Mean Sq2']
  
  #get variance between groups SEE Lessels and Boag (1987)
  Ng<-length(unique(geno))
  ns<-table(geno)
  n0<-(1/(Ng-1)) * (sum(ns)-(sum(ns^2)/sum(ns)))
  
  sa<-(suma['Mean Sq1']-suma['Mean Sq2'])/n0
  
  
  #Observed Heritability
  H2.0<-sa/(sa+sw)
  return(H2.0)
}

jk.h2<-function(pheno,geno)
{
  #make sure geno is factor
  geno<-as.factor(geno)
  
  #get genotype IDs
  geno.ids<-unique(geno)
  
  #how many genotypes?
  ngeno<-length(geno.ids)
  
  #set up output vector
  H2.p<-numeric(length=ngeno)
  
  #get estimate
  afit<-aov(pheno~geno)
  suma<-unlist(summary(afit))
  #get variance within groups
  sw<-suma['Mean Sq2']
  
  #get variance between groups SEE Lessels and Boag (1987)
  Ng<-length(unique(geno))
  ns<-table(geno)
  n0<-(1/(Ng-1)) * (sum(ns)-(sum(ns^2)/sum(ns)))
  
  sa<-(suma['Mean Sq1']-suma['Mean Sq2'])/n0
  
  
  #Observed Heritability
  H2.0<-sa/(sa+sw)
  
  
  #do drop one jackknife
  for (i in 1:ngeno)
  {
    pheno.s<-pheno[geno!=geno.ids[i]]
    geno.s<-geno[geno!=geno.ids[i]]
    #Get new estimate
    afit<-aov(pheno.s~geno.s)
    suma<-unlist(summary(afit))
    #get variance within groups
    sw<-suma['Mean Sq2']
    
    #get variance between groups SEE Lessels and Boag (1987)
    Ng<-length(unique(geno.s))
    ns<-table(geno.s)
    n0<-(1/(Ng-1)) * (sum(ns)-(sum(ns^2)/sum(ns)))
    
    sa<-(suma['Mean Sq1']-suma['Mean Sq2'])/n0
    
    #New Heritability
    H2<-sa/(sa+sw)
    
    #Pseudovalue
    H2.p[i]<-(nril*H2.0) - ((nril-1)*H2)
    cat(i, '\n')
  }
  output<-list(c(mean(H2.p),sd(H2.p)))
  return(output)
}



jk.h2.lme<-function(pheno,geno)
{
  
  #make sure geno is factor
  geno<-as.factor(geno)
  
  #get genotype IDs
  geno.ids<-unique(geno)
  
  #how many genotypes?
  ngeno<-length(geno.ids)
  
  #set up output vector
  H2.p<-rep(NA, ngeno)

  
  #get estimate
  afit<-lme(pheno ~ 1 , random= ~1 | geno)
  vars<-VarCorr(afit)
  
  va<-as.numeric(vars[1])
  vw<-as.numeric(vars[2])
  
  #Heritability
  H2.0 <- va/(va+vw)
  
  
  #do drop one jackknife
  for (i in 1:ngeno)
  {
    tryCatch({
    pheno.s<-pheno[geno!=geno.ids[i]]
    geno.s<-geno[geno!=geno.ids[i]]
    #Get new estimate
    afit<-lme(pheno.s ~ 1 , random= ~1 | geno.s)
    vars<-VarCorr(afit)
    
    va<-as.numeric(vars[1])
    vw<-as.numeric(vars[2])
    
    #Heritability
    H2 <- va/(va+vw)
    
    #Pseudovalue
    H2.p[i]<-(ngeno*H2.0) - ((ngeno-1)*H2)
    cat(i, '\n')
  }, error = function(e){cat("Error",conditionMessage(e), "\n")})
  }
  
  return(H2.p)
}

quant.norm<-function(y)
{
  nn<-qqnorm(y,plot.it=FALSE)
  return(nn$x)
}






