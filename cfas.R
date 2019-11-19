#setup, load libraries
options(scipen=999) #disables scientific notation
options(max.print=1000000) #enable long outputs
library(semTools)

#load & setup data
datatg <- read.csv("datatg.csv");datatg <- datatg[,2:17]
datato <-datatg; datato[,1:16] <- lapply(datatg[,1:16], ordered)

######
######
#cfas#
######
######

#############################
####associationist models####
#############################
####1 general lv (atti)
aso1 <- '
atti =~
risk01r + risk02r + risk03i + risk04i +
bene01d + bene02r + bene03r + bene04d +
trus01r + trus02r + trus03d + trus04d +
acce01r + acce02r + acce03d + acce04d
'
####4 first order lvs (risk,bene,trus,acce), 1 second order lv (atti) 
aso2 <- '
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
trus =~ trus01r + trus02r + trus03d + trus04d
acce =~ acce01r + acce02r + acce03d + acce04d
atti =~ risk + bene + trus + acce
'
####aso2 w bifactor method lv
aso3 <- '
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
trus =~ trus01r + trus02r + trus03d + trus04d
acce =~ acce01r + acce02r + acce03d + acce04d
atti =~ risk + bene + trus + acce
meth =~
-1*risk01r + -1*risk02r + 1*risk03i + 1*risk04i +
 1*bene01d + -1*bene02r + -1*bene03r + 1*bene04d +
-1*trus01r + -1*trus02r + 1*trus03d + 1*trus04d +
-1*acce01r + -1*acce02r + 1*acce03d + 1*acce04d
meth ~~ 0*atti
'
#####################
####causal models####
#####################
####4 first order lvs (risk, bene, trus, acce)
cau1 <- '
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
trus =~ trus01r + trus02r + trus03d + trus04d
acce =~ acce01r + acce02r + acce03d + acce04d
'
####cau1 bifactor with method lv
cau2 <-'
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
trus =~ trus01r + trus02r + trus03d + trus04d
acce =~ acce01r + acce02r + acce03d + acce04d
meth =~
-1*risk01r + -1*risk02r + 1*risk03i + 1*risk04i +
 1*bene01d + -1*bene02r + -1*bene03r + 1*bene04d +
-1*trus01r + -1*trus02r + 1*trus03d + 1*trus04d +
-1*acce01r + -1*acce02r + 1*acce03d + 1*acce04d
meth ~~ 0*risk;meth ~~ 0*bene;meth ~~ 0*trus;meth ~~ 0*acce
'
##############################
####trust-heuristic models####
##############################
#### 2 oblique lvs (atti, trus)
heu1 <- '
atti =~
risk01r + risk02r + risk03i + risk04i +
bene01d + bene02r + bene03r + bene04d +
acce01r + acce02r + acce03d + acce04d
trus =~ trus01r + trus02r + trus03d + trus04d
'
####3 first order lvs (risk, bene, acce), 1 second order lv (atti), and 1 oblique first order lv (trus)
heu2 <- '
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
acce =~ acce01r + acce02r + acce03d + acce04d
atti =~ risk + bene + acce
trus =~ trus01r + trus02r + trus03d + trus04d
'
####heu2 bifactor with method lv 
heu3 <- '
risk =~ risk01r + risk02r + risk03i + risk04i
bene =~ bene01d + bene02r + bene03r + bene04d
acce =~ acce01r + acce02r + acce03d + acce04d
atti =~ risk + bene + acce
trus =~ trus01r + trus02r + trus03d + trus04d
meth =~
  -1*risk01r + -1*risk02r + 1*risk03i + 1*risk04i +
  1*bene01d + -1*bene02r + -1*bene03r + 1*bene04d +
  -1*trus01r + -1*trus02r + 1*trus03d + 1*trus04d +
  -1*acce01r + -1*acce02r + 1*acce03d + 1*acce04d
meth ~~ 0*atti; meth~~0*trus
'

######################
####running models####
######################
fitaso1 <-cfa(aso1, data=datatg, estimator="MLR")
fitaso2 <-cfa(aso2, data=datatg, estimator="MLR")
fitaso3 <-cfa(aso3, data=datatg, estimator="MLR")
fitcau1 <-cfa(cau1, data=datatg, estimator="MLR")
fitcau2 <-cfa(cau2, data=datatg, estimator="MLR")
fitheu1 <-cfa(heu1, data=datatg, estimator="MLR")
fitheu2 <-cfa(heu2, data=datatg, estimator="MLR")
fitheu3 <-cfa(heu3, data=datatg, estimator="MLR")
#alternative estimator
fitaso1 <-cfa(aso1, data=datato, estimator="DWLS")
fitaso2 <-cfa(aso2, data=datato, estimator="DWLS")
fitaso3 <-cfa(aso3, data=datato, estimator="DWLS")
fitcau1 <-cfa(cau1, data=datato, estimator="DWLS")
fitcau2 <-cfa(cau2, data=datato, estimator="DWLS")
fitheu1 <-cfa(heu1, data=datato, estimator="DWLS")
fitheu2 <-cfa(heu2, data=datato, estimator="DWLS")
fitheu3 <-cfa(heu3, data=datato, estimator="DWLS")

#########################
####results:model fit####
#########################
fitMeasures(fitaso1, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitaso2, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitaso3, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitcau1, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitcau2, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitheu1, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitheu2, c("chisq","df","pvalue","cfi","tli"))
fitMeasures(fitheu3, c("chisq","df","pvalue","cfi","tli"))
#or
fitMeasures(fitaso1, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitaso2, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitaso3, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitcau1, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitcau2, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitheu1, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitheu2, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))
fitMeasures(fitheu3, c("chisq.scaled","df","pvalue.scaled","cfi.scaled","tli.scaled"))

#####################################################
####results: standardized loadings & correlations####
#####################################################
standardizedSolution(fitaso1)
standardizedSolution(fitaso2)
standardizedSolution(fitaso3)
standardizedSolution(fitcau1)
standardizedSolution(fitcau2)
standardizedSolution(fitheu1)
standardizedSolution(fitheu2)
standardizedSolution(fitheu3)

##################################
####results: model comparisons####
##################################
anova(fitaso2,fitcau1,fitheu2)
anova(fitaso3,fitcau2,fitheu3)
net(fitaso2,fitcau1)
net(fitaso3,fitcau2)

