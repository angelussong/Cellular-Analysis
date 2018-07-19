# This is a function to omit NA in certain columns.
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

Ephys_ori=read.csv("/Users/hsong1/Desktop/MSN ephys with animal names.csv",header=TRUE)
Ephys=completeFun(Ephys_ori,"Frequency.Hz.")
# str(Ephys)
Ep_D1=subset(Ephys,Ephys[,2]==1|Ephys[,2]==3)
Ep_D2=subset(Ephys,Ephys[,2]==2|Ephys[,2]==4)
library(Rmisc)
library(plyr)
library(dplyr)
library(nlme)
library(lme4)
library(lmerTest)
# as.factor read in categorical variable, not numerics

fcelltype=as.factor(Ep_D2$Cell.Type)
fAnimals=Ep_D2$Animals
fexp_freq=Ep_D2$Frequency.Hz.
#mod1=lme(fRN ~ fcelltype, random=~1|fAnimals,method="REML")
mod1 = lmer(fexp_freq ~ fcelltype + (1|fAnimals),REML=TRUE)
Anova(mod1)
anova(mod1)

mod1.fixed=gls(fexp_freq ~ fcelltype,method="REML")
anova(mod1,mod1.fixed)
rand(mod1)

mod2=lm(fexp_freq ~ fcelltype)
library(car)
Anova(mod2,type="II")
anova(mod2)
summary(mod1)
summary(mod2)
