Ephys=read.csv("/Users/hsong1/Desktop/MSN ephys with animal names.csv",header=TRUE)
Ep_D1=subset(Ephys,Ephys[,2]==1|Ephys[,2]==3)
Ep_D2=subset(Ephys,Ephys[,2]==2|Ephys[,2]==4)
library(Rmisc)
library(plyr)
library(dplyr)

Ephys$fcelltype=as.factor(Ephys$Cell.Type)
Ephys$fNeurons=as.factor(Ephys$Neurons)
Ephys$fExcFreq=as.factor(Ephys$Frequency.Hz.)

summary(mod1)
anova(mod1)

mod2=lm(FR~fgenotype+fcelltype+Rheobase+RN+True_i+fgenotype*True_i+fcelltype*True_i+
          fgenotype*Rheobase+fcelltype*Rheobase+fgenotype*RN+
          fcelltype*RN,data=Ephys)
summary(mod2)
anova(mod2)

res <- stack(data.frame(Observed = Ephys$FR))
res <- cbind(res, x = rep(Ephys$Rheobase, 2))
head(res)
require("lattice")
xyplot(values ~ x, data=res,group=ind,auto.key=TRUE)

dat_sub=subset(Ephys,Ephys$True_i==180)
dat_sub
dat_sub$Rheobase

fitted(mod2)[Ephys$True_i==180]
res <- stack(data.frame(Observed = dat_sub$FR, Predicted=fitted(mod2)[Ephys$True_i==180]))
res <- cbind(res, x = rep(dat_sub$Rheobase, 2))
head(res)
require("lattice")
xyplot(values ~ x, data=res,group=ind,auto.key=TRUE)
RN=[Ephys$RN,Ephys$RN>=0]
mean(RN)

mod_lm=lm(FR~Rheobase,data=dat_sub)
summary(mod_lm)
anova(mod_lm)
res <- stack(data.frame(Observed = dat_sub$FR, Predicted=fitted(mod_lm)))
res <- cbind(res, x = rep(dat_sub$Rheobase, 2))
require("lattice")
xyplot(dat_sub$FR ~ dat_sub$Rheobase, data = dat_sub, type = c("p","r"), col.line = "red")

dat_sub$FR
