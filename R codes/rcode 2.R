FI_all=read.csv("/Users/cweaver2/Dropbox/CHDI_hsongsim/Goodliffe et al/FI_all.csv",header=TRUE)
str(FI_all)
FI_all$fgenotype=as.factor(FI_all$Genotype)
FI_all$fcelltype=as.factor(FI_all$Celltype)

mod1=lm(FR~fgenotype+fcelltype+True_i+fgenotype*True_i+fcelltype*True_i,data=FI_all)
summary(mod1)
anova(mod1)

mod2=lm(FR~fgenotype+fcelltype+Rheobase+RN+True_i+fgenotype*True_i+fcelltype*True_i+
          fgenotype*Rheobase+fcelltype*Rheobase+fgenotype*RN+
          fcelltype*RN,data=FI_all)
summary(mod2)
anova(mod2)

mod3=lm(FR~fgenotype+fcelltype+Rheobase+RN+True_i+fgenotype*True_i+fcelltype*True_i+
          fgenotype*Rheobase+fcelltype*Rheobase+fgenotype*fcelltype*Rheobase+fgenotype*RN+
          fcelltype*RN,data=FI_all)
summary(mod3)
anova(mod3)

plot(FI_all$FR,FI_all)

plot(mod1)
