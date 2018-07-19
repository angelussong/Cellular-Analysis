inputData <- read.csv("/Users/hsong1/Desktop/FEB_2018/msnHETD1_Nov4IR2g/param_summary.csv")
library(partykit)
library(party)
library(relaimpo)
library(mlbench)
library(caret)
# x is all the parameters, y is the err 5, corresponding to the delay
x <- as.matrix(inputData[,2:18])
# let's see how it fits for only five parameters

rheo<- as.matrix(inputData[,19])

n <- nrow(inputData)
test_rows <-sample(1:n,0.75*n)
train_rows <- inputData[sample(nrow(inputData),0.75*n),]

# 2:18 is the whole set of parameters
# the following is setup of the train and test data at 75% for training and 25% for testing. 

x.train <- as.matrix(train_rows[,2:18])
x.test <- x[-test_rows, ]

rheo_train<- as.matrix(train_rows[,19])

rheo_test<- rheo[-test_rows]

#First, let's run some PCA for spikelet or not 
log.x <-log(inputData[,2:18])
x.rheo <- inputData[,19]

x.pca <- prcomp(log.x,center = TRUE,scale. = TRUE)
summary(x.pca)

library(ggbiplot)
library(modeltools)
g <- ggbiplot(x.pca, obs.scale = 1, var.scale=1, groups=x.rheo, ellipse = TRUE, circle = TRUE)

g <- g+theme(legend.direction='horizontal',legend.position='top')
print(g)

#Ridge Regression creates a linear regression model 
#that is penalized with the L2-norm which is the sum of the squared coefficients.
library(glmnet)
fit_ridge <- glmnet(x.train,rheo_train,family="gaussian", alpha=1, lambda=0.001)
summary(fit_ridge)
predictions_ridge <-predict(fit_ridge,x.test,type="link")
mse_ridge <-mean((rheo_test-predictions_ridge)^2)
print(mse_ridge)

#Least Absolute Shrinkage and Selection Operator (LASSO) creates a regression model 
# that is penalized with the L1-norm which is the sum of the absolute coefficients. 
library(lars)
fit_lasso <- lars(x.train, rheo_train, type="lasso")
summary(fit_lasso)
best_step <- fit_lasso$df[which.min(fit_lasso$RSS)]
predictions_lasso <- predict(fit_lasso, x, s=best_step, type="fit")$fit
mse_lasso <- mean((rheo - predictions_lasso)^2)
print(mse_lasso)
plot(fit_lasso)
#Every step of iteration, the coefficients for fitting
coef(fit_lasso)
# cross validation
cv.fit_lasso <- cv.glmnet(x, rheo, alpha=1)
#cross validation with corresponding lambda values
plot(cv.fit_lasso)
# This is the best fitting for lambda (usually the last iteration)
(best.lambda <- cv.fit_lasso$lambda.min)
# this 
plot(fit_lasso$lambda)
write.table(predictions_lasso, "/Users/hsong1/Desktop/FEB_2018/msnHETD1_Nov4IR2g/freq_predicted.txt", sep="\t")

library(earth)
marsModel <- earth(rheo ~ ., data=inputData[,2:18])
ev <-evimp (marsModel)
plot(ev)
ev

#stepwise regression for delay
base.mod <- lm(delay ~ 1, inputData[,2:18])
all.mod <- lm(delay ~., inputData[,2:18])
stepMod <-step(base.mod, scope = list(lower=base.mod, upper=all.mod), direction="both", trace=0, steps=1000)
shortlistedVars <- names(unlist(stepMod[[1]]))
shortlistedVars <-shortlistedVars[!shortlistedVars %in% "(Intercept)"]
print(shortlistedVars)
#stepwise regression for frequency
base.mod <- lm(freq ~ 1, inputData[,2:18])
all.mod <- lm(freq ~., inputData[,2:18])
stepMod <-step(base.mod, scope = list(lower=base.mod, upper=all.mod), direction="both", trace=0, steps=1000)
shortlistedVars <- names(unlist(stepMod[[1]]))
shortlistedVars <-shortlistedVars[!shortlistedVars %in% "(Intercept)"]
print(shortlistedVars)

# Boruta importance for delay & frequency
library(Boruta)
boruta_output <-Boruta(inputData[,19] ~., data=na.omit(inputData[,2:18]), doTrace=2)
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
print(boruta_signif)
par(mar=c(1,1,1,1))
plot(boruta_output, cex.axis=1, las=2, main="Variable Importance")


#relative importance for delay & frequency
library(relaimpo)
lmMod_delay <- lm(inputData[,21] ~ ., data=na.omit(inputData[,2:18]))
relIm_delay <-calc.relimp(lmMod_delay,type="lmg",rela=TRUE)
sort(relIm_delay$lmg,decreasing=TRUE)


lmMod_fre <- lm(inputData[,19] ~ ., data=na.omit(inputData[,2:18]))
relIm_fre <-calc.relimp(lmMod_fre,type="lmg",rela=TRUE)
sort(relIm_fre$lmg,decreasing=TRUE)





