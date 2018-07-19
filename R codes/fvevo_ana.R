inputData <- read.csv("/Users/hsong1/Desktop/FEB_2018/WTD1_Apr26IR2c.csv")
library(partykit)
library(party)
library(relaimpo)
library(mlbench)
library(caret)
# x is all the parameters, y is the err 5, corresponding to the delay
x <- as.matrix(inputData[,2:18])
# let's see how it fits for only five parameters

peak<- as.matrix(inputData[,19])
spk_let <- as.matrix(inputData[,20])
delay <- as.matrix(inputData[,21])
freq <- as.matrix(inputData[,22])
n <- nrow(inputData)
test_rows <-sample(1:n,0.75*n)
train_rows <- inputData[sample(nrow(inputData),0.75*n),]

# 2:18 is the whole set of parameters

x.train <- as.matrix(train_rows[,2:18])
x.test <- x[-test_rows, ]

peak_train<- as.matrix(train_rows[,19])
spk_let_train <- as.matrix(train_rows[,20])
delay_train <- as.matrix(train_rows[,21])
freq_train <- as.matrix(train_rows[,22])

peak_test<- as.matrix(train_rows[,19])
spk_let_test <- as.matrix(train_rows[,20])
delay_test <- as.matrix(train_rows[,21])
freq_test <- as.matrix(train_rows[,22])

#Ridge Regression creates a linear regression model 
#that is penalized with the L2-norm which is the sum of the squared coefficients.
library(glmnet)
fit_ridge <- glmnet(x.train,y.train,family="gaussian", alpha=1, lambda=0.001)
summary(fit_ridge)
predictions_ridge <-predict(fit_ridge,x.test,type="link")
mse_ridge <-mean((y.test-predictions_ridge)^2)
print(mse_ridge)

#Least Absolute Shrinkage and Selection Operator (LASSO) creates a regression model 
# that is penalized with the L1-norm which is the sum of the absolute coefficients. 
library(lars)
fit_lasso <- lars(x.train,delay_train, type="lasso")
summary(fit_lasso)
best_step <- fit_lasso$df[which.min(fit_lasso$RSS)]
predictions_lasso <- predict(fit_lasso, x.test, s=best_step, type="fit")$fit
mse_lasso <- mean((delay_test - predictions_lasso)^2)
print(mse_lasso)
plot(fit_lasso)
#Every step of iteration, the coefficients for fitting
coef(fit_lasso)
# cross validation
cv.fit_lasso <- cv.glmnet(x.train, delay_train, alpha=1)
#cross validation with corresponding lambda values
plot(cv.fit_lasso)
# This is the best fitting for lambda (usually the last iteration)
(best.lambda <- cv.fit_lasso$lambda.min)
# this 
plot(fit_lasso$lambda)

#Elastic Net creates a regression model that is penalized with both the L1-norm and L2-norm. 

library(glmnet)
fit_elas <- glmnet(x.train, delay_train, family="gaussian", alpha=0.5, lambda=0.001)
summary(fit_elas)
predictions_elas <- predict(fit_elas, x.test, type="link")
mse_elas <- mean((delay_test - predictions_elas)^2)
print(mse_elas)

library(earth)
marsModel <- earth(inputData[,22] ~ ., inputData[,2:18])
ev <-evimp (marsModel)
plot(ev)
library(relaimpo)
lmMod <-lm(inputData[,21]~., data=inputData[,2:18])
relImportance <- calc.relimp(lmMod,type="lmg",rela=TRUE)
sort(relImportance$lmg,decreasing=TRUE)
#stepwise regression for delay
base.mod <- lm(y ~ 1, inputData[,3:19])
all.mod <- lm(y ~., inputData[,3:19])
stepMod <-step(base.mod, scope = list(lower=base.mod, upper=all.mod), direction="both", trace=0, steps=1000)
shortlistedVars <- names(unlist(stepMod[[1]]))
shortlistedVars <-shortlistedVars[!shortlistedVars %in% "(Intercept)"]
print(shortlistedVars)
#stepwise regression for frequency
base.mod <- lm(err_fre ~ 1, inputData[,3:19])
all.mod <- lm(err_fre ~., inputData[,3:19])
stepMod <-step(base.mod, scope = list(lower=base.mod, upper=all.mod), direction="both", trace=0, steps=1000)
shortlistedVars <- names(unlist(stepMod[[1]]))
shortlistedVars <-shortlistedVars[!shortlistedVars %in% "(Intercept)"]
print(shortlistedVars)

# Boruta importance for delay & frequency
library(Boruta)
boruta_output <-Boruta(inputData[,21] ~., data=na.omit(inputData[,3:19]), doTrace=2)
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
print(boruta_signif)
par(mar=c(1,1,1,1))
plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

borutafre_output <-Boruta(inputData[,24] ~., data=na.omit(inputData[,3:19]), doTrace=2)
borutafre_signif <- names(borutafre_output$finalDecision[borutafre_output$finalDecision %in% c("Confirmed", "Tentative")]) 
print(borutafre_signif)
par(mar=c(1,1,1,1))
plot(borutafre_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

#relative importance for delay & frequency
library(relaimpo)
lmMod_delay <- lm(inputData[,21] ~ ., data=na.omit(inputData[,3:19]))
relIm_delay <-calc.relimp(lmMod_delay,type="lmg",rela=TRUE)
sort(relIm_delay$lmg,decreasing=TRUE)


lmMod_fre <- lm(inputData[,21] ~ ., data=na.omit(inputData[,3:19]))
relIm_fre <-calc.relimp(lmMod_fre,type="lmg",rela=TRUE)
sort(relIm_fre$lmg,decreasing=TRUE)

#use cross validation to compute importance 
library(earth)
mars_delay <- earth(inputData[,24] ~ ., data=na.omit(inputData[,3:19])) # build model
ev_delay <- evimp (mars_delay)
plot(ev,las=2)

mars_fre <- earth(inputData[,25] ~ ., data=na.omit(inputData[,3:19]))
ev_Fre <-evimp(mars_fre)

mars_shape1 <- earth(inputData[,22] ~., data=na.omit(inputData[,3:19]))
ev_shape1 <-evimp(mars_shape)

mars_shape2 <- earth(inputData[,26] ~., data=na.omit(inputData[,3:19]))
ev_shape2 <-evimp(mars_shape)

mars_amp <-earth(inputData[,23] ~., data=na.omit(inputData[,3:19]))
ev_amp <-evimp(mars_amp)

# nonlinear regression try
x <- as.matrix(inputData[,3:19])
y <- as.matrix(inputData[,24])
err_fre <- as.matrix(inputData[,25])
n <- nrow(inputData)
test_rows <-sample(1:n,0.75*n)
train_rows <- inputData[sample(nrow(inputData),0.75*n),]
x.train <- as.matrix(train_rows[,3:19])
x.test <- x[-test_rows, ]
y.train <- as.matrix(train_rows[,24])
y.test <-y[-test_rows]
fre.train <- as.matrix(train_rows[,25])
fre.test <-err_fre[-test_rows]
library(kernlab)
fit <- ksvm(y.train~., x.train)
summary(fit)
predictions <- predict(fit, x.train)
mse <- mean((longley$Employed - predictions)^2)
print(mse)
