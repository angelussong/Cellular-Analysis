rm(list = ls())
FSI=read.csv("/Users/hsong1/Desktop/fs_internrn_neuron/fsi1.csv",header=TRUE)
library(partykit)
library(party)
library(relaimpo)
library(mlbench)
library(caret)
library(lars)
library(glmnet)
library(lpSolve)
#How to optimize parameters now that we have the fitting functions for the variables. 
#In all neurons, the amplitude is only dependent on the gGABA so we need to solve 
# the 1-D equation for the conductance value using the target and the intercept. 

#Then since we value the decay and halfwidth (shape) more than rise, we could use the decay
# and half width as constraints, and the rise function minus the target value as the target. Also remember, tau_rise
# and tau_decay should both be above 0. 
num <- as.matrix()
x <- as.matrix(FSI[,1:3])
amp <- as.matrix(FSI[,4])
rise <- as.matrix(FSI[,5])
decay <- as.matrix(FSI[,6])
hfw <-as.matrix(FSI[,7])



# lasso

fit_amp <- lars(x, amp, type="lasso",intercept = TRUE)
#summary(fit_amp)
best_step0 <- fit_amp$df[which.min(fit_amp$RSS)]
predictions_lasso0 <- predict(fit_amp, x, s=best_step0, type="fit")$fit
mse_lasso0 <- mean((amp - predictions_lasso0)^2)
cv.fit_amp <- cv.glmnet(x,amp, alpha=1)
#plot(cv.fit_amp)
#(best.lambda1 <- cv.fit_lasso1$lambda.min)
#plot(fit_amp$lambda)
coef(cv.fit_amp)

#fit_lasso1 is the fitting function for rise
fit_lasso1 <- lars(x, rise, type="lasso",intercept = TRUE)
#summary(fit_lasso1)
best_step1 <- fit_lasso1$df[which.min(fit_lasso1$RSS)]
predictions_lasso1 <- predict(fit_lasso1, x, s=best_step1, type="fit")$fit
mse_lasso1 <- mean((rise - predictions_lasso1)^2)
#print(mse_lasso1)
#coef(fit_lasso1)
cv.fit_lasso1 <- cv.glmnet(x,rise, alpha=1)
#plot(cv.fit_lasso1)
#(best.lambda1 <- cv.fit_lasso1$lambda.min)
#plot(fit_lasso1$lambda)
coef(cv.fit_lasso1)

#fit_lasso2 is the fitting function for decay
fit_lasso2 <- lars(x, decay, type="lasso",intercept = TRUE)
#summary(fit_lasso2)
best_step2 <- fit_lasso2$df[which.min(fit_lasso2$RSS)]
predictions_lasso2 <- predict(fit_lasso2, x, s=best_step2, type="fit")$fit
mse_lasso2 <- mean((decay - predictions_lasso2)^2)
#print(mse_lasso2)
#coef(fit_lasso2)
cv.fit_lasso2 <- cv.glmnet(x, decay, alpha=1)
#plot(cv.fit_lasso2)
#(best.lambda <- cv.fit_lasso2$lambda.min)
#plot(fit_lasso2$lambda)
coef(cv.fit_lasso2)

#fit_lasso3 is the fitting function for half width
fit_lasso3 <- lars(x, hfw, type="lasso",intercept = TRUE)
#summary(fit_lasso3)
best_step3 <- fit_lasso3$df[which.min(fit_lasso3$RSS)]
predictions_lasso3 <- predict(fit_lasso3, x, s=best_step3, type="fit")$fit
mse_lasso3 <- mean((hfw - predictions_lasso3)^2)
#print(mse_lasso3)
#coef(fit_lasso3)
cv.fit_lasso3 <- cv.glmnet(x, hfw, alpha=1)
#plot(cv.fit_lasso3)
#(best.lambda <- cv.fit_lasso2$lambda.min)
#plot(fit_lasso3$lambda)
coef(cv.fit_lasso3)

tar_r=11.362
tar_d=5.438
tar_h=5.409

#solve for gGABA, 1-D equation
tar_amp=23.88
Coef_amp1=coef(cv.fit_amp)[4]
b_amp<-tar_amp-coef(cv.fit_amp)[1]
Ma_amp <- matrix(c(Coef_amp1),1,1,byrow=TRUE)
colnames(Ma_amp)<-paste0('g',1)
Ma_bamp<-c(b_amp)
gGABA <- solve(Ma_amp,Ma_bamp)


#Linear Optimization
Coef_gr=coef(cv.fit_lasso1)[4]
tar_rf = tar_r-coef(cv.fit_lasso1)[4]*gGABA
tar_df = tar_d-coef(cv.fit_lasso2)[4]*gGABA
tar_hf = tar_h-coef(cv.fit_lasso3)[4]*gGABA

Coef_r1 <- coef(cv.fit_lasso1)[2]
Coef_r2 <- coef(cv.fit_lasso1)[3]
Coef_d1 <- coef(cv.fit_lasso2)[2]
Coef_d2 <- coef(cv.fit_lasso2)[3]
Coef_h1 <- coef(cv.fit_lasso3)[2]
Coef_h2 <- coef(cv.fit_lasso3)[3]
library(lpSolve)
objective.in <- c(Coef_r1, Coef_r2)
const.mat <- matrix(c(1,0,0,1,Coef_d1,Coef_d2,Coef_h1,Coef_h2), nrow=4, byrow=TRUE)
const.rhs <- c(1,2,tar_df,tar_hf)
const.dir <- c(">=",">=",">=",">=")
optimum <- lp(direction="min",objective.in, const.mat,const.dir,const.rhs)
optimum$solution
optimum$objval

#This section is to solve the equations for taus
#library(matlib)
#A <- matrix(c(Coef_11, Coef_12,
#              Coef_21, Coef_22), 2, 2, byrow=TRUE)
#colnames(A) <- paste0('x', 1:2)
#b_1=tar_d-coef(cv.fit_lasso2)[1]
#b_2=tar_h-coef(cv.fit_lasso3)[1]
#b <- c(b_1, b_2)
#showEqn(A, b)
#solve(A, b)
gGABA
Coef_d1*optimum$solution[1]+Coef_d2*optimum$solution[2]
Coef_h1*optimum$solution[1]+Coef_h2*optimum$solution[2]
optimum$solution
optimum$objval
