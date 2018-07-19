rm(list = ls())

D1PCA=read.csv("/Users/hsong1/Desktop/MAR_2018/D1PCAcombined.csv",header=TRUE)

head(D1PCA,3)
D1.data <-D1PCA[,1:17]
D1.type <-D1PCA[,18]
D1.pca <- prcomp(D1.data, center=TRUE, scale.=TRUE)
print(D1.pca)
plot(D1.pca,type="l")
library(devtools)
library(ggbiplot)
g <-ggbiplot(D1.pca,obs.scale=1,var.scale=1,groups=D1.type,ellipse=TRUE,circle=TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
