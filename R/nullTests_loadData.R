#Load in saved data from Null Distrubution Test already computed


rm(list = ls())  #clear workspace

#Add dependencies
library(ggplot2)
library(lme4)
#lme4 library for mixed model, with rat as random effect.
#if comparing models that differ in fixed effects, can't use REML. REML's method of estimate unbaises the estimation
# of the residual variance
#http://users.stat.umn.edu/~gary/classes/5303/handouts/REML.pdf

setwd("C:/Users/Jared/Dropbox/NVanalysis_data/allCh_2months_OneminWinFeats/data_R")


load("gammaNull.Rda")
load("LL_nullTest.Rda")


##### Gamma Power
hist(gammaNull$permChi,breaks = 75,main = expression("Gamma Power Model Null Distribution of" ~ chi^{2} ~ "Values"), 
     xlab=expression(chi^{2}), ylab="Frequency")   #,xlim=c(0,gammaNull$chiVal[1])) 
#abline(v = chiVal)



##### Line Length
hist(LLnull$permChi,breaks = 75,main= expression("Line Length Model Null Distribution of" ~ chi^{2} ~ "Values"), 
     xlab=expression(chi^{2}), ylab="Frequency")#,xlim=c(0,LLnull$chiVal[1])) 
#abline(v = LLnull$chiVal[1])

