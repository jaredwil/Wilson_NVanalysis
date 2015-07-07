rm(list = ls())  #clear workspace

setwd("C:/Users/Jared/Dropbox/NVanalysis_data/allCh_2months_OneminWinFeats")

featData = read.table("all_ptFeat.csv",header=FALSE, sep=",")
names(featData) = c('Hospital','PatientID','Time','AvgLL','StdLL','AvgEnergy','StdEnergy')
#make patient ID unique since multiple per hospital
featData$PatientID = featData$Hospital*featData$PatientID 
featData$Hospital = as.factor(featData$Hospital)
featData$PatientID = as.factor(featData$PatientID)
foo = 1:nrow(featData)

#Downsample by M = 60
i = 1
idxToKeep = c(1,foo[1:(i+59)==(i+59)]) #take every 60th sample
#newData = featData[idxToKeep,]         
newData = featData
LogData = featData[idxToKeep,]         
fullData = featData



#Take log of avgLL and avgEnergy
logData= logData[logData$AvgLL!=0,]   #remove idices where avgLL is zero
logData= logData[logData$AvgEnergy!=0,]   #remove indices where avgEnergy is zero
logData$Time = logData$Time/60/60/24  #normalize all times to be 0-1 instead of seconds
logData$AvgLL = log(logData$AvgLL)
logData$AvgEnergy = log(logData$AvgEnergy)

#Normal Data
newData= newData[newData$AvgLL!=0,]   #remove idices where avgLL is zero
newData= newData[newData$AvgEnergy!=0,]   #remove indices where avgEnergy is zero
newData$Time = newData$Time/60/60/24  #normalize all times to be 0-1 instead of seconds

###################### Visualize  #################################
# library(ggplot2)
# 
# a <- ggplot(data = newData, aes(x = Time, y = AvgLL))
# a <- a + geom_point()
# a <- a + facet_wrap(PatientID)
# a <- a + xlab("Time") + ylab("AvgLL") + ggtitle("Overview")
# a
# 
# b <- ggplot(data = newData, aes(x = Time, y = AvgEnergy))
# b <- b + geom_point()
# b <- b + facet_wrap(PatientID)
# b <- b + xlab("Time") + ylab("AvgEnergy") + ggtitle("Overview")
# b


##################### Generate Mixed Models #################
library(lme4)
#lme4 library for mixed model, with rat as random effect.
#if comparing models that differ in fixed effects, can't use REML. REML's method of estimate unbaises the estimation
# of the residual variance
#http://users.stat.umn.edu/~gary/classes/5303/handouts/REML.pdf

#mod = lmer(AvgLL ~ Time + (1|PatientID) + (1|Hospital),REML=TRUE,data=newData)
#mod.null = lmer(AvgLL ~ (1|PatientID) + (1|Hospital),REML=TRUE,data=newData)
mod = lmer(AvgLL ~ Time + (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgLL ~ (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT

hist(residuals(mod),breaks=1000)
boxplot(residuals(mod))
plot(residuals(mod),newData$AvgLL)
plot(fitted(mod),newData$AvgLL)

#Energy
mod = lmer(AvgEnergy ~ Time + (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgEnergy ~ (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT


#random slope and intercept
#mod = lmer(AvgLL ~ Time + (1|PatientID) + (Time-1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
#mod.null = lmer(AvgLL ~ (1|PatientID) + (Time-1|PatientID) + (1|Hospital),REML=FALSE,data=newData)

mod = lmer(AvgLL ~ Time + (1|PatientID) + (Time-1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgLL ~ (1|PatientID) + (Time-1|PatientID) + (1|Hospital),REML=FALSE,data=newData)

anova(mod.null,mod) #LRT 

