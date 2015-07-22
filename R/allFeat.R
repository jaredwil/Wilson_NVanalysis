rm(list = ls())  #clear workspace

setwd("C:/Users/Jared/Dropbox/NVanalysis_data/allCh_2months_OneminWinFeats/data_R")

featData = read.table("LLEnergy_allPt.csv",header=FALSE, sep=",")
gammaData = read.table("gammaBP_allPt.csv",header=FALSE, sep=",")
areaData = read.table("AreaV2_allPt.csv",header=FALSE, sep=",")

featData = cbind(featData,gammaData[,4:5],areaData[,4:5])

names(featData) = c('Hospital','PatientID','Time','AvgLL','StdLL','AvgEnergy','StdEnergy','AvgGamma','StdGamma','AvgArea','StdArea')
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
logData = featData[idxToKeep,]         
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

####Line Length

# mod = lmer(AvgLL ~ Time + (1|PatientID) + (1|Hospital),REML=TRUE,data=newData)
# mod.null = lmer(AvgLL ~ (1|PatientID) + (1|Hospital),REML=TRUE,data=newData)
mod = lmer(AvgLL ~ Time + (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgLL ~ (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT
 
hist(residuals(mod),breaks=1000)
boxplot(residuals(mod))
plot(residuals(mod),newData$AvgLL)
plot(fitted(mod),newData$AvgLL)

plot(newData$Time,residuals(mod))#ylim=c(-10, 10))
plot(newData$Time,fitted(mod)) 


#####Energy

mod = lmer(AvgEnergy ~ Time + (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgEnergy ~ (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT

hist(residuals(mod),breaks=1000)
boxplot(residuals(mod))
plot(residuals(mod),newData$AvgEnergy)
plot(fitted(mod),newData$AvgEnergy)

plot(newData$Time,residuals(mod))#ylim=c(-10, 10))
plot(newData$Time,fitted(mod)) 


####Gamma Power

mod = lmer(AvgGamma ~ Time + (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
mod.null = lmer(AvgGamma ~ (1|PatientID) + (1|Hospital),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT

hist(residuals(mod),breaks=1000)
boxplot(residuals(mod),ylim=c(-40,75))
plot(residuals(mod),newData$AvgGamma)
plot(fitted(mod),newData$AvgGamma)
plot(fitted(mod),residuals(mod),ylim=c(-50,50))  #hertroskadisticity

plot(newData$Time,residuals(mod),ylim=c(-50, 50))
plot(newData$Time,fitted(mod)) 



### Area 
mod = lmer(AvgArea ~ Time + (1|PatientID) ,REML=FALSE,data=newData)
mod.null = lmer(AvgArea ~ (1|PatientID) ,REML=FALSE,data=newData)
anova(mod.null,mod) #LRT

hist(residuals(mod),breaks=1000)
boxplot(residuals(mod),ylim=c(-40,75))
plot(residuals(mod),newData$AvgGamma)
plot(fitted(mod),newData$AvgGamma)
plot(fitted(mod),residuals(mod),ylim=c(-50,50))  #hertroskadisticity

plot(newData$Time,residuals(mod),ylim=c(-50, 50))
plot(newData$Time,fitted(mod)) 

# #random slope and intercept 
# #resulting anova p = 0.172
# mod = lmer(AvgLL ~ Time + (1|PatientID) + (Time-1|PatientID),REML=FALSE,data=newData)
# mod.null = lmer(AvgLL ~ (1|PatientID) + (Time-1|PatientID),REML=FALSE,data=newData)
# #resulting anova p = 0.1732  very similar to above (I think the syntax is saying the same thing)
# mod = lmer(AvgLL ~ Time + (1+Time|PatientID),REML=FALSE,data=newData)
# mod.null = lmer(AvgLL ~ (1+Time|PatientID),REML=FALSE,data=newData)

mod = lmer(AvgLL ~ Time + (1|PatientID),REML=FALSE,data=newData)
mod.null = lmer(AvgLL ~ (1|PatientID),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT 

#Gamma
mod = lmer(AvgGamma ~ Time + (1|PatientID) + (Time-1|PatientID),REML=FALSE,data=newData)
mod.null = lmer(AvgGamma ~ (1|PatientID) + (Time-1|PatientID),REML=FALSE,data=newData)

anova(mod.null,mod) #LRT 

hist(residuals(mod),breaks=1000)
boxplot(residuals(mod))
plot(residuals(mod),newData$AvgGamma)
plot(fitted(mod),newData$AvgGamma,ylim=c(0, 200))

plot(newData$Time,residuals(mod),ylim=c(-30, 30))
plot(newData$Time,fitted(mod)) 









