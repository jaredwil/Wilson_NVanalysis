#Create a Null Distrubution of the data using permutations to verify model results > chance
rm(list = ls())  #clear workspace

#Add dependencies
library(ggplot2)
library(lme4)
#lme4 library for mixed model, with rat as random effect.
#if comparing models that differ in fixed effects, can't use REML. REML's method of estimate unbaises the estimation
# of the residual variance
#http://users.stat.umn.edu/~gary/classes/5303/handouts/REML.pdf

setwd("C:/Users/Jared/Dropbox/NVanalysis_data/allCh_2months_OneminWinFeats/csv")

featData = read.table("LLEnergy_allPt.csv",header=FALSE, sep=",")
gammaData = read.table("gammaBP_allPt.csv",header=FALSE, sep=",")

featData = cbind(featData,gammaData[,4:5])

names(featData) = c('Hospital','PatientID','Time','AvgLL','StdLL','AvgEnergy','StdEnergy','AvgGamma','StdGamma')
#make patient ID unique since multiple per hospital
featData$PatientID = featData$Hospital*featData$PatientID 
featData$Hospital = as.factor(featData$Hospital)
featData$PatientID = as.factor(featData$PatientID)
foo = 1:nrow(featData)


#Run through this 10000 times in order to create null distribution
t = (1:max(featData$Time))/60/60/24
M = 2000
permChi <- vector(mode = "integer",length = M)
permP <- vector(mode = "integer",length = M)

for (i in 1:M ) {
  # ptm <- proc.time()  #START CLOCK
  
  permData = featData
  
  #Reassigne time values for each pt.
  for (pt in 1:14){
    timePerm = sample(t,length(t));
    permData$Time[((pt-1)*length(t)+1):(pt*length(t))] = timePerm;
  }
  
  permData= permData[permData$AvgLL!=0,]   #remove idices where avgLL is zero
  permData= permData[permData$AvgEnergy!=0,]   #remove indices where avgEnergy is zero
  permData= permData[permData$AvgGamma!=0,]   #remove indices where avgEnergy is zero
  
  mod = lmer(AvgGamma ~ Time + (1|PatientID),REML=FALSE,data=permData)
  mod.null = lmer(AvgGamma ~ (1|PatientID),REML=FALSE,data=permData)
  permP[i] = anova(mod.null,mod)$`Pr(>Chisq)`[2] #LRT
  permChi[i]   = anova(mod.null,mod)$Chisq[2]

  # proc.time() - ptm   #STOP CLOCK
  
}


#Downsample by M = 60
# i = 1
# idxToKeep = c(1,foo[1:(i+59)==(i+59)]) #take every 60th sample
#newData = featData[idxToKeep,]         
newData = featData


#Normal Data
newData= newData[newData$AvgLL!=0,]   #remove idices where avgLL is zero
newData= newData[newData$AvgEnergy!=0,]   #remove indices where avgEnergy is zero
newData= newData[newData$AvgGamma!=0,]   #remove indices where avgEnergy is zero
newData$Time = newData$Time/60/60/24  #normalize all times to be 0-1 instead of seconds


####Gamma Power TEST

mod = lmer(AvgGamma ~ Time + (1|PatientID),REML=FALSE,data=newData)
mod.null = lmer(AvgGamma ~ (1|PatientID),REML=FALSE,data=newData)
anova(mod.null,mod) #LRT
pVal     = anova(mod.null,mod)$`Pr(>Chisq)`[2] #LRT
chiVal   = anova(mod.null,mod)$Chisq[2]



# r = anova(mod.null,mod)

hist(permChi,breaks=500) 
abline(v = chiVal)


# hist(residuals(mod),breaks=1000)
# boxplot(residuals(mod),ylim=c(-40,75))
# plot(residuals(mod),newData$AvgGamma)
# plot(fitted(mod),newData$AvgGamma)
# plot(fitted(mod),residuals(mod),ylim=c(-50,50))  #hertroskadisticity
# 
# plot(newData$Time,residuals(mod),ylim=c(-50, 50))
# plot(newData$Time,fitted(mod)) 

