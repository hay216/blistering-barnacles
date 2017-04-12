rm(list=ls())



# source("I:\\ACERA\\Biofouling\\Biofouling\\papers\\R script New_testLoc.txt", echo=TRUE)





############## Data #############

# samplesdata<-read.table("I:/ACERA/Biofouling/Biofouling/data/Data/samples data/samples.csv", sep=",")

# propzerodata<-read.table("I:\\ACERA\\Biofouling\\Biofouling\\data\\Data\\samples data\\aa.csv", sep=",")

# covars<-read.table("I:\\ACERA\\Biofouling\\Biofouling\\data\\Data\\samples data\\covars.csv", sep=",")

samplesdata <- read.csv("../data-raw/samples.csv")
propzerodata <- read.csv("../data-raw/aa.csv")
covars <- read.csv("../data-raw/covars.csv")


 nonzerosamples<-samplesdata[samplesdata$wetWeight>.5, ]

# nonzerosamples$logwetWeight<-log(nonzerosamples$wetWeight)



## There are two outliers, so the analysis of LMM is done by removing those outliers


 nonzerosamples1<-subset(nonzerosamples, boatID!=24 & boatID!=51)





########  samples only from Wet Hull ####################

Loc<-c("Anchor", "Deck", "Dry Hull", "Fishing gear", "Internal spaces", "Propeller", "Rudder", "Wet Hull")

samplesdata_WH<-samplesdata[samplesdata$Location==Loc[8], ]

nonzerosamples_WH<-samplesdata_WH[samplesdata_WH$wetWeight>.5, ]

# nonzerosamples_WH$logwetWeight<-log(nonzerosamples_WH$wetWeight)



# There are two outliers, so the analysis of LMM is done by removing those outliers


nonzerosamples1_WH<-subset(nonzerosamples_WH, boatID!=24 & boatID!=51)


###############################   HA Only  ###############################################



testLoc<-c("IB", "HA", "PJ", "PB", "HD", "HAH", "HP", "HH", "HK")



samplesdata_HA<-samplesdata[samplesdata$LocID==testLoc[2], ]





nonzerosamples_HA<-samplesdata_HA[samplesdata_HA$wetWeight>.5, ]



nonzerosamples1_HA<-subset(nonzerosamples_HA, boatID!=24 & boatID!=51)


##################  LocID (IB, HA, .., HK) ##################################################

# samplesdata_testLoc<-subset(samplesdata,LocID %in% c(testLoc[1], testLoc[2], testLoc[3], 
# testLoc[4], testLoc[5], testLoc[6], testLoc[7], testLoc[8], testLoc[9]))


samplesdata_testLoc<-subset(samplesdata,LocID %in% testLoc[1:9])

nonzerosamples_testLoc<-subset(nonzerosamples,LocID %in% testLoc[1:9])

##########################  propzerodata_IB<-propzerodata[propzerodata$LocID==testLoc[1], ]
##########################  


propzerodata_testLoc<-subset(propzerodata, LocID %in% testLoc[1:9])





################ Graphical Summary ###########
################### Histogram Figure 1 ##############

# Three Histograms summarizing the number of days since the vessel’s last trip, 
# the number of days since it was last cleaned or antifouled, and the median 
# number of trips per annum.


hist(covars$days1, prob=F, xlab="Days since last trip", ylab="Frequency", main="")
	text(0,36,"a")

hist(covars$days2, prob=F, xlab="Days since last cleaned or antifouled", ylab="Frequency", main="")
	text(0,24,"b")

hist(covars$midTrips, prob=F, xlab="Median number of trips per annum", ylab="Frequency", main="")
	text(0,37,"c") 

########### Histogram Figure 2 ####################

# Histograms of sample wet weight (grams) for the whole sample, and the log of sample 
# wet weight (grams) for the whole sample.

# Histogram of wet weight shows that samples are zero-inflated. However, histogram of 
# log transformed data shows that samples are approximately normally distributed.


hist(samplesdata$wetWeight, prob=F, xlab="Wet Weight of fouling biomass", ylab="Frequency", main="")
	text(0,1510,"a")

hist(log(samplesdata$wetWeight), prob=F, xlab="log Wet Weight of fouling biomass", ylab="Frequency", main="")
	text(0,370,"b")


######### Histogram Figure 3 #####################

# Histograms of sample wet weight (grams) for the non-zero (biomass) sample, and the log of 
# sample wet weight (grams) for the non-zero (biomass) samples.

# The wet weight fouling biomass of non-zero samples is approximately log normally distributed 
# but the log of these non-zero biomasses is normally distributed.  

nonzerosamples<-samplesdata[samplesdata$wetWeight>.5, ]

hist(nonzerosamples$wetWeight, prob=F, xlab="Wet weight of fouling biomass of 
	non-zero samples", ylab="Frequency", main="" ) 

text(0,770, "a")

hist(log(nonzerosamples$wetWeight), prob=F, xlab="Log of the Wet weight of 
	fouling biomass of non-zero samples", ylab="Frequency", main="" ) 

text(0,200, "b")


################### Smooth Plot Figures 4-6 #################

# Smooth plots showing the log of the wet weight (grams) of biomass among 
# the levels of paint type against days1, days2 and midTrips.

# These smooth plots show that there is a variation of log of the wet weight
# among the levels of paint type factor.

require(nlme)
require(lattice) 
library(ggplot2)
qplot(x=days2, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
qplot(x=days1, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
qplot(x=midTrips, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)





require(nlme)
require(lattice) 
library(ggplot2)

qplot(x=days2, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)

qplot(x=days1, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)


qplot(x=midTrips, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)


# The smooth plots also show there is a variation of log of the wet weight among different testLoc.




################### Scatter plot matrix Figure7 ####################

# Sactterplot matrix of the predictor variables

# The scatter plot matrix contains all the pairwise scatter plots of the variables on a single page
# in a matrix format. Here, there are 4 variables (days1, days2, midTrips and paitType), the scatter 
# plot matrix will have 4 rows and 4 columns and the 2nd row and 3rd column of this matrix is a plot 
# of days2 and midTrips.

# The sactter plot matrix shows that there seems to be some outliers in the data. However, it does
# not seem any paiwise correlation between the vaiables.Thus for further analysis, we will keep all
# four variables.


boat<-unique(samplesdata[, c("days1", "days2", "midTrips", "paintType")])
pairs(unique(boat[, c("days1", "days2", "midTrips", "paintType")]))

############## GLMM ##################

# For GLMM, we consider data from LocID "IB"  "HA"  "PJ"  "PB"  "HD"  "HAH" "HP"  "HH"  "HK".
# Thus, here the dataset is samplesdata_testLoc. 


######## Selecting a correct model ########################

# Generalized linear mixed-effect model (GLMM) considering all predictors (and their interactions)
# and boatID as a random effect.



library(lme4)

samples.logit_all2wayglmer1<-glmer(wetWeight1~days1+days2+midTrips+paintType+factor(LocID)+days1:days2+
	days1:midTrips+ days1:paintType+days2:midTrips+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)
			
summary(samples.logit_all2wayglmer1)



## Comments: 

# The summary of samples.logit_all2wayglmer1 suggests dropping days2:midTrips. 



samples.logit_all2wayglmer2<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips)
			
summary(samples.logit_all2wayglmer2)



# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer1,samples.logit_all2wayglmer2, test="F")


# Which suggests to accept the reduced model (p-value 0.5496). In the reduced model
# AIC is less than the full model. Thus in terms of 
# AIC reduced model (samples.logit_all2wayglmer2) is the better model.However, we perform a
# parametric bootstrap test to take a final decision.  


##############  Parametric Bootstrap test #############

pboot<-function(m0,m1){
s<-simulate(m0)
L0<-logLik(refit(m0,s))
L1<-logLik(refit(m1,s))
2*(L1-L0)
}

samples.logit_all2wayglmer.obsdev1.2<-c(2*(logLik(samples.logit_all2wayglmer1)-logLik(samples.logit_all2wayglmer2)))


set.seed(1001)

samples.logit_all2wayglmer.PB1.2<-replicate(500,pboot(samples.logit_all2wayglmer2, samples.logit_all2wayglmer1))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB1.2>samples.logit_all2wayglmer.obsdev1.2)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)



## PB shows a p-value 0.57, indicating reduced model (samples.logit_all2wayglmer2) is beter model. 






################ Residual Plot  ####################


## Full Model (samples.logit_all2wayglmer1)


scatter.smooth(fitted(samples.logit_all2wayglmer1), residuals(samples.logit_all2wayglmer1, type="pearson"), 
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")


## Reduced Model (samples.logit_all2wayglmer2)


scatter.smooth(fitted(samples.logit_all2wayglmer2), residuals(samples.logit_all2wayglmer2, type="pearson"), 
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")





# The residual shows that there are some outliers. There is a curved pattern that indicates 
# the residuals are non-linear. 


########## Q-Q plot for Residuals ###########


## Full Model (samples.logit_all2wayglmer1)

qqnorm(residuals(samples.logit_all2wayglmer1, type="pearson"), main="Q-Q Normal-Residuals")

qqline(residuals(samples.logit_all2wayglmer1, type="pearson"), col="red")



## Reduced Model (samples.logit_all2wayglmer2)

qqnorm(residuals(samples.logit_all2wayglmer2, type="pearson"), main="Q-Q Normal-Residuals")

qqline(residuals(samples.logit_all2wayglmer2, type="pearson"), col="red")

# q-q plot for residual shows a problem with the Gaussian assumption for both the model.


###########  Q-Q Plot for Random effect ######

## Full Model (samples.logit_all2wayglmer1)


samples.logit_all2wayglmer1r<-lmer(wetWeight1~days1+days2+midTrips+paintType+factor(LocID)+days1:days2+
	days1:midTrips+ days1:paintType+days2:midTrips+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)

			
summary(samples.logit_all2wayglmer1r)


fit_samples.logit_all2wayglmer1r<-fitted(samples.logit_all2wayglmer1r)

#library(glmfun)

res_samples.logit_all2wayglmer1r<-residuals(samples.logit_all2wayglmer1r)



raf<-ranef(samples.logit_all2wayglmer1r)

rand.samples.logit_all2wayglmer1r<-raf[[1]][, "(Intercept)"]

qqnorm(rand.samples.logit_all2wayglmer1r, main="Q-Q plot for random intercept")
qqline(rand.samples.logit_all2wayglmer1r)



## Reduced Model (samples.logit_all2wayglmer2)


samples.logit_all2wayglmer2r<-lmer(wetWeight1~days1+days2+midTrips+paintType+factor(LocID)+days1:days2+
	days1:midTrips+ days1:paintType+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)

			
summary(samples.logit_all2wayglmer2r)


fit_samples.logit_all2wayglmer2r<-fitted(samples.logit_all2wayglmer2r)

#library(glmfun)

res_samples.logit_all2wayglmer2r<-residuals(samples.logit_all2wayglmer2r)



raf<-ranef(samples.logit_all2wayglmer2r)

rand.samples.logit_all2wayglmer2r<-raf[[1]][, "(Intercept)"]

qqnorm(rand.samples.logit_all2wayglmer2r, main="Q-Q plot for random intercept")
qqline(rand.samples.logit_all2wayglmer2r)






# Q-Q plot for random intercept both for reduced model looks good.  



## More comments: From Andrew








# Dropping midTrips:paintType


samples.logit_all2wayglmer3<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-midTrips:paintType)


			
summary(samples.logit_all2wayglmer3)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer3,samples.logit_all2wayglmer2, test="F")



# Which suggests no reason not to accept the reduced model (p-value=0.1852). AIC of the reduced model is slightly less in 
# the reduced  than the full model. We perform a parametric bootstrap test to take a final decision.  


##############  Parametric Bootstrap test #############


samples.logit_all2wayglmer.obsdev2.3<-c(2*(logLik(samples.logit_all2wayglmer2)-logLik(samples.logit_all2wayglmer3)))


set.seed(1001)

samples.logit_all2wayglmer.PB2.3<-replicate(500,pboot(samples.logit_all2wayglmer3, samples.logit_all2wayglmer2))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB2.3>samples.logit_all2wayglmer.obsdev2.3)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.3, df=2, lower.tail=FALSE)



## PB shows a p-value 0.038, indicating the model (samples.logit_all2wayglmer2) is better model.As PB is the 
## most powerful test, we keep midTrips:paintType and the model is samples.logit_all2wayglmer2.  









# Dropping LocID


samples.logit_all2wayglmer4<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-factor(LocID))



summary(samples.logit_all2wayglmer4)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer2,samples.logit_all2wayglmer4, test="F")



# Which suggests not to accept the reduced model.AIC is also greater for reduced model (551 compared to 583). We perform a
# parametric bootstrap test to take a final decision.  



samples.logit_all2wayglmer.obsdev2.4<-c(2*(logLik(samples.logit_all2wayglmer2)-logLik(samples.logit_all2wayglmer4)))


set.seed(1001)

samples.logit_all2wayglmer.PB2.4<-replicate(500,pboot(samples.logit_all2wayglmer4, samples.logit_all2wayglmer2))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB2.4>samples.logit_all2wayglmer.obsdev2.4)

## Classical test: pchisq(samples.logit_all2wayglmer.obsdev3.4, df=1, lower.tail=FALSE)



## PB shows a p-value 0, indicating the model (samples.logit_all2wayglmer2) is better model. 







## Dropping midTrips




samples.logit_all2wayglmer5<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-midTrips)



			
summary(samples.logit_all2wayglmer5)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer2,samples.logit_all2wayglmer5, test="F")



# Which suggests to acept the reduced model (p-value=.155). We perform a
# parametric bootstrap test to take a final decision.  



samples.logit_all2wayglmer.obsdev2.5<-c(2*(logLik(samples.logit_all2wayglmer2)-logLik(samples.logit_all2wayglmer5)))


set.seed(1001)

samples.logit_all2wayglmer.PB2.5<-replicate(500,pboot(samples.logit_all2wayglmer5, samples.logit_all2wayglmer2))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB2.5>samples.logit_all2wayglmer.obsdev2.5)

## Classical test: pchisq(samples.logit_all2wayglmer.obsdev4.5, df=1, lower.tail=FALSE)




## PB shows a p-value 0.098, indicating the reduced model (samples.logit_all2wayglmer5) is better model. 



samples.logit_all2wayglmer5r<-lmer(wetWeight1~days1+days2+paintType+factor(LocID)+days1:days2+
	days1:midTrips+ days1:paintType+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)



summary(samples.logit_all2wayglmer5r)









## Dropping days1:days2




samples.logit_all2wayglmer6<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-midTrips-days1:days2)



			
summary(samples.logit_all2wayglmer6)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer5,samples.logit_all2wayglmer6, test="F")



# Which suggests to accept the reduced model (p-value=.2019). We perform a
# parametric bootstrap test to take a final decision.  



samples.logit_all2wayglmer.obsdev5.6<-c(2*(logLik(samples.logit_all2wayglmer5)-logLik(samples.logit_all2wayglmer6)))


set.seed(1001)

samples.logit_all2wayglmer.PB5.6<-replicate(500,pboot(samples.logit_all2wayglmer6, samples.logit_all2wayglmer5))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB5.6>samples.logit_all2wayglmer.obsdev5.6)

## Classical test: pchisq(samples.logit_all2wayglmer.obsdev4.5, df=1, lower.tail=FALSE)




## PB shows a p-value 0.054, indicating the the reduced model (samples.logit_all2wayglmer6) is better model. AIC
## for both models are almost same. 




samples.logit_all2wayglmer6r<-lmer(wetWeight1~days1+days2+paintType+factor(LocID)+
	days1:midTrips+ days1:paintType+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)



summary(samples.logit_all2wayglmer6r)






## Dropping days1:paintType2




samples.logit_all2wayglmer7<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-midTrips-days1:days2-days1:paintType)



			
summary(samples.logit_all2wayglmer7)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer6,samples.logit_all2wayglmer7, test="F")



# Which suggests not to accept the reduced model (p-value=.00118). We perform a
# parametric bootstrap test to take a final decision.  



samples.logit_all2wayglmer.obsdev6.7<-c(2*(logLik(samples.logit_all2wayglmer6)-logLik(samples.logit_all2wayglmer7)))


set.seed(1001)

samples.logit_all2wayglmer.PB6.7<-replicate(500,pboot(samples.logit_all2wayglmer7, samples.logit_all2wayglmer6))


## parametric bootstrap-based test (PB))

mean(samples.logit_all2wayglmer.PB6.7>samples.logit_all2wayglmer.obsdev6.7)




## PB shows a p-value 0.274, indicating the the reduced model (samples.logit_all2wayglmer7) is better model.
# Howevre, AIC for model samples.logit_all2wayglmer6 is less (551 compared to 560) than samples.logit_all2wayglmer7.
# Thus we consider samples.logit_all2wayglmer6 is better model.





## Dropping days2:paintType2




samples.logit_all2wayglmer8<-update(samples.logit_all2wayglmer1,.~.-days2:midTrips-midTrips-days1:days2-days2:paintType)



			
summary(samples.logit_all2wayglmer8)




# The reduced model can be tested against the full model by using following anova


anova(samples.logit_all2wayglmer6,samples.logit_all2wayglmer8, test="F")



# Which suggests not to accept the reduced model (p-value=.001999). Also AIC for samples.logit_all2wayglmer6
# is still less than samples.logit_all2wayglmer8 (551 compared to 559).Thus we consider samples.logit_all2wayglmer6 is better model.


## Thus for GLMM, samples.logit_all2wayglmer6 is the best model which consists of the variabls: days1, days2,
## paintType, LocID, days1:midTrps, days1:paintType, days2:paintType and midtrips:paintType.






## No more dropping. The model samples.logit_all2wayglmer6 is the final model.



# So the model samples.logit_all2wayglmer6
# is probably the best model to use, and this model suggest that days1, days2 and paintType is statistically significant 
# for fouling biomass. LocID is also significent. 


#### The summary of full model suggests dropping days2:midTrips. Refitting without
#### this term, and then examining the significance of terms in the resulting model, suggests dropping midTrips.
#### Continuing in the same way, we conclude that samples.logit_all2wayglmer6 is probably the best
#### model to use. 




################ Residual Plot  ####################



## Reduced Model (samples.logit_all2wayglmer6)


scatter.smooth(fitted(samples.logit_all2wayglmer6), residuals(samples.logit_all2wayglmer6, type="pearson"), 
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")





# The residual shows that there are some outliers. There is a curved pattern that indicates 
# the residuals are non-linear. It is likely that plot of of this model seems to be better than the full model.   


########## Q-Q plot for Residuals ###########




## Reduced Model (samples.logit_all2wayglmer6)

qqnorm(residuals(samples.logit_all2wayglmer6, type="pearson"), main="Q-Q Normal-Residuals")

qqline(residuals(samples.logit_all2wayglmer6, type="pearson"), col="red")

# q-q plot for residual shows with the Gaussian assumption. It seems that the samples.logit_all2wayglmer6 is good.
# It shows that Gaussian assumption is more appropriate for samples.logit_all2wayglmer6 than the full model.

 


###########  Q-Q Plot for Random effect ######




## Reduced Model (samples.logit_all2wayglmer6)


samples.logit_all2wayglmer6r<-lmer(wetWeight1~days1+days2+paintType+factor(LocID)+
	days1:midTrips+ days1:paintType+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc, na.action=na.exclude)



summary(samples.logit_all2wayglmer6r)




fit_samples.logit_all2wayglmer6r<-fitted(samples.logit_all2wayglmer6r)

#library(glmfun)

res_samples.logit_all2wayglmer6r<-residuals(samples.logit_all2wayglmer6r)



raf<-ranef(samples.logit_all2wayglmer6r)

rand.samples.logit_all2wayglmer6r<-raf[[1]][, "(Intercept)"]

qqnorm(rand.samples.logit_all2wayglmer6r, main="Q-Q plot for random intercept")
qqline(rand.samples.logit_all2wayglmer6r)





# Q-Q plot for random intercept both for the reduced model looks good (better than the full model).  







###########  GAMM ################

# There may be a case that the relationship between the variables is expected to be of a complex form and
# not easily fitted by standard linear or non-linear models. GAM do not involve strong assumptions about the
# relationship that is implicit in standard parametic regression. Infact, such assumptions may force the fitted
# relationship away from its natural path at critical points. 

# Here we use a GAMM model considering b-splines of days1 and days2.Here we use thin plate regression and cubic splines.


######### Here we use the model that we select possibly the best model (samples.logit_all2wayglmer4). We check necessary
######### for the model and finally draw conclusion.


library(mgcv)


samplesdata_testLoc.gam1<-gam(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=samplesdata_testLoc)



gam.check(samplesdata_testLoc.gam1)


# gam.check of the model samplesdata_testLoc.gam1 shows that QQ plot is very close to a straight line, suggesting the distributional
# assumption is reasonable.The residual plot suggests that variance is approximately constant (Not sure!). The Histogram of 
# residual is approximately normal. The response against fitted values may be problematic!!



samplesdata_testLoc.gamm1<-gamm(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=samplesdata_testLoc, random=list(boatID=~1))




summary(samplesdata_testLoc.gamm1$gam)

anova(samplesdata_testLoc.gamm1$gam)


plot(samplesdata_testLoc.gamm1$lme) 

# Problematic (not sure!)

#################### selecting a best GAMM model ####################


# drop midTrips:paintType

samplesdata_testLoc.gam2<-gam(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType,
					data=samplesdata_testLoc)

gam.check(samplesdata_testLoc.gam2)

# as before

samplesdata_testLoc.gamm2<-gamm(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType, 
					data=samplesdata_testLoc, random=list(boatID=~1))

summary(samplesdata_testLoc.gamm2$gam)

anova(samplesdata_testLoc.gamm2$gam)

plot(samplesdata_testLoc.gamm2$lme) 

# As before


# The anova(samplesdata_testLoc.gamm1$lme, samplesdata_testLoc.gamm2$lme) (p-value=0.482) and the AIC suggests (548 compared to 551)
# that the reduced model samplesdata_testLoc.gamm2 is the better than the previous model.

anova(samplesdata_testLoc.gamm2$gam)

# Drop days1:paintType


samplesdata_testLoc.gam3<-gam(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:midTrips+days2:paintType,
					data=samplesdata_testLoc)

gam.check(samplesdata_testLoc.gam3)

# as before


samplesdata_testLoc.gamm3<-gamm(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:midTrips+days2:paintType, 
					data=samplesdata_testLoc, random=list(boatID=~1))

summary(samplesdata_testLoc.gamm3$gam)

anova(samplesdata_testLoc.gamm3$gam)

plot(samplesdata_testLoc.gamm3$lme) 



# As before


# The anova(samplesdata_testLoc.gamm2$lme, samplesdata_testLoc.gamm3$lme) (p-value=0.4679) and the AIC (546 compared to 578)
# suggests that the reduced model samplesdata_testLoc.gamm3 is better than the previous model.

anova(samplesdata_testLoc.gamm3$gam)


# Drop days2:midTrips




samplesdata_testLoc.gam4<-gam(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:paintType,
					data=samplesdata_testLoc)

gam.check(samplesdata_testLoc.gam4)

# as before


samplesdata_testLoc.gamm4<-gamm(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:paintType, 
					data=samplesdata_testLoc, random=list(boatID=~1))

summary(samplesdata_testLoc.gamm4$gam)

anova(samplesdata_testLoc.gamm4$gam)

plot(samplesdata_testLoc.gamm4$lme) 



# As before


# The anova(samplesdata_testLoc.gamm3$lme, samplesdata_testLoc.gamm4$lme) (p-value=0.3086) and the AIC (545 compared to 546)
# suggests that the reduced model samplesdata_testLoc.gamm4 is better than the previous model.

anova(samplesdata_testLoc.gamm4$gam)

# Drop midTrips




samplesdata_testLoc.gam5<-gam(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:paintType,
					data=samplesdata_testLoc)

gam.check(samplesdata_testLoc.gam5)

# as before


samplesdata_testLoc.gamm5<-gamm(wetWeight1~s(days1, bs="tp")+s(days2, bs="cr")+paintType+factor(LocID)+days1:days2+
				days1:midTrips+days2:paintType, 
					data=samplesdata_testLoc, random=list(boatID=~1))

summary(samplesdata_testLoc.gamm5$gam)

anova(samplesdata_testLoc.gamm5$gam)

plot(samplesdata_testLoc.gamm5$lme) 



# As before


# The anova(samplesdata_testLoc.gamm4$lme, samplesdata_testLoc.gamm5$lme) (p-value=0.7841) and the AIC (543 compared to 545)
# suggests that the reduced model samplesdata_testLoc.gamm5 is better than the previous model.

anova(samplesdata_testLoc.gamm5$gam)





# No more dropping. samplesdata_testLoc.gamm5 is possibly the best model to be used.The variables included in the model are
# days1, days2, paintType, LocID, days1:days2, days1:midTrips and days2:paintType

## All the variables included in the model are significant.












################# Figure 9-11 #############

# relationship between the proportions of zero samples in each location and the three potential explanatory variables, 
# days1, days2 and midTrips. 

## data set  propzerodata_testLoc

# These figures suggest that days1 and midTrips have very little predictive power in this context. In this 
# initial inspection the predictive power of days2 appeared to be more promising.   


newy1<-plot(propzerodata_testLoc$days1, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Number of days since the vessel was last used", 
	ylab="Proportion of zero samples", main="")
fit1<-lsfit(propzerodata_testLoc$days1, propzerodata_testLoc$n0/propzerodata_testLoc$n)
abline(fit1)



newy2<-plot(propzerodata_testLoc$days2, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Number of days since the vessel was last cleaned ot antifouled", 
	ylab="Proportion of zero samples", main="")

fit2<-lsfit(propzerodata_testLoc$days2, propzerodata_testLoc$n0/propzerodata_testLoc$n)
abline(fit2)


newy3<-plot(propzerodata_testLoc$midTrips, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Median number of trips per annum", 
	ylab="Proportion of zero samples", main="")
fit3<-lsfit(propzerodata_testLoc$midTrips, propzerodata_testLoc$n0/propzerodata_testLoc$n)
abline(fit3)

############# Quasi-likelihood Analysis ################

# Quasi likelihood analysis to test the predictive power of days1, days2 and midTrips.


glm_quasi<-glm(n0/n~days1+days2+midTrips, data=propzerodata_testLoc, family=quasibinomial, method="glm.fit", weights = n, trace=F) 
summary(glm_quasi)

# Our initial deduction is also supported by the quasi likelihood analysis that the predictive power of days2 is more convincing as
# days2 is statistically significant for the proportion of non-zero samples. 

############# Figure 12 #################

# Actual and predicted proportion of zero samples based on quasi-likelihood analysis.

glm_quasi_days2<-glm(n0/n~days2, data=propzerodata_testLoc, family=quasibinomial, weights=n, method="glm.fit", trace=F)

#bitmap("Zero_samples_GLM.tif", type="tiffpack", width=8, height=8, res=600)
newx<-seq(0, 1000, len=500)
newy<-predict(glm_quasi_days2, newdata=data.frame(days2=newx), type="link", se.fit=T, conf.level=0.95)
plot(propzerodata_testLoc$days2, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Days since vessel last cleaned", ylab="Proportion of clean samples")
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]]+2*newy$se.fit), lty=2)
lines(newx, quasibinomial()$linkinv(newy[[1]]-2*newy$se.fit), lty=2)

# The graph shows that the proportion of clean samples diminishes very quickly after vessel has been clean or antifouled. 



########################### Linear Mixed effect Model (LMM)##################

# Nonzero wet weight taken only from HA and removing outliers

# Dataset nonzerosamples1_HA


################################### Fitting of MRM, this not a part of Analysis

## lww.lm1_1<-lm(logwetWeight~days1+days2+midTrips+factor(paintType)+days1:days2+days1:midTrips+days1:factor(paintType)+
##	days2:midTrips+days2:factor(paintType)+midTrips:factor(paintType), data=nonzerosamples1)
## summary(lww.lm1_1)

# Days2 and midTrips are statistically significant for fouling biomass. PaintType itself is not a significant but it’s 
# interaction with days2 and midTrips are statistically significant. If the vessel is not cleaned one day, then the log of the wet weight
# is increased by .004 grams given an ablative antifouling paint is used and other variables are constant. On the other hand, similar
# if ablative antifouling paint is used the log of the wet weight is decreased by .036 grams, if median number of trips is 
# increased by one day when other variables are constant. 
# paintType influence the effect of vessel activity on the log of the wet weight of biofouling. On an avaege days2, self-polishing
# paint is significant on fouling biomass campared to ablative paint.
# On the other hand, on an average midTrips, hard paint and self-polishing both are significant for fouling biomass compared to 
# ablative paint.  


# anova(lww.lm1_1)


################### Residual plot (without outliers) ######

# scatter.smooth(fitted(lww.lm1_1), residuals(lww.lm1_1, type="pearson"), xlab="Fitted Values", ylab="Residuals", main="Residuals Plot")
# abline(h=0, col="red")


# The graph shows that constant variance in the vertical. There does not seem to be heteroscedascity, non-linearity and outliers. 
# Thus no change in the model are nesessary and the graphic supports the use of MRM. 


####################################### LMM ################

# LMM is used as considering boatID as a random effect. The analysis is done by removing the outliers.
# We first fit a LMM considering all fixed effects and their interaction. However, we do not consider
# the interaction of LocID with other fixed efects as LocID has nine levels: "IB"  "HA"  "PJ"  "PB"  "HD"  "HAH" "HP"  "HH"  "HK".

 

library(lme4)
require(nlme)
nonzerosamples1_HA$boatID_id<-factor(nonzerosamples1_HA$boatID)

nonzerosamples1_HA_new<-groupedData(wetWeight~days2|boatID_id, data=nonzerosamples1_HA)

lww_HA.lmm1<-lme(log(wetWeight)~days1+days2+midTrips+paintType+days1:days2+days1:midTrips+days1:paintType+
	days2:midTrips+days2:paintType+midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm1)

anova(lww_HA.lmm1)




# ANOVA suggests that days2 is statistically significant. However it also shows that midTrips is significant at 5% level
# of significance. 





############### Residual plot, for the model (lww.lmm1) ##########

# plot(lww_HA.lmm1)	# More convenient.


scatter.smooth(fitted(lww_HA.lmm1), residuals(lww_HA.lmm1, type="pearson"), xlab="Fitted Values", 
	ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")



# variance-adjusted, outermost residuals have approximately constant variance as the residuals
# do not have any relationship. Theare is no curved patern that indicates that the residual are
# normlly distributed. There are no substantial outliers. Thus it can be said that LMM model fits the data
# well. 


############ Q-Q plot Figure 14 ##############

## Q-Q plot for Random effect

ref.boat<-ranef(lww_HA.lmm1)[[1]]
ref.var.boat<-tapply(residuals(lww_HA.lmm1,type="pearson", level=1), nonzerosamples1_HA_new$boatID_id, var)
qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
qqline(ref.boat, col="red")


## Q-Q plot for Residuals ################

qqnorm(residuals(lww_HA.lmm1, type="pearson"), main="Q-Q Normal-Residuals")
qqline(residuals(lww_HA.lmm1, type="pearson"), col="red")





# The Q-Q plot for residuals follows approximately in the line, that means two distributions agree 
# after linearly transforming the values. However, Q-Q plot for for random effect slightly differ 
# from straight line.  





######## Selecting the best model by using LMM ###########

# lww_HA.lmm1 suggests dropping paintType. We refitt the model without this term and compare with the full model by using AIC. 


lww_HA.lmm2<-lme(log(wetWeight)~days1+days2+midTrips+days1:days2+days1:midTrips+days1:paintType+
	days2:midTrips+days2:paintType+midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm2)

anova(lww_HA.lmm2)



# The anova(lww_HA.lmm1, lww_HA.lmm2) test is problematic here (p-value 0.1232, suggesting to accept reduced model model). 
# However, AIC(lww_HA.lmm1, lww_HA.lmm2) is also the same (539) for both models. We perform a bootstrap test for final decision. 
 



# Parametric Bootstrap test (PB) ####


lww_HA.lmm1r<-lmer(log(wetWeight)~days1+days2+midTrips+paintType+days1:days2+days1:midTrips+days1:paintType+
	days2:midTrips+days2:paintType+midTrips:paintType+ (1|boatID), 
		data=nonzerosamples1_HA, na.action=na.exclude)

lww_HA.lmm2r<-update(lww_HA.lmm1r,.~.-paintType)



lww_HA.lmm1.2.obsdev<-c(2*(logLik(lww_HA.lmm1r)-logLik(lww_HA.lmm2r)))


set.seed(1001)

lww_HA.lmm1.2.PB1<-replicate(500,pboot(lww_HA.lmm2r, lww_HA.lmm1r))

mean(lww_HA.lmm1.2.PB1> lww_HA.lmm1.2.obsdev)




### PB shows a p-value 0.62, indicating the model (lww_HA.lmm2r) is better model. 

# Thus we keep the model lww_HA.lmm2 and try testing dropping day2:paintType



lww_HA.lmm3<-lme(log(wetWeight)~days1+days2+midTrips+days1:days2+days1:midTrips+days1:paintType+
	days2:midTrips+midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm3)

anova(lww_HA.lmm3)



# The anova(lww_HA.lmm2, lww_HA.lmm3) test is problematic here (p-value 0.0004, suggesting to accept lww_HA.lmm2). 
# But AIC for lww_HA.lmm3 is less than lww_HA.lmm2) (520 compared to 539). We perform a bootstrap test for final decision. 
 



# Parametric Bootstrap test (PB) ####



lww_HA.lmm3r<-update(lww_HA.lmm1r,.~.-paintType-days2:paintType)



lww_HA.lmm2.3.obsdev<-c(2*(logLik(lww_HA.lmm2r)-logLik(lww_HA.lmm3r)))


set.seed(1001)

lww_HA.lmm2.3.PB1<-replicate(500,pboot(lww_HA.lmm3r, lww_HA.lmm2r))

mean(lww_HA.lmm2.3.PB1> lww_HA.lmm2.3.obsdev)




### PB shows a p-value 0.132, indicating the model (lww_HA.lmm3r) is better model. 

# Thus we keep the model lww_HA.lmm3 and try testing dropping day2:midTrips



lww_HA.lmm4<-lme(log(wetWeight)~days1+days2+midTrips+days1:days2+days1:midTrips+days1:paintType+
	midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm4)

anova(lww_HA.lmm4)



# The anova(lww_HA.lmm3, lww_HA.lmm4) test is problematic here (p-value 0.0001, suggesting to accept lww_HA.lmm3). 
# But AIC for lww_HA.lmm4 is less than lww_HA.lmm3) (497 compared to 520). We perform a bootstrap test for final decision. 
 



# Parametric Bootstrap test (PB) ####



lww_HA.lmm4r<-update(lww_HA.lmm1r,.~.-paintType-days2:paintType-days2:midTrips)



lww_HA.lmm3.4.obsdev<-c(2*(logLik(lww_HA.lmm3r)-logLik(lww_HA.lmm4r)))


set.seed(1001)

lww_HA.lmm3.4.PB1<-replicate(500,pboot(lww_HA.lmm4r, lww_HA.lmm3r))

mean(lww_HA.lmm3.4.PB1> lww_HA.lmm3.4.obsdev)




### PB shows a p-value 0.744, indicating the model (lww_HA.lmm4r) is better model. 

# Thus we keep the model lww_HA.lmm4 and try testing dropping days1:day2



lww_HA.lmm5<-lme(log(wetWeight)~days1+days2+midTrips+days1:midTrips+days1:paintType+
	midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm5)

anova(lww_HA.lmm5)



# The anova(lww_HA.lmm4, lww_HA.lmm5) test is problematic here (p-value 0.0001, suggesting to accept lww_HA.lmm4). 
# But AIC for lww_HA.lmm5 is less than lww_HA.lmm4) (475 compared to 497). We perform a bootstrap test for final decision. 
 



# Parametric Bootstrap test (PB) ####



lww_HA.lmm5r<-update(lww_HA.lmm1r,.~.-paintType-days2:paintType-days2:midTrips-days1:days2)



lww_HA.lmm4.5.obsdev<-c(2*(logLik(lww_HA.lmm4r)-logLik(lww_HA.lmm5r)))


set.seed(1001)

lww_HA.lmm4.5.PB1<-replicate(500,pboot(lww_HA.lmm5r, lww_HA.lmm4r))

mean(lww_HA.lmm4.5.PB1> lww_HA.lmm4.5.obsdev)




### PB shows a p-value 0.642, indicating the model (lww_HA.lmm5r) is better model. 

# Thus we keep the model lww_HA.lmm5 and try testing dropping days1:paintType 




lww_HA.lmm6<-lme(log(wetWeight)~days1+days2+midTrips+days1:midTrips+
	midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm6)

anova(lww_HA.lmm6)



# The anova(lww_HA.lmm5, lww_HA.lmm6) test is problematic here (p-value 0.0005, suggesting to accept lww_HA.lmm5). 
# But AIC for lww_HA.lmm6 is less than lww_HA.lmm5) (455 compared to 475). We perform a bootstrap test for final decision. 
 



# Parametric Bootstrap test (PB) ####



lww_HA.lmm6r<-update(lww_HA.lmm1r,.~.-paintType-days2:paintType-days2:midTrips-days1:days2-days1:paintType)



lww_HA.lmm5.6.obsdev<-c(2*(logLik(lww_HA.lmm5r)-logLik(lww_HA.lmm6r)))


set.seed(1001)

lww_HA.lmm5.6.PB1<-replicate(500,pboot(lww_HA.lmm6r, lww_HA.lmm5r))

mean(lww_HA.lmm5.6.PB1> lww_HA.lmm5.6.obsdev)




### PB shows a p-value 0.124, indicating the model (lww_HA.lmm6r) is better model. 

# Thus we keep the model lww_HA.lmm6 and try testing dropping midTrips:paintType 




lww_HA.lmm7<-lme(log(wetWeight)~days1+days2+midTrips+days1:midTrips, random=~1|boatID_id,
			data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm7)

anova(lww_HA.lmm7)



# The anova(lww_HA.lmm6, lww_HA.lmm7) test is problematic here. 
# But AIC for lww_HA.lmm7 is greater than lww_HA.lmm6) (512 compared to 455). 
# Thus we keep the term midTrips:paintType and consider lww_HA.lmm6 is the better model.


# Consider testing dropping days1   
 



lww_HA.lmm8<-lme(log(wetWeight)~days2+midTrips+days1:midTrips+
	midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm8)

anova(lww_HA.lmm8)



# The anova(lww_HA.lmm8, lww_HA.lmm6) test is problematic here (p-value 0.0046, suggesting to accept lww_HA.lmm6). 
# But AIC for lww_HA.lmm8 is less than lww_HA.lmm6) (445 compared to 455). We perform a bootstrap test for final decision. 
 






# Parametric Bootstrap test (PB) ####


lww_HA.lmm8r<-update(lww_HA.lmm1r,.~.-paintType-days2:paintType-days2:midTrips-days1:days2-days1:paintType-days1)




lww_HA.lmm6.8.obsdev<-c(2*(logLik(lww_HA.lmm6r)-logLik(lww_HA.lmm8r)))


set.seed(1001)

lww_HA.lmm6.8.PB1<-replicate(500,pboot(lww_HA.lmm8r, lww_HA.lmm6r))

mean(lww_HA.lmm6.8.PB1> lww_HA.lmm6.8.obsdev)




### PB shows a p-value 0.114, indicating the model (lww_HA.lmm8r) is better model. 

# Thus we keep the model lww_HA.lmm8 and try testing dropping days1:midTrips 







lww_HA.lmm9<-lme(log(wetWeight)~days2+midTrips+
	midTrips:paintType, random=~1|boatID_id, 
		data=nonzerosamples1_HA_new, na.action=na.exclude)

summary(lww_HA.lmm9)

anova(lww_HA.lmm9)



# The anova(lww_HA.lmm8, lww_HA.lmm9) test is problematic here as both model do not consider the same number 
# of observations. But AIC for lww_HA.lmm9 is less than that of lww_HA.lmm8 (440 compared to 445).We can not perform 
# PB test here as of the reason described above. Considering the overall situation, lww_HA.lmm8 would be the best model to be used. 



# Thus lww_HA.lmm8 would be the best model to be used.The variables included in the model are
# days2, midTrips, days1:midTrips and midTrips:paintType. 





# Days2 is strongly significant, whereas midTrips are significant at 5% level of significance. PaintType itself 
# is not a significant but it’s interaction with midTrips are statistically significant. If the vessel is not 
# cleaned one day, then the log of the wet weight is increased by .004 grams given an ablative antifouling paint is used
# and other variables are constant. On the other hand, similar if ablative antifouling paint is used the log of the wet 
# weight is decreased by .029 grams, if median number of trips is increased by one day when other variables are constant. 
# paintType influence the effect of vessel activity on the log of the wet weight of biofouling. On an avaege midTrips, hard  
# paint is significant on fouling biomass campared to ablative paint. 





############### Residual plot, for the model (lww_HA.lmm8) ##########




scatter.smooth(fitted(lww_HA.lmm8), residuals(lww_HA.lmm8, type="pearson"), xlab="Fitted Values", 
	ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")



# Residual plot looks good. Thus the model is reasonable. 


############ Q-Q plot Figure 14 ##############

## Q-Q plot for Random effect

ref.boat<-ranef(lww_HA.lmm8)[[1]]
ref.var.boat<-tapply(residuals(lww_HA.lmm8,type="pearson", level=1), nonzerosamples1_HA_new$boatID_id, var)
qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
qqline(ref.boat, col="red")


## Q-Q plot for Residuals ################

qqnorm(residuals(lww_HA.lmm8, type="pearson"), main="Q-Q Normal-Residuals")
qqline(residuals(lww_HA.lmm8, type="pearson"), col="red")





# The Q-Q plot for residuals follows approximately in the line, that means two distributions agree 
# after linearly transforming the values. However, Q-Q plot for for random effect slightly differ 
# from straight line but more reasonable than before and thus the model is appropriate.  
















###########  GAMM ################

# We fit a GAMM model considering smooths of days1 and days2, and the parametric term of other variables.


library(mgcv)



b_HA.gam1<-gam(log(wetWeight)~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


## gam.check of model b_HA.gam1 shows that QQ plot is very close to a straight line, suggesting the distributional assumption is reasonable. 
## The residual plot suggests that variance is approximately constant as mean increases. The Histogram of residuals is approximtely
## consistent with normality. The plot of response against fitted values shows a positive relationship with a moderate scatter.
## So the model is reasonable.   
   


b_HA.gamm1<-gamm(log(wetWeight)~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+days1:days2+
				days1:midTrips+days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))


summary(b_HA.gamm1$gam)

anova(b_HA.gamm1$gam)


plot(b_HA.gamm1$lme)

# Residual plot looks fine.


############   Selecting a best GAMM model ####################

# The model b_HA.gamm1 suggests dropping smooth tem of days1:midTrips. 



b_HA.gam2<-gam(log(wetWeight)~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+days1:days2+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


gam.check(b_HA.gam2)    # Ok same as before, the model is reasonable



b_HA.gamm2<-gamm(log(wetWeight)~s(days1, bs="tp")+s(days2, bs="cr")+midTrips+paintType+days1:days2+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))



summary(b_HA.gamm2$gam)



plot(b_HA.gamm2$lme)   ## Ok, residual plots looks good



# The anova(b_HA.gamm1$lme, b_HA.gamm2$lme) (p-value=0.9387) and the AIC (386 compared to 388) suggests that the reduced model 
# b_HA.gamm2 is better than b_HA.gamm1.

anova(b_HA.gamm2$gam)

# Drop days1:days2


b_HA.gam3<-gam(log(wetWeight)~s(days1, bs="ts")+s(days2, bs="cr")+midTrips+paintType+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


gam.check(b_HA.gam3)    # Ok same as before, the model is reasonable



b_HA.gamm3<-gamm(log(wetWeight)~s(days1, bs="ts")+s(days2, bs="cr")+midTrips+paintType+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))

# We used thin plate regression splines with shrinkage (TPRS with shrinkage, bs=ts) as we got an error massges in R
# "Singularity in backsolve at level 0, block 1", when we used bs=tp

summary(b_HA.gamm3$gam)



plot(b_HA.gamm3$lme)   ## Ok, residual plots looks good



# The anova(b_HA.gamm2$lme, b_HA.gamm3$lme) (p-value=0.8486) and the AIC (384 compared to 386) suggests that the reduced model 
# b_HA.gamm3 is better than b_HA.gamm1.

anova(b_HA.gamm3$gam)



# dropping days1


b_HA.gam4<-gam(log(wetWeight)~ s(days2, bs="cr")+midTrips+paintType+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


gam.check(b_HA.gam4)    # Ok same as before, the model is reasonable



b_HA.gamm4<-gamm(log(wetWeight)~s(days2, bs="cr")+midTrips+paintType+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))

summary(b_HA.gamm4$gam)



plot(b_HA.gamm4$lme)   ## Ok, residual plots looks good



# The anova(b_HA.gamm3$lme, b_HA.gamm4$lme) (p-value=0.9998) and the AIC (382 compared to 384) suggests that the reduced model 
# b_HA.gamm4 is better than b_HA.gamm1.

anova(b_HA.gamm4$gam)


# Drop paintType



b_HA.gam5<-gam(log(wetWeight)~ s(days2, bs="cr")+midTrips+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


gam.check(b_HA.gam5)    # Ok same as before, the model is reasonable



b_HA.gamm5<-gamm(log(wetWeight)~s(days2, bs="cr")+midTrips+
				days1:paintType+days2:midTrips+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))

summary(b_HA.gamm5$gam)



plot(b_HA.gamm5$lme)   ## Ok, residual plots looks good



# The anova(b_HA.gamm4$lme, b_HA.gamm5$lme) (p-value=0.3034) and the AIC (380 compared to 382) suggests that the reduced model 
# b_HA.gamm5 is better than b_HA.gamm1.

anova(b_HA.gamm5$gam)



# Dropping days2:midTrips




b_HA.gam6<-gam(log(wetWeight)~ s(days2, bs="cr")+midTrips+
				days1:paintType+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA)


gam.check(b_HA.gam6)    # Ok same as before, the model is reasonable



b_HA.gamm6<-gamm(log(wetWeight)~s(days2, bs="cr")+midTrips+
				days1:paintType+days2:paintType+
					midTrips:paintType,data=nonzerosamples1_HA, random=list(boatID=~1))





# we got an error massges in R "Singularity in backsolve at level 0, block 1", when we used bs=cr. 
# Thus we stop dropping and consider b_HA.gamm5 is the best model. 


# No more dropping. b_HA.gamm5 is the best model to be used. The variables included in the model are
# days2, midTrips, days1:paintType, days2:midTrips, days2:paintType, and midTrips:paintType.



# b_HA.gamm5 may be the final model to be used.The model suggests that smooth term days2 and midTrips are strongly significant.
# However the interaction term midTrips:paintType are also significant wherewas days2:midTrips are significant at 5% level of significance.






























































