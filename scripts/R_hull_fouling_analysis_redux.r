# ==============================================================================
#   DEWHA EV project hull-fouling samples analysis
#   Andrew Robinson, Enamul Kabir, and Keith Hayes
# ==============================================================================
## Clean up workspace
rm(list=ls())

## Load libraries
library(nlme)
library(lattice)
library(ggplot2)
library(lme4)
library(mgcv)

require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(boot)
library(Amelia)


require(Hmisc)

source("functions.r")

# ========================================================================
#                           CODE
# ========================================================================
#
#
# ========================================================================
#                           Read in data
# ========================================================================

samplesdata <- read.csv("../data-raw/samples.csv")
propzerodata <- read.csv("../data-raw/aa.csv")
covars <- read.csv("../data-raw/covars.csv")


## Identify non-zero samples and exclude two outliers boat 24 and 51
nonzerosamples <- samplesdata[samplesdata$wetWeight>.5, ]
nonzerosamples1 <- subset(nonzerosamples, boatID!=24 & boatID!=51)

## Identify broad and specific sample locations
## KRH: but note that the sample units in the testLoc list are very different
## The only samples based on the same sampling unit (quadrats) are:
## Hull quadrats (HA), Rudder surface (PJ), Keel fixed (HP)
Loc <- c("Anchor", "Deck", "Dry Hull", "Fishing gear", "Internal spaces",
       "Propeller", "Rudder", "Wet Hull")
testLoc <- c("IB", "HA", "PJ", "PB", "HD", "HAH", "HP", "HH", "HK")
testLoc1 <- c("HA", "PJ", "HP")

## Identify non-zero samples from wet hull
samplesdata_WH <- samplesdata[samplesdata$Location==Loc[8], ]
nonzerosamples_WH <- samplesdata_WH[samplesdata_WH$wetWeight>.5, ]
nonzerosamples1_WH <- subset(nonzerosamples_WH, boatID!=24 & boatID!=51)

## Identify non-zero samples from HA (Hull quadrats)
samplesdata_HA <- samplesdata[samplesdata$LocID=="HA", ]
nonzerosamples_HA <- samplesdata_HA[samplesdata_HA$wetWeight>.5, ]
nonzerosamples1_HA <- subset(nonzerosamples_HA, boatID!=24 & boatID!=51)

## Identify non-zero samples in test locations
samplesdata_testLoc <- subset(samplesdata,LocID %in% testLoc[1:9])
nonzerosamples_testLoc <- subset(nonzerosamples,LocID %in% testLoc[1:9])

## Proportion of zero data in each of the test locations
propzerodata_testLoc <- subset(propzerodata, LocID %in% testLoc[1:9])

## KRH: same data sets based on standard sample units
samplesdata_testLoc1 <- subset(samplesdata,LocID %in% testLoc1[1:3])
nonzerosamples_testLoc1 <- subset(nonzerosamples,LocID %in% testLoc1[1:3])

# EK: non-zero samples based on standard sampling units and removing outliers
nonzerosamples1_testLoc1 <- subset(nonzerosamples1,LocID %in% testLoc1[1:3])


propzerodata_testLoc1 <- subset(propzerodata, LocID %in% testLoc1[1:3])




# ===========================================================================================
#              Graphical summary and data exploration
# ===========================================================================================

## Three Histograms summarizing the number of days since the vessel's last trip,
## the number of days since it was last cleaned or antifouled, and the median
## number of trips per annum.
# par(mfrow=c(2,2))
# hist(covars$days1, prob=F, xlab="Days since last trip", ylab="Frequency", main="")
# 	text(0,36,"a")
# hist(covars$days2, prob=F, xlab="Days since last cleaned or antifouled", ylab="Frequency", main="")
# 	text(0,24,"b")
# hist(covars$midTrips, prob=F, xlab="Median number of trips per annum", ylab="Frequency", main="")
#	text(0,37,"c")

####### FIGURE 1 ################
## The scatter plot matrix contains all the pairwise scatter plots of the variables on a single page
## in a matrix format. Here, there are 4 variables (days1, days2, midTrips and paintType), the scatter
## plot matrix will have 4 rows and 4 columns and the 2nd row and 3rd column of this matrix is a plot
## of days2 and midTrips.

## The scatter plot matrix shows that there seems to be some outliers in the data. However, it does
## not seem any pairwise correlation between the variables. Thus for further analysis, we will keep all
## four variables.

## KRH: I have added to the scatter plot the histograms on the diagonal and the correlation coefficients
## in the upper triangle, this combines all the information from the previous histograms into a
## single figure

## APR: Beautiful!

boat<-unique(samplesdata[, c("days1", "days2", "midTrips", "paintType")])

pairs(x=unique(boat[, c("days1", "days2", "midTrips", "paintType")]), diag.panel = panel.hist,
      upper.panel = panel.cor, lower.panel=points, pch = 20, ps = 8, cex.labels = 1.3)

######### FIGURE 2 ###############
## Histograms of sample wet weight (grams) for the whole sample, and the log of sample
## wet weight (grams) for the whole sample. Histogram of wet weight shows that samples are
## zero-inflated. However, histogram of log transformed data shows that samples are approximately
## normally distributed.
par(mfrow=c(2,2))
hist(samplesdata$wetWeight, prob=F, xlab="Wet-weight of fouling biomass (grams)", ylab="Frequency", main="")
	text(9000,1490,"a")
hist(log(samplesdata$wetWeight), prob=F, xlab="Log wet-weight of fouling biomass", ylab="Frequency", main="")
	text(9,330,"b")

## Histograms of sample wet weight (grams) for the non-zero (biomass) sample, and the log of
## sample wet weight (grams) for the non-zero (biomass) samples. The wet weight fouling biomass of
## non-zero samples is approximately log normally distributed but the log of these non-zero biomasses
## is normally distributed.
hist(nonzerosamples$wetWeight, prob=F, xlab="Wet-weight of non-zero samples (grams)",
     ylab="Frequency", main="" )
text(9000,740, "c")
hist(log(nonzerosamples$wetWeight), prob=F, xlab="Log wet-weight of non-zero samples",
     ylab="Frequency", main="" )
text(9,200, "d")

###### FIGURES 3-5 #############
## Smooth plots showing the log of the wet weight (grams) of biomass among
## the levels of paint type against days1, days2 and midTrips.

## These smooth plots show that there is a variation of log of the wet-weight
## among the levels of paint-type factor. The smooth plots also show there is a variation
## of log of the wet-weight among different testLoc.

## KRH: the variation between the testLoc is not surprising - they are very different sample
## Have added a figure using only those samples based on quadrats.

# EK: Make sense. But we do not need to consider qplot for nonzerosamples_HA.

# qplot(x=days2, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
# qplot(x=days2, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
# qplot(x=days1, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
# qplot(x=midTrips, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"), data=nonzerosamples)
# qplot(x=days2, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)
# qplot(x=days1, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)
# qplot(x=midTrips, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"), data=nonzerosamples_testLoc)

qplot(x=days2, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)


#qplot(x=days2, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      #data=nonzerosamples_HA)

qplot(x=days1, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)

#qplot(x=days1, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      #data=nonzerosamples_HA)

qplot(x=midTrips, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)

#qplot(x=midTrips, y=log(wetWeight), group=paintType, colour=paintType, geom=c("point", "smooth"),
      #data=nonzerosamples_HA)

qplot(x=days2, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)

qplot(x=days1, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)

qplot(x=midTrips, y=log(wetWeight), group=LocID, colour=factor(LocID), geom=c("point", "smooth"),
      data=nonzerosamples_testLoc1)

## KRH: Box plot to see how the mean and variance of the wetweight varies between locations
## (so do we need a random effect at the location level?). Seems to suggest that we don't.
# qplot(factor(LocID), log(wetWeight), geom=c("boxplot", "jitter"), data = nonzerosamples_testLoc1)
qplot(factor(LocID), log(wetWeight), geom=c("boxplot"), data = nonzerosamples_testLoc1)

## EK: Excellent!

# ===============================================================================================================================
#                                   GLMM for fouled/not fouled response (data: samplesdata_testLoc)
# ===============================================================================================================================
## For GLMM, we consider data from LocID "IB"  "HA"  "PJ"  "PB"  "HD"  "HAH" "HP"  "HH"  "HK".
## Thus, here the dataset is samplesdata_testLoc.

####### Selecting a correct model ###################
## Generalized linear mixed-effect model (GLMM) considering all predictors (and their interactions)
## and boatID as a random effect.
## The summary of samples.logit_all2wayglmer1 suggests dropping days2:midTrips.
## KRH: OK here the use of all the locations for a logistic regression of fouled or not fouled
## (wetWeight1) makes sense


## APR Scale continuous parameters

samplesdata_testLoc.s <- droplevels(samplesdata_testLoc)

#for (columns in c("days1","days2","midTrips"))
 # samplesdata_testLoc.s[,columns] <- scale(samplesdata_testLoc.s[,columns])




#samples.logit_all2wayglmer1<-glmer(wetWeight1~days1+days2+midTrips+paintType+factor(LocID)+days1:days2+
#	days1:midTrips+ days1:paintType+days2:midTrips+days2:paintType+
#		midTrips:paintType+(1|boatID),family=binomial,
#			data=samplesdata_testLoc, na.action=na.exclude)

#summary(samples.logit_all2wayglmer1)


#######################################################################################

par(mfrow = c(3,2))
for (columns in c("days1","days2","midTrips")) {
  plot(density(na.omit(samplesdata_testLoc.s[,columns])))
  samplesdata_testLoc.s[,columns] <-
    scale(log(samplesdata_testLoc.s[,columns] + 0.1))
  plot(density(na.omit(samplesdata_testLoc.s[,columns])))
}


#################################################################################################

# Address missingness - test

lapply(samplesdata_testLoc.s,
       function(x) sum(is.na(x)))

naplot(naclus(samplesdata_testLoc.s))

# Substantial

samplesdata_testLoc.s$LocID <- factor(samplesdata_testLoc.s$LocID)

set.seed(300)

bf.out <-
    amelia(samplesdata_testLoc.s[,c("wetWeight1","days1","days2","midTrips","paintType",
                                    "LocID","boatID")],
           m = 500,   ### 100 is just a guess
           cs = "boatID",
           noms = c("paintType","LocID"))

# http://stackoverflow.com/questions/16571580/multi-level-regression-model-on-multiply-imputed-data-set-in-r-amelia-zelig-l

#################################################################################################


#samples.logit_all2wayglmer1.s<-glmer(wetWeight1~days1+days2+midTrips+paintType+factor(LocID)+days1:days2+
	#days1:midTrips+ days1:paintType+days2:midTrips+days2:paintType+
		#midTrips:paintType+(1|boatID),family=binomial,
			#data=samplesdata_testLoc, na.action=na.exclude)

### bobyqa optimizar --- BASE MODEL

samples.logit_all2wayglmer1.s <-
  glmer(wetWeight1 ~ days1 + days2 + midTrips + paintType + factor(LocID) +
        days1:days2 + days1:midTrips + days1:paintType +
        days2:midTrips + days2:paintType +
        midTrips:paintType + (1|boatID),
        family = binomial,
        data = samplesdata_testLoc.s,
        nAGQ = 10,
        control = glmerControl(optimizer="bobyqa"),
        na.action = na.exclude)

####### Residual Plot ####################
## Full Model (samples.logit_all2wayglmer1.s)
par(mfrow = c(1,2))
scatter.smooth(fitted(samples.logit_all2wayglmer1.s),
               residuals(samples.logit_all2wayglmer1.s, type="pearson"),
               xlab="Fitted Values", ylab="Residuals", main="Residuals Plot GLMER1")
abline(h=0, col="red")

## Reduced Model (samples.logit_all2wayglmer2.s)
scatter.smooth(fitted(samples.logit_all2wayglmer2.s),
               residuals(samples.logit_all2wayglmer2.s, type="pearson"),
               xlab="Fitted Values", ylab="Residuals", main="Residuals Plot GLMER2")
abline(h=0, col="red")

## The residual shows that there are some outliers.
## There is a curved pattern that indicates the residuals are non-linear.

######## Q-Q plot for Residuals #################
## Full Model (samples.logit_all2wayglmer1.s)
par(mfrow = c(1,2))
qqnorm(residuals(samples.logit_all2wayglmer1.s, type="pearson"),
       main="Q-Q Normal-Residuals GLMER1")
qqline(residuals(samples.logit_all2wayglmer1.s, type="pearson"), col="red")

## Reduced Model (samples.logit_all2wayglmer2.s)
qqnorm(residuals(samples.logit_all2wayglmer2.s, type="pearson"),
       main="Q-Q Normal-Residuals GLMER2")
qqline(residuals(samples.logit_all2wayglmer2.s, type="pearson"), col="red")

## Q-Q plot for residual shows a problem with the Gaussian assumption for both models.


raf <- ranef(samples.logit_all2wayglmer1.s)
rand.samples.logit_all2wayglmer1.s <- raf[[1]][, "(Intercept)"]
qqnorm(rand.samples.logit_all2wayglmer1.s,
       main="Q-Q plot, random intercept GLMER1", cex.main = 0.75)
qqline(rand.samples.logit_all2wayglmer1.s)

## Now start model selection.  Only get these diagnostics for the final model.

#####################

## Re-run the analysis but drop days1:paintType

samples.logit_all2wayglmer2.s <- update(samples.logit_all2wayglmer1.s,.~.-days1:paintType)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer1.s, samples.logit_all2wayglmer2.s, test="F")

## Which suggests to accept the reduced model (p-value 0.00249). In the reduced model
## AIC is less than the full model. Thus in terms of
## AIC reduced model (samples.logit_all2wayglmer2.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

##########################################

## Now introduce imputation test instead of PB, which chokes because
## of the missingness.

imputeLRT.c(samples.logit_all2wayglmer2.s,  #### SMALLER MODEL FIRST
            samples.logit_all2wayglmer1.s,
            bf.out)

# This says p = 0.308.  Not a significant difference.

##########################################

## Parametric bootstrap-based test (PB) for comparison --- using reduced data

red.1 <- update(samples.logit_all2wayglmer1.s,
                data = na.omit(samplesdata_testLoc.s[,c("wetWeight1","days1","days2",
                    "midTrips","paintType","LocID","boatID")]))
red.2 <- update(samples.logit_all2wayglmer2.s,
                data = na.omit(samplesdata_testLoc.s[,c("wetWeight1","days1","days2",
                    "midTrips","paintType","LocID","boatID")]))

samples.logit_all2wayglmer.PB1.2.s <- replicate(200, pboot(red.2, red.1))

(samples.logit_all2wayglmer.obsdev1.2.s <- 2*(logLik(red.1)-logLik(red.2)))

mean(samples.logit_all2wayglmer.PB1.2.s > samples.logit_all2wayglmer.obsdev1.2.s)

## This mean is the p-value

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)

## PB shows a p-value 0.838, indicating reduced model (samples.logit_all2wayglmer2.s)
## is better model.







































