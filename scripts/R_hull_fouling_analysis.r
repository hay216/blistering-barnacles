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

require(Hmisc)

# ==============================================================================
#                           FUNCTIONS
# ==============================================================================
## Function to produce histrogram in a pairs plot
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col='blue', cex.main = 0.8,...)
}

## Function to produce kernel density estimate in a pairs plot
panel.density <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  den <- density(x, na.rm = T)
  rug(x, ticksize = 0.06, lwd = 1, col = 'red')
  y <- den$y; y <- y/max(y)
  x <- den$x
  lines(x, y)
}

## Function to produce correlation coefficients in a pairs plot
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1 <- abs(cor(x, y, use ='pairwise.complete.obs', method = 'pearson'))
  r2 <- abs(cor(x, y, use ='pairwise.complete.obs', method = 'kendall'))
  txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
  txt2 <- format(c(r2, 0.123456789), digits=digits)[1]
  txt1 <- paste('Pearson = ', txt1, sep='')
  txt2 <- paste('Kendall = ',txt2, sep='')
  txt <- paste(prefix, txt1, txt2, sep=" ")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = 1)
  text(0.5, 0.6, txt1)
  text(0.5, 0.4, txt2)
}

## Function to calculate the deviance (2*likelihood ratio) of two mixed models
pboot<-function(m0,m1){
  s <- simulate(m0)
  L0<-logLik(refit(m0,s))
  L1<-logLik(refit(m1,s))
  2*(L1-L0)
}

imputeLRT <- function(h0, h1, imputed.data.object) {
    if (!(length(fixef(h0)) <  length(fixef(h1))))
        stop("The first argument needs the smaller model.")
    imputed.data.list <- imputed.data.object$imputations
    m <- length(imputed.data.list)
    h0.models <- lapply(imputed.data.list,
                        function (dataset) {
                            update(h0, data = dataset)
                        })
    h0.dev.funs <- lapply(imputed.data.list,
                          function (dataset) {
                              update(h0, data = dataset,
                                     devFunOnly = TRUE)
                          })
    Q.bar.0 <-
        colMeans(do.call(rbind,
                         lapply(h0.models,
                                function (fit) {
                                    c(getME(fit, "theta"), fixef(fit))
                                })))
    h1.models <- lapply(imputed.data.list,
                        function (dataset) {
                            update(h1, data = dataset)
                        })
    h1.dev.funs <- lapply(imputed.data.list,
                          function (dataset) {
                              update(h1, data = dataset,
                                     devFunOnly = TRUE)
                          })
    Q.bar.1 <-
        colMeans(do.call(rbind,
                         lapply(h1.models,
                                function (fit) {
                                    c(getME(fit, "theta"), fixef(fit))
                                })))
    d.prime.m.bar <- mean(unlist(lapply(1:m,
                                        function(i) {
                                            anova(h0.models[[i]],
                                                  h1.models[[i]])$Chisq[2]
                                        })))
    d.L.bar <- mean(unlist(lapply(1:m,
                                  function(i) {
                                      h0.dev.funs[[i]](Q.bar.0) -
                                          h1.dev.funs[[i]](Q.bar.1)
                                  })))
    p0 <- length(Q.bar.0)
    p1 <- length(Q.bar.1)
    k <- p1 - p0
    rL <- (m + 1) / ((m - 1) * k) * (d.prime.m.bar - d.L.bar) # 3.8
    D.L <- d.L.bar / ((1 + rL) * k)                           # 3.7
    v <- k * (m - 1)
    w <- ifelse(v > 4,                                        # 2.7
                4 + (v - 4) * (1 + (1 - v/2) / rL)^2,
                v / 2 * (1 + 1/k) * (1 + 1/rL)^2)
    Pval <- 1 - pf(D.L, k, w)
                                        #  browser()
    return(list(
        D.L = D.L,
        Pval = Pval,
        rL = rL,
        d.L.bar = d.L.bar,
        d.prime.m.bar = d.prime.m.bar,
#    Q.bar.0 = Q.bar.0,
#    Q.bar.1 = Q.bar.1,
        k = k,
        w = w,
        p0 = p0,
        p1 = p1,
        m = m))
}

imputeLRT.c <- cmpfun(imputeLRT)

# ========================================================================================
#                           CODE
# ========================================================================================
## Set working directory
# source("I:\\ACERA\\Biofouling\\Biofouling\\papers\\R script New_testLoc.txt", echo=TRUE)
# setwd("C:/data/Publications/Journal articles/Hull fouling risk factors")


# ========================================================================================
#                           Read in data
# ========================================================================================

samplesdata <- read.csv("../data/samples.csv")

propzerodata <- read.csv("../data/aa.csv")

covars <- read.csv("../data/covars.csv")

# samplesdata<-read.table("C:\\Enamul\\ACERA\\Biofouling\\Biofouling\\Analysis_sep_12\\samples.csv", sep=",")

# propzerodata<-read.table("C:\\Enamul\\ACERA\\Biofouling\\Biofouling\\Analysis_sep_12\\aa.csv", sep=",")

# covars<-read.table("C:\\Enamul\\ACERA\\Biofouling\\Biofouling\\Analysis_sep_12\\covars.csv", sep=",")





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

library(Amelia)

samplesdata_testLoc.s$LocID <- factor(samplesdata_testLoc.s$LocID)

set.seed(300)

bf.out <-
    amelia(samplesdata_testLoc.s[,c("wetWeight1","days1","days2","midTrips","paintType",
                                    "LocID","boatID")],
           m = 500,
           cs = "boatID",
           noms = c("paintType","LocID"))

# http://stackoverflow.com/questions/16571580/multi-level-regression-model-on-multiply-imputed-data-set-in-r-amelia-zelig-l

str(bf.out)

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

summary(samples.logit_all2wayglmer1.s)

print(samples.logit_all2wayglmer1.s, corr=FALSE)

### Table of estimates with 95% CI

se1<-sqrt(diag(vcov(samples.logit_all2wayglmer1.s)))

tab1<-cbind(Est=fixef(samples.logit_all2wayglmer1.s), LL=fixef(samples.logit_all2wayglmer1.s)-1.96*se1, UL=fixef(samples.logit_all2wayglmer1.s)+1.96*se1)

exp(tab1)

#####################

tvalue1<-fixef(samples.logit_all2wayglmer1.s)/se1
pvalue1<-2*(1-pnorm(abs(tvalue1)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue1<-sprintf("%.4f", pvalue1)
tab1.1<-cbind(tab1, tvalue1, pvalue1)

tab1.1

## Re-run the analysis but drop days1:paintType
samples.logit_all2wayglmer2.s<-update(samples.logit_all2wayglmer1.s,.~.-days1:paintType)
summary(samples.logit_all2wayglmer2.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer1.s, samples.logit_all2wayglmer2.s, test="F")

## Which suggests to accept the reduced model (p-value 0.0.5111). In the reduced model
## AIC is less than the full model. Thus in terms of
## AIC reduced model (samples.logit_all2wayglmer2.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

print(samples.logit_all2wayglmer2.s, corr=FALSE)

### Table of estimates with 95% CI

se2<-sqrt(diag(vcov(samples.logit_all2wayglmer2.s)))

tab2<-cbind(Est=fixef(samples.logit_all2wayglmer2.s),
            LL=fixef(samples.logit_all2wayglmer2.s)-1.96*se2,
            UL=fixef(samples.logit_all2wayglmer2.s)+1.96*se2)

exp(tab2)

#####################

tvalue2<-fixef(samples.logit_all2wayglmer2.s)/se2
pvalue2<-2*(1-pnorm(abs(tvalue2)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue2<-sprintf("%.4f", pvalue2)
tab2.1<-cbind(tab2, tvalue2, pvalue2)


tab2.1

##########################################

## Now introduce imputation test instead of PB, which chokes because
## of the missingness.

imputeLRT.c(samples.logit_all2wayglmer2.s,  #### SMALLER MODEL FIRST
            samples.logit_all2wayglmer1.s,
            bf.out)

##########################################

## Parametric bootstrap-based test (PB) for comparison --- using reduced data

red.1 <- update(samples.logit_all2wayglmer1.s,
                data = na.omit(samplesdata_testLoc.s[,c("wetWeight1","days1","days2",
                    "midTrips","paintType","LocID","boatID")]))
red.2 <- update(samples.logit_all2wayglmer2.s,
                data = na.omit(samplesdata_testLoc.s[,c("wetWeight1","days1","days2",
                    "midTrips","paintType","LocID","boatID")]))

samples.logit_all2wayglmer.PB1.2.s <- replicate(500, pboot(red.2, red.1))

(samples.logit_all2wayglmer.obsdev1.2.s <- 2*(logLik(red.1)-logLik(red.2)))

mean(samples.logit_all2wayglmer.PB1.2.s > samples.logit_all2wayglmer.obsdev1.2.s)

## This mean is the p-value

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)

#PB shows a p-value 0.838, indicating reduced model (samples.logit_all2wayglmer2.s)
#is better model.

####### Residual Plot ####################
## Full Model (samples.logit_all2wayglmer1.s)
par(mfrow = c(1,2))
scatter.smooth(fitted(samples.logit_all2wayglmer1.s), residuals(samples.logit_all2wayglmer1.s, type="pearson"),
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot GLMER1")
abline(h=0, col="red")

## Reduced Model (samples.logit_all2wayglmer2.s)
scatter.smooth(fitted(samples.logit_all2wayglmer2.s), residuals(samples.logit_all2wayglmer2.s, type="pearson"),
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot GLMER2")
abline(h=0, col="red")

## The residual shows that there are some outliers.
## There is a curved pattern that indicates the residuals are non-linear.

######## Q-Q plot for Residuals #################
## Full Model (samples.logit_all2wayglmer1.s)
par(mfrow = c(1,2))
qqnorm(residuals(samples.logit_all2wayglmer1.s, type="pearson"), main="Q-Q Normal-Residuals GLMER1")
qqline(residuals(samples.logit_all2wayglmer1.s, type="pearson"), col="red")

## Reduced Model (samples.logit_all2wayglmer2.s)
qqnorm(residuals(samples.logit_all2wayglmer2.s, type="pearson"), main="Q-Q Normal-Residuals GLMER2")
qqline(residuals(samples.logit_all2wayglmer2.s, type="pearson"), col="red")

## Q-Q plot for residual shows a problem with the Gaussian assumption for both models.

######## Q-Q Plot for Random effect ######
## Full Model (samples.logit_all2wayglmer1)
## KRH: Why do we need to switch between glmer and lmer??
## APR: accident of Enamul's code --- doesn't make any difference.
## EK: Thanks Andrew.

raf<-ranef(samples.logit_all2wayglmer1.s)
rand.samples.logit_all2wayglmer1.s<-raf[[1]][, "(Intercept)"]
qqnorm(rand.samples.logit_all2wayglmer1.s, main="Q-Q plot, random intercept GLMER1", cex.main = 0.75)
qqline(rand.samples.logit_all2wayglmer1.s)

## Reduced Model (samples.logit_all2wayglmer2.s)
raf<-ranef(samples.logit_all2wayglmer2.s)
rand.samples.logit_all2wayglmer2.s<-raf[[1]][, "(Intercept)"]
qqnorm(rand.samples.logit_all2wayglmer2.s, main="Q-Q plot, random intercept GLMER2", cex.main = 0.75)
qqline(rand.samples.logit_all2wayglmer2.s)

## Q-Q plot for random intercept both for reduced model looks good.
## More comments: From Andrew





######### Further model selection ############

# Dropping midTrips:paintType

samples.logit_all2wayglmer3.s<-update(samples.logit_all2wayglmer1.s,.~.-days1:paintType-midTrips:paintType)
summary(samples.logit_all2wayglmer3.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer2.s,samples.logit_all2wayglmer3.s, test="F")

## Which suggests to accept the reduced model (p-value 0.3027). In the reduced model
## AIC is less than the full model. Thus in terms of
## AIC the model (samples.logit_all2wayglmer3.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

## Parametric bootstrap-based test (PB)

samples.logit_all2wayglmer.obsdev2.3.s <-
  c(2*(logLik(samples.logit_all2wayglmer2.s)-logLik(samples.logit_all2wayglmer3.s)))
set.seed(1001)
samples.logit_all2wayglmer.PB2.3.s<-replicate(500,pboot(samples.logit_all2wayglmer3.s,
                                                      samples.logit_all2wayglmer2.s))
mean(samples.logit_all2wayglmer.PB2.3.s>samples.logit_all2wayglmer.obsdev2.3.s)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)

#PB shows a p-value 0.034, indicating model (samples.logit_all2wayglmer2.s)
#is better model.







 print(samples.logit_all2wayglmer3.s, corr=FALSE)

### Table of estimates with 95% CI


se3<-sqrt(diag(vcov(samples.logit_all2wayglmer3.s)))

tab3<-cbind(Est=fixef(samples.logit_all2wayglmer3.s), LL=fixef(samples.logit_all2wayglmer3.s)-1.96*se3, UL=fixef(samples.logit_all2wayglmer3.s)+1.96*se3)

exp(tab3)

#####################

tvalue3<-fixef(samples.logit_all2wayglmer3.s)/se3
pvalue3<-2*(1-pnorm(abs(tvalue3)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue3<-sprintf("%.4f", pvalue3)
tab3.1<-cbind(tab3, tvalue3, pvalue3)


tab3.1













# Dropping days2:midTrips

samples.logit_all2wayglmer4.s<-update(samples.logit_all2wayglmer1.s,.~.-days1:paintType-midTrips:paintType-days2:midTrips)
summary(samples.logit_all2wayglmer4.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer2.s,samples.logit_all2wayglmer4.s, test="F")

## Which suggests to accept the reduced model (p-value 0.4853). In the reduced model
## AIC is less than the full model. Thus in terms of
## AIC reduced model (samples.logit_all2wayglmer4.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

## Parametric bootstrap-based test (PB)

samples.logit_all2wayglmer.obsdev3.4.s <-
  c(2*(logLik(samples.logit_all2wayglmer3.s)-logLik(samples.logit_all2wayglmer4.s)))
set.seed(1001)
samples.logit_all2wayglmer.PB3.4.s<-replicate(500,pboot(samples.logit_all2wayglmer4.s,
                                                      samples.logit_all2wayglmer3.s))
mean(samples.logit_all2wayglmer.PB3.4.s>samples.logit_all2wayglmer.obsdev3.4.s)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)



#PB shows a p-value 0.384, indicating reduced model (samples.logit_all2wayglmer4.s)
#is better model.





 print(samples.logit_all2wayglmer4.s, corr=FALSE)

### Table of estimates with 95% CI


se4<-sqrt(diag(vcov(samples.logit_all2wayglmer4.s)))

tab4<-cbind(Est=fixef(samples.logit_all2wayglmer4.s), LL=fixef(samples.logit_all2wayglmer4.s)-1.96*se4, UL=fixef(samples.logit_all2wayglmer4.s)+1.96*se4)

exp(tab4)

#####################

tvalue4<-fixef(samples.logit_all2wayglmer4.s)/se4
pvalue4<-2*(1-pnorm(abs(tvalue4)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue4<-sprintf("%.4f", pvalue4)
tab4.1<-cbind(tab4, tvalue4, pvalue4)


tab4.1










# Dropping LocID

samples.logit_all2wayglmer5.s<-update(samples.logit_all2wayglmer1.s,.~.-days1:paintType-midTrips:paintType-days2:midTrips-factor(LocID))


summary(samples.logit_all2wayglmer5.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer4.s,samples.logit_all2wayglmer5.s, test="F")

## Which suggests not to accept the reduced model (p-value  0). In the reduced model
## AIC is more than the full model. Thus in terms of
## AIC the model (samples.logit_all2wayglmer4.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

## Parametric bootstrap-based test (PB)

samples.logit_all2wayglmer.obsdev4.5.s <-
  c(2*(logLik(samples.logit_all2wayglmer4.s)-logLik(samples.logit_all2wayglmer5.s)))
set.seed(1001)
samples.logit_all2wayglmer.PB4.5.s<-replicate(500,pboot(samples.logit_all2wayglmer5.s,
                                                      samples.logit_all2wayglmer4.s))
mean(samples.logit_all2wayglmer.PB4.5.s>samples.logit_all2wayglmer.obsdev4.5.s)



#PB shows a p-value 0.122. As PB is most powerful, we conclude to accept the reduced model (samples.logit_all2wayglmer5.s)









# Dropping days2:paintType

samples.logit_all2wayglmer6.s<-update(samples.logit_all2wayglmer1.s,.~.-days1:paintType-midTrips:paintType-days2:midTrips-days2:paintType)


summary(samples.logit_all2wayglmer6.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer4.s,samples.logit_all2wayglmer6.s, test="F")

## Which suggests merginally not to accept the reduced model (p-value 0.0216). In the reduced model
## AIC is more than the full model. Thus in terms of
## AIC the model (samples.logit_all2wayglmer4.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

## Parametric bootstrap-based test (PB)

samples.logit_all2wayglmer.obsdev4.6.s <-
  c(2*(logLik(samples.logit_all2wayglmer4.s)-logLik(samples.logit_all2wayglmer6.s)))
set.seed(1001)
samples.logit_all2wayglmer.PB4.6.s<-replicate(500,pboot(samples.logit_all2wayglmer6.s,
                                                      samples.logit_all2wayglmer4.s))
mean(samples.logit_all2wayglmer.PB4.6.s>samples.logit_all2wayglmer.obsdev4.6.s)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)












#PB shows a p-value 0.118, indicating reduced model (samples.logit_all2wayglmer6.s)
#is better model


# Dropping days1:days2

samples.logit_all2wayglmer7.s<-update(samples.logit_all2wayglmer1.s,.~.-days2:midTrips-days1:midTrips-midTrips:paintType-days2:paintType-days1:days2)


summary(samples.logit_all2wayglmer7.s)

## The reduced model can be tested against the full model by using following anova
anova(samples.logit_all2wayglmer6.s,samples.logit_all2wayglmer7.s, test="F")

## Which suggests to accept the reduced model (p-value 0.2976). In the reduced model
## AIC is less than the full model. Thus in terms of
## AIC the model (samples.logit_all2wayglmer7.s) is the better model.
## However, we perform a parametric bootstrap test to take a final decision.

## Parametric bootstrap-based test (PB)

samples.logit_all2wayglmer.obsdev6.7.s <-
  c(2*(logLik(samples.logit_all2wayglmer6.s)-logLik(samples.logit_all2wayglmer7.s)))
set.seed(1001)
samples.logit_all2wayglmer.PB6.7.s<-replicate(500,pboot(samples.logit_all2wayglmer7.s,
                                                      samples.logit_all2wayglmer6.s))
mean(samples.logit_all2wayglmer.PB6.7.s>samples.logit_all2wayglmer.obsdev6.7.s)

## Classical test: chisq(samples.logit_all2wayglmer.obsdev1.2, df=63, lower.tail=FALSE)












#PB shows a p-value 0.35, indicating reduced model (samples.logit_all2wayglmer7.s)
#is better model




## No more dropping. The model samples.logit_all2wayglmer7.s is the final model.














## Thus for GLMM, samples.logit_all2wayglmer4.s is the best model which consists of the
## variabls: days1, days2, midTrips, paintType, LocID, days1:days2, days1:midTrips, days2:paintType.



# So the model samples.logit_all2wayglmer4.s
# is probably the best model to use, and this model suggest that days2 is statistically significant
# for fouling biomass. midTrips is significant at 5% level. LocID is also significent.

## The summary of full model suggests dropping days2:midTrips. Refitting without
## this term, and then examining the significance of terms in the resulting model,
## suggests dropping days1:paintType. Continuing in the same way, we conclude that
## samples.logit_all2wayglmer4.s is probably the best model to use.





########## Residual Plot  ####################
## Reduced Model (samples.logit_all2wayglmer4.s)
par(mfrow = c(1,3))
scatter.smooth(fitted(samples.logit_all2wayglmer4.s), residuals(samples.logit_all2wayglmer4.s, type="pearson"),
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot for the final model", cex.main = 0.8)
abline(h=0, col="red")

## The residual shows that there are some outliers. There is a curved pattern that indicates
## the residuals are non-linear. It is likely that plot of of this model seems to be better
## than the full model.

######## Q-Q plot for Residuals ###########
## Reduced Model (samples.logit_all2wayglmer4.s)
qqnorm(residuals(samples.logit_all2wayglmer4.s, type="pearson"), main="Q-Q Normal-Residuals for the final model", cex.main = 0.8)
qqline(residuals(samples.logit_all2wayglmer4.s, type="pearson"), col="red")

## Q-Q plot for residual shows with the Gaussian assumption.
## It seems that the samples.logit_all2wayglmer4.s is good.
## The Gaussian assumption is more appropriate for samples.logit_all2wayglmer4.s than the full model.

#######  Q-Q Plot for Random effect #######
## Reduced Model (samples.logit_all2wayglmer4.s)
raf<-ranef(samples.logit_all2wayglmer4.s)
rand.samples.logit_all2wayglmer4.s<-raf[[1]][, "(Intercept)"]
qqnorm(rand.samples.logit_all2wayglmer4.s, main="Q-Q plot for random intercept (boat) for the final model", cex.main = 0.8)
qqline(rand.samples.logit_all2wayglmer4.s)

## Q-Q plot for random intercept both for the reduced model looks good
## (better than the full model).



# =======================================================================================================
#                         Non-linear Relationships
# =======================================================================================================



# There is no guarantee that the relationship between the response variable and the continuous predictors (days1 and days2)
# should be linear, in fact, such a claim would tend to lack credibility. We need a model that will permit the
# detection of curvature, but not enforce it. The best approach is to use splines. Here, we use a cubic B-spline.


library(splines)


######### for Full Model ################

bs_samples.logit_all2wayglmer1.s<-glmer(wetWeight1~bs(days1)+bs(days2)+bs(midTrips)+paintType+factor(LocID)+days1:days2+
	days1:midTrips+ days1:paintType+days2:midTrips+days2:paintType+
		midTrips:paintType+(1|boatID),family=binomial,
			data=samplesdata_testLoc.s, na.action=na.exclude)



summary(bs_samples.logit_all2wayglmer1.s)







 print(bs_samples.logit_all2wayglmer1.s, corr=FALSE)

### Table of estimates with 95% CI


bs_se1<-sqrt(diag(vcov(bs_samples.logit_all2wayglmer1.s)))

bs_tab1<-cbind(Est=fixef(bs_samples.logit_all2wayglmer1.s), LL=fixef(bs_samples.logit_all2wayglmer1.s)-1.96*bs_se1, UL=fixef(bs_samples.logit_all2wayglmer1.s)+1.96*bs_se1)

exp(bs_tab1)

#####################

bs_tvalue1<-fixef(bs_samples.logit_all2wayglmer1.s)/bs_se1
bs_pvalue1<-2*(1-pnorm(abs(bs_tvalue1)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

bs_pvalue1<-sprintf("%.4f", bs_pvalue1)
bs_tab1.1<-cbind(bs_tab1, bs_tvalue1, bs_pvalue1)


bs_tab1.1




# Tried to fit the model with bs(days1) and bs(days2) with the interaction terms but got an error massage in summary
# Error in asMethod(object) : matrix is not symmetric [1,2]



########### for best model ###############


bs_samples.logit_all2wayglmer4.s<-glmer(wetWeight1~bs(days1)+bs(days2)+midTrips+factor(paintType)+factor(LocID)+bs(days1):bs(days2)+bs(days1):midTrips+bs(days2):paintType
				+(1|boatID),family=binomial,
					data=samplesdata_testLoc.s, na.action=na.exclude)

summary(bs_samples.logit_all2wayglmer4.s)






 print(bs_samples.logit_all2wayglmer4.s, corr=FALSE)

### Table of estimates with 95% CI


bs_se4<-sqrt(diag(vcov(bs_samples.logit_all2wayglmer4.s)))

bs_tab4<-cbind(Est=fixef(bs_samples.logit_all2wayglmer4.s), LL=fixef(bs_samples.logit_all2wayglmer4.s)-1.96*bs_se4, UL=fixef(bs_samples.logit_all2wayglmer4.s)+1.96*bs_se4)

exp(bs_tab4)

#####################

bs_tvalue4<-fixef(bs_samples.logit_all2wayglmer4.s)/bs_se4
bs_pvalue4<-2*(1-pnorm(abs(bs_tvalue4)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

bs_pvalue4<-sprintf("%.4f", bs_pvalue4)
bs_tab4.1<-cbind(bs_tab4, bs_tvalue4, bs_pvalue4)


bs_tab4.1







########## Residual Plot  ####################

## best Model (bs_samples.logit_all2wayglmer4.s)
par(mfrow = c(1,2))
scatter.smooth(fitted(bs_samples.logit_all2wayglmer4.s), residuals(bs_samples.logit_all2wayglmer4.s, type="pearson"),
	xlab="Fitted Values", ylab="Residuals", main="Residuals Plot GLMER4", cex.main = 0.8)
abline(h=0, col="red")


######## Q-Q plot for Residuals ###########
## best Model (bs_samples.logit_all2wayglmer4.s)

qqnorm(residuals(bs_samples.logit_all2wayglmer4.s, type="pearson"), main="Q-Q Normal-Residuals GLMER4", cex.main = 0.8)
qqline(residuals(bs_samples.logit_all2wayglmer4.s, type="pearson"), col="red")


#######  Q-Q Plot for Random effect #######
## best Model (bs_samples.logit_all2wayglmer4.s)
bs_raf<-ranef(bs_samples.logit_all2wayglmer4.s)
bs_rand.samples.logit_all2wayglmer4.s<-bs_raf[[1]][, "(Intercept)"]
qqnorm(bs_rand.samples.logit_all2wayglmer4.s, main="Q-Q plot for random intercept GLMER4", cex.main = 0.8)
qqline(bs_rand.samples.logit_all2wayglmer4.s)






################# Figure 9-11 #############

## Relationship between the proportions of zero samples in each location and the three potential explanatory variables,
## days1, days2 and midTrips. Data set  propzerodata_testLoc

## These figures suggest that days1 and midTrips have very little predictive power in this context. In this
## initial inspection the predictive power of days2 appeared to be more promising.
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

## KRH: Repeating this analysis with standard sampling units doesn't change the message

## EK: yes, they are same. I think we can do for testLoc and make a conclusion that they are same with standard sampling units.


newy1a<-plot(propzerodata_testLoc1$days1, propzerodata_testLoc1$n0/propzerodata_testLoc1$n, xlab="Number of days since the vessel was last used",
  ylab="Proportion of zero samples", main="")
fit1a<-lsfit(propzerodata_testLoc1$days1, propzerodata_testLoc1$n0/propzerodata_testLoc1$n)
abline(fit1a)

newy2a<-plot(propzerodata_testLoc1$days2, propzerodata_testLoc1$n0/propzerodata_testLoc1$n, xlab="Number of days since the vessel was last cleaned ot antifouled",
	ylab="Proportion of zero samples", main="")
fit2a<-lsfit(propzerodata_testLoc1$days2, propzerodata_testLoc1$n0/propzerodata_testLoc1$n)
abline(fit2a)

newy3a<-plot(propzerodata_testLoc1$midTrips, propzerodata_testLoc1$n0/propzerodata_testLoc1$n, xlab="Median number of trips per annum",
	ylab="Proportion of zero samples", main="")
fit3a<-lsfit(propzerodata_testLoc1$midTrips, propzerodata_testLoc1$n0/propzerodata_testLoc1$n)
abline(fit3a)



# =======================================================================================================
#                        ####### Quasi-likelihood Analysis ################
# =======================================================================================================


## Quasi likelihood analysis to test the predictive power of days1, days2 and midTrips.
glm_quasi<-glm(n0/n~days1+days2+midTrips, data=propzerodata_testLoc, family=quasibinomial, method="glm.fit", weights = n, trace=F)
summary(glm_quasi)

## Our initial deduction is also supported by the quasi likelihood analysis that the predictive power of days2 is
## more convincing as days2 is statistically significant for the proportion of non-zero samples.


# EK:Repeating this analysis with standard sampling units doesn't change the message


glm_quasia<-glm(n0/n~days1+days2+midTrips, data=propzerodata_testLoc1, family=quasibinomial, method="glm.fit", weights = n, trace=F)
summary(glm_quasia)





############# Figure 12 #################
## Actual and predicted proportion of zero samples based on quasi-likelihood analysis.
par(mfrow=c(2,1))

glm_quasi_days2 <-
  glm(n0/n~days2, data=propzerodata_testLoc, family=quasibinomial, weights=n, method="glm.fit", trace=F)
newx<-seq(0, 1000, len=500)
newy<-predict(glm_quasi_days2, newdata=data.frame(days2=newx), type="link", se.fit=T, conf.level=0.95)
plot(propzerodata_testLoc$days2, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Days since vessel last cleaned", ylab="Proportion of clean samples")
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]]+2*newy$se.fit), lty=2)
lines(newx, quasibinomial()$linkinv(newy[[1]]-2*newy$se.fit), lty=2)

## The graph shows that the proportion of clean samples diminishes very quickly after vessel has been clean
## or antifouled.
## KRH: But is this result driven largely by just three vessels ??
## APR: If we remove those three vessels we get very similar results.
## EK: yes, so we can use glm_quasi_days2 and make a conclusion that we get similar results if we remove those three vassels.

glm_quasi_days3 <-
  glm(n0/n~days2, data=propzerodata_testLoc, subset = days2 < 600,
      family=quasibinomial, weights=n, method="glm.fit", trace=F)
newx<-seq(0, 1000, len=500)
newy<-predict(glm_quasi_days3, newdata=data.frame(days2=newx), type="link", se.fit=T, conf.level=0.95)
plot(propzerodata_testLoc$days2, propzerodata_testLoc$n0/propzerodata_testLoc$n, xlab="Days since vessel last cleaned", ylab="Proportion of clean samples")
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]]+2*newy$se.fit), lty=2)
lines(newx, quasibinomial()$linkinv(newy[[1]]-2*newy$se.fit), lty=2)


# =========================================================================================================
#           Linear Mixed effect Model (LMM) for non-zero wet weight taken only from standard sampling units
#                           Dataset is nonzerosamples1_testLoc1, outliers removed
# =========================================================================================================

## KRH: Note also this analysis could be re-done with all the standard sampling units
## EK: Done

#####   LMM ######
## LMM is used as considering boatID as a random effect. The analysis is done by removing the outliers.






## Histograms of sample wet weight (grams) for the non-zero (biomass) sample from all the standrad sampling units (HA, PJ, HP), and the log of
## sample wet weight (grams) for the non-zero (biomass) samples from those locations. The wet weight fouling biomass of
## non-zero samples from standard sampling units is approximately log normally distributed but the log of these non-zero biomasses
## is normally distributed.


par(mfrow=c(1,2))
hist(nonzerosamples1_testLoc1$wetWeight, prob=F, xlab="Wet-weight of non-zero samples (grams)",
     ylab="Frequency", main="" )
##text(7000,100, "a")
hist(log(nonzerosamples_testLoc1$wetWeight), prob=F, xlab="Log wet-weight of non-zero samples",
     ylab="Frequency", main="" )
##text(9,80, "b")








## We first fit a LMM considering all fixed effects and their interaction.

#nonzerosamples1_testLoc1$boatID_id <- factor(nonzerosamples1_testLoc1$boatID)
#nonzerosamples1_testLoc1_new <- groupedData(wetWeight~days2|boatID_id, data=nonzerosamples1_testLoc1)


#lww_testLoc1.lmm1 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+days1:paintType+
	#days2:midTrips+days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)






lww_testLoc1.lmm1 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+days1:paintType+
	 days2:midTrips+days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)




summary(lww_testLoc1.lmm1)
anova(lww_testLoc1.lmm1)







 print(lww_testLoc1.lmm1, corr=FALSE)

### Table of estimates with 95% CI


se1<-sqrt(diag(vcov(lww_testLoc1.lmm1)))

tab1<-cbind(Est=fixef(lww_testLoc1.lmm1), LL=fixef(lww_testLoc1.lmm1)-1.96*se1, UL=fixef(lww_testLoc1.lmm1)+1.96*se1)

exp(tab1)

#####################

tvalue1<-fixef(lww_testLoc1.lmm1)/se1
pvalue1<-2*(1-pnorm(abs(tvalue1)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue1<-sprintf("%.4f", pvalue1)
tab1.1<-cbind(tab1, tvalue1, pvalue1)


tab1.1








## ANOVA suggests that LocID is statistically significant. However days2 is
## significant at 5% level of significance.









############### Residual plot for the model (lww_testLoc1.lmm1) ##########

## plot(lww_testLoc1.lmm1)	# More convenient.

scatter.smooth(fitted(lww_testLoc1.lmm1), residuals(lww_testLoc1.lmm1, type="pearson"), xlab="Fitted Values",
	ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")


## Variance-adjusted, outermost residuals have approximately constant variance as the residuals
## do not have any relationship. Theare is no curved patern that indicates that the residual are
## normlly distributed. There are no substantial outliers. Thus it can be said that LMM model fits the
## data well.

############ Q-Q plot Figure 14 ##############
## Q-Q plot for Random effect


#ref.boat<-ranef(lww_testLoc1.lmm1)[[1]]
#ref.var.boat<-tapply(residuals(lww_testLoc1.lmm1,type="pearson", level=1), nonzerosamples1_testLoc1_new$boatID, var)
#qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
#qqline(ref.boat, col="red")




ref.boat<-ranef(lww_testLoc1.lmm1)[[1]]
qqnorm(ranef(lww_testLoc1.lmm1)$boatID$"(Intercept)", main="Q-Q plot for the random effect" )
qqline(ref.boat, col="red")



## Q-Q plot for Residuals ################

qqnorm(residuals(lww_testLoc1.lmm1, type="pearson"), main="Q-Q Normal-Residuals")
qqline(residuals(lww_testLoc1.lmm1, type="pearson"), col="red")

## The Q-Q plot for residuals follows approximately in the line, that means two distributions agree
## after linearly transforming the values. However, Q-Q plot for for random effect slightly differ
## from straight line.


######## Selecting the best model by using LMM ###########



## lww_testLoc1.lmm1 suggests dropping days2:midTrips. We refitt the model without this term and compare with
## the full model by using AIC.

#lww_testLoc1.lmm2 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	#days2:midTrips+days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)




lww_testLoc1.lmm2 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+days1:paintType+
	 +days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)



summary(lww_testLoc1.lmm2)
anova(lww_testLoc1.lmm2)










 print(lww_testLoc1.lmm2, corr=FALSE)

### Table of estimates with 95% CI


se2<-sqrt(diag(vcov(lww_testLoc1.lmm2)))

tab2<-cbind(Est=fixef(lww_testLoc1.lmm2), LL=fixef(lww_testLoc1.lmm2)-1.96*se2, UL=fixef(lww_testLoc1.lmm2)+1.96*se2)

exp(tab2)

#####################

tvalue2<-fixef(lww_testLoc1.lmm2)/se2
pvalue1<-2*(1-pnorm(abs(tvalue2)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue2<-sprintf("%.4f", pvalue2)
tab2.1<-cbind(tab2, tvalue2, pvalue2)


tab2.1



## The anova(lww_testLoc1.lmm1, lww_testLoc1.lmm2) shows p-value 0.9085, suggesting to accept reduced
## model). However, AIC(lww_testLoc1.lmm1, lww_testLoc1.lmm2) is greater for full model (679 compared to 677).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm1r <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+days1:days2+days1:midTrips+days1:paintType+
	#days2:midTrips+days2:paintType+midTrips:paintType+ (1|boatID),
		#data=nonzerosamples1_HA, na.action=na.exclude)

#lww_testLoc1.lmm2r <- update(lww_testLoc1.lmm1r,.~.-days1:paintType)





lww_testLoc1.lmm1.2.obsdev <- c(2*(logLik(lww_testLoc1.lmm1)-logLik(lww_testLoc1.lmm2)))
set.seed(1001)
lww_testLoc1.lmm1.2.PB1 <- replicate(500,pboot(lww_testLoc1.lmm2, lww_testLoc1.lmm1))
mean(lww_testLoc1.lmm1.2.PB1> lww_testLoc1.lmm1.2.obsdev)





## PB shows a p-value 0.748, indicating the model (lww_testLoc1.lmm2) is better model.






## Thus we keep the model lww_testLoc1.lmm2 and try testing dropping day2:LocID





lww_testLoc1.lmm3 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+days1:paintType+
	 +days2:paintType+midTrips:paintType+LocID:days1+LocID:midTrips+LocID:paintType  +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)



# lww_testLoc1.lmm3 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	# days2:midTrips+days2:paintType+midTrips:paintType+LocID:days1+LocID:midTrips+LocID:paintType,
		# random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)



summary(lww_testLoc1.lmm3)
anova(lww_testLoc1.lmm3)




## The anova(lww_testLoc1.lmm2, lww_testLoc1.lmm3) suggests p-value  0.003495, suggesting to reject reduced
## model). AIC(lww_testLoc1.lmm2, lww_testLoc1.lmm3) is greater for reduced model (684 compared to 677).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

# lww_testLoc1.lmm3r <- update(lww_testLoc1.lmm1r,.~.-days1:paintType-LocID:days2)


lww_testLoc1.lmm2.3.obsdev <- c(2*(logLik(lww_testLoc1.lmm2)-logLik(lww_testLoc1.lmm3)))
set.seed(1001)
lww_testLoc1.lmm2.3.PB1 <- replicate(500,pboot(lww_testLoc1.lmm3, lww_testLoc1.lmm2))
mean(lww_testLoc1.lmm2.3.PB1> lww_testLoc1.lmm2.3.obsdev)


## PB shows a p-value 0.004, indicating the model (lww_testLoc1.lmm2) is better model.


## Thus we keep the model lww_testLoc1.lmm2 and try testing dropping days2:paintType






#lww_testLoc1.lmm4 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	#days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)




lww_testLoc1.lmm4 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+days1:paintType+
	 midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)







summary(lww_testLoc1.lmm4)
anova(lww_testLoc1.lmm4)






 print(lww_testLoc1.lmm4, corr=FALSE)

### Table of estimates with 95% CI


se4<-sqrt(diag(vcov(lww_testLoc1.lmm4)))

tab4<-cbind(Est=fixef(lww_testLoc1.lmm4), LL=fixef(lww_testLoc1.lmm4)-1.96*se4, UL=fixef(lww_testLoc1.lmm4)+1.96*se4)

exp(tab4)

#####################

tvalue4<-fixef(lww_testLoc1.lmm4)/se4
pvalue4<-2*(1-pnorm(abs(tvalue4)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue4<-sprintf("%.4f", pvalue4)
tab4.1<-cbind(tab4, tvalue4, pvalue4)


tab4.1






## The anova(lww_testLoc1.lmm2, lww_testLoc1.lmm4) suggests p-value 0.1087, suggesting to accept reduced
## model). AIC(lww_testLoc1.lmm2, lww_testLoc1.lmm3) is smilar.
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm4r <- update(lww_testLoc1.lmm2r,.~.-days2:midTrips)


lww_testLoc1.lmm2.4.obsdev <- c(2*(logLik(lww_testLoc1.lmm2)-logLik(lww_testLoc1.lmm4)))
set.seed(1001)
lww_testLoc1.lmm2.4.PB1 <- replicate(500,pboot(lww_testLoc1.lmm4, lww_testLoc1.lmm2))
mean(lww_testLoc1.lmm2.4.PB1> lww_testLoc1.lmm2.4.obsdev)


## PB shows a p-value 0.002, indicating the model (lww_testLoc1.lmm2) is better model.


## Thus we keep the model lww_testLoc1.lmm4 and try testing dropping days1:paintType






#lww_testLoc1.lmm5 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	#days2:paintType+midTrips:paintType+LocID:days1+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)




lww_testLoc1.lmm5 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	 midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)



summary(lww_testLoc1.lmm5)
anova(lww_testLoc1.lmm5)





 print(lww_testLoc1.lmm5, corr=FALSE)

### Table of estimates with 95% CI


se5<-sqrt(diag(vcov(lww_testLoc1.lmm5)))

tab5<-cbind(Est=fixef(lww_testLoc1.lmm5), LL=fixef(lww_testLoc1.lmm5)-1.96*se5, UL=fixef(lww_testLoc1.lmm5)+1.96*se5)

exp(tab5)

#####################

tvalue5<-fixef(lww_testLoc1.lmm5)/se5
pvalue5<-2*(1-pnorm(abs(tvalue5)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue5<-sprintf("%.4f", pvalue5)
tab5.1<-cbind(tab5, tvalue5, pvalue5)


tab5.1





## The anova(lww_testLoc1.lmm4, lww_testLoc1.lmm5) suggests p-value .6221, indicating to accept reduced
## model). AIC is also less for the reduced model (674 compared to 677).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

# lww_testLoc1.lmm5r <- update(lww_testLoc1.lmm4r,.~.-LocID:days2)

lww_testLoc1.lmm4.5.obsdev <- c(2*(logLik(lww_testLoc1.lmm4)-logLik(lww_testLoc1.lmm5)))
set.seed(1001)
lww_testLoc1.lmm4.5.PB1 <- replicate(500,pboot(lww_testLoc1.lmm5, lww_testLoc1.lmm4))
mean(lww_testLoc1.lmm4.5.PB1> lww_testLoc1.lmm4.5.obsdev)




## PB shows a p-value .744, indicating the model (lww_testLoc1.lmm5) is better model.


## Thus we keep the model lww_testLoc1.lmm5 and try testing dropping days1:LocID






#lww_testLoc1.lmm6 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	#days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)




lww_testLoc1.lmm6 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:days2+days1:midTrips+
	 midTrips:paintType+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)




summary(lww_testLoc1.lmm6)
anova(lww_testLoc1.lmm6)




## The anova(lww_testLoc1.lmm5, lww_testLoc1.lmm6) suggests p-value 0, indicating to reject reduced
## model). AIC for reduced model is also more than the full model (674 compared to 690).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm6r <- update(lww_testLoc1.lmm4r,.~.-days1:days2)

lww_testLoc1.lmm5.6.obsdev <- c(2*(logLik(lww_testLoc1.lmm5)-logLik(lww_testLoc1.lmm6)))
set.seed(1001)
lww_testLoc1.lmm5.6.PB1 <- replicate(500,pboot(lww_testLoc1.lmm6, lww_testLoc1.lmm5))
mean(lww_testLoc1.lmm5.6.PB1> lww_testLoc1.lmm5.6.obsdev)





## PB shows a p-value 0.93, indicating the model (lww_testLoc1.lmm6) is better model.


## Thus we keep the model lww_testLoc1.lmm5 and try testing dropping days1:days2






lww_testLoc1.lmm7 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	 midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)






#lww_testLoc1.lmm7 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	#days2:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)




summary(lww_testLoc1.lmm7)
anova(lww_testLoc1.lmm7)





 print(lww_testLoc1.lmm7, corr=FALSE)

### Table of estimates with 95% CI


se7<-sqrt(diag(vcov(lww_testLoc1.lmm7)))

tab7<-cbind(Est=fixef(lww_testLoc1.lmm7), LL=fixef(lww_testLoc1.lmm7)-1.96*se7, UL=fixef(lww_testLoc1.lmm7)+1.96*se7)

exp(tab7)

#####################

tvalue7<-fixef(lww_testLoc1.lmm7)/se7
pvalue7<-2*(1-pnorm(abs(tvalue7)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue7<-sprintf("%.4f", pvalue7)
tab7.1<-cbind(tab7, tvalue7, pvalue7)


tab7.1






## The anova(lww_testLoc1.lmm5, lww_testLoc1.lmm7) suggests p-value .6912, indicating marginally  accept reduced
## model). However, AIC(lww_testLoc1.lmm6, lww_testLoc1.lmm7) is less for reduced model (672 compared to 674).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm7r <- update(lww_testLoc1.lmm6r,.~.-midTrips:paintType)

lww_testLoc1.lmm5.7.obsdev <- c(2*(logLik(lww_testLoc1.lmm5)-logLik(lww_testLoc1.lmm7)))
set.seed(1001)
lww_testLoc1.lmm5.7.PB1 <- replicate(2000,pboot(lww_testLoc1.lmm5, lww_testLoc1.lmm7))
mean(lww_testLoc1.lmm5.7.PB1> lww_testLoc1.lmm5.7.obsdev)



## For 500 time replication PB shows a p-value 0.112, indicating to accept reduced model (lww_testLoc1.lmm7). As AIC is slightly more, we replicate the modeling
## 2000 times and shows a p- value 0.109. Thus we conclude that the reduced model (lww_testLoc1.lmm7) is better model.



## Thus we keep the model lww_testLoc1.lmm7 and try testing dropping paintType:LocID



#lww_testLoc1.lmm8 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	#days2:paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)


lww_testLoc1.lmm8 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	 midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)




summary(lww_testLoc1.lmm8)
anova(lww_testLoc1.lmm8)




## The anova(lww_testLoc1.lmm7, lww_testLoc1.lmm8) suggests to reject the reduced model p-value 0.
## However, AIC(lww_testLoc1.lmm7, lww_testLoc1.lmm8) is greater for reduced model (672 compared to 684).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm8r <- update(lww_testLoc1.lmm6r,.~.-LocID:midTrips)

lww_testLoc1.lmm7.8.obsdev <- c(2*(logLik(lww_testLoc1.lmm7)-logLik(lww_testLoc1.lmm8)))
set.seed(1001)
lww_testLoc1.lmm7.8.PB1 <- replicate(2000,pboot(lww_testLoc1.lmm8, lww_testLoc1.lmm7))
mean(lww_testLoc1.lmm7.8.PB1> lww_testLoc1.lmm7.8.obsdev)




## PB shows a p-value 0.182, indicating the model (lww_testLoc1.lmm8) is better model.




## Thus we keep the model lww_testLoc1.lmm7 and try testing dropping days1:midTrips


#lww_testLoc1.lmm9 <- lme(log(wetWeight)~days1+days2+midTrips+paintType+LocID+days1:midTrips+
	#midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)





lww_testLoc1.lmm9 <- lmer(log(wetWeight)~days1+days2+midTrips+paintType+LocID+
	 midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType+
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)





summary(lww_testLoc1.lmm9)
anova(lww_testLoc1.lmm9)







 print(lww_testLoc1.lmm9, corr=FALSE)

### Table of estimates with 95% CI


se9<-sqrt(diag(vcov(lww_testLoc1.lmm9)))

tab9<-cbind(Est=fixef(lww_testLoc1.lmm9), LL=fixef(lww_testLoc1.lmm9)-1.96*se9, UL=fixef(lww_testLoc1.lmm9)+1.96*se9)

exp(tab9)

#####################

tvalue9<-fixef(lww_testLoc1.lmm9)/se9
pvalue9<-2*(1-pnorm(abs(tvalue9)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

pvalue9<-sprintf("%.4f", pvalue9)
tab9.1<-cbind(tab9, tvalue9, pvalue9)


tab9.1

















## The anova(lww_testLoc1.lmm7, lww_testLoc1.lmm9) suggests to accept the reduced model (p-value 0.2865),
## although AIC is same for reduced model (672 compared to 672).
## We perform a bootstrap test for final decision.

## Parametric Bootstrap test (PB) ####

#lww_testLoc1.lmm9r <- update(lww_testLoc1.lmm6r,.~.-days2:paintType)

lww_testLoc1.lmm7.9.obsdev <- c(2*(logLik(lww_testLoc1.lmm7)-logLik(lww_testLoc1.lmm9)))
set.seed(1001)
lww_testLoc1.lmm7.9.PB1 <- replicate(5000,pboot(lww_testLoc1.lmm9, lww_testLoc1.lmm7))
mean(lww_testLoc1.lmm7.9.PB1> lww_testLoc1.lmm7.9.obsdev)









## PB shows a p-value 0.076, indicating the model (lww_testLoc1.lmm9) is better model.



## Thus we keep the model lww_testLoc1.lmm9. Some of the interactions are insignificant but PB test
## suggests to include them in the model. paintType is not significant but we can not drop it as its
## interaction with other variables are included in the model.




#####################################################  lww_testLoc1.lmm9 would be the best model to be used.
## Thus lww_testLoc1.lmm9 would be the best model to be used.The variables included in the model are
## days1, days2, midTrips, paintType, LocID, LocID:days1, LocID:days2,
## LocID:midTrips, LocID:paintType,





## Days2 and LocID are strongly significant.
## PaintType itself is not a significant but its interaction with LocID is statistically significant.
## If the vessel is not cleaned one day, then the log of the wet weight is increased by .003 grams in testLoc1 given
## an ablative antifouling paint is used and other variables are constant.


## On the other hand, if ablative
## antifouling paint is used the log of the wet weight is decreased by .0356 grams in testLoc1, if median number of
## trips is increased by one day when other variables are constant. paintType influence the effect of
## vessel activity on the log of the wet weight of biofouling. On an average midTrips, hard
## paint is significant on fouling biomass campared to ablative paint.



############### Residual plot, for the model (lww_testLoc1.lmm9) ##########
scatter.smooth(fitted(lww_testLoc1.lmm9), residuals(lww_testLoc1.lmm9, type="pearson"), xlab="Fitted Values",
	ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")

## Residual plot looks good. Thus the model is reasonable.
## KRH: same plot for lw_HA.lmm9 is virtually identical

## scatter.smooth(fitted(lww_HA.lmm9), residuals(lww_HA.lmm9, type="pearson"), xlab="Fitted Values",
## ylab="Residuals", main="Residuals Plot")
## abline(h=0, col="red")


############ Q-Q plot Figure 14 ##############

## Q-Q plot for Random effect

# ref.boat<-ranef(lww_testLoc1.lmm10)[[1]]
# ref.var.boat<-tapply(residuals(lww_testLoc1.lmm10,type="pearson", level=1), nonzerosamples1_testLoc1_new$boatID_id, var)
# qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
# qqline(ref.boat, col="red")




ref.boat<-ranef(lww_testLoc1.lmm9)[[1]]
qqnorm(ranef(lww_testLoc1.lmm9)$boatID$"(Intercept)", main="Q-Q plot for the random effect" )
qqline(ref.boat, col="red")







## KRH And for model 9 - slightly better fit?
## ref.boat<-ranef(lww_HA.lmm9)[[1]]
## ref.var.boat<-tapply(residuals(lww_HA.lmm9,type="pearson", level=1), nonzerosamples1_HA_new$boatID_id, var)
## qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
## qqline(ref.boat, col="red")


##### Q-Q plot for Residuals ################

qqnorm(residuals(lww_testLoc1.lmm9, type="pearson"), main="Q-Q Normal-Residuals")
qqline(residuals(lww_testLoc1.lmm9, type="pearson"), col="red")

## KRH And for model 9, again a slightly better fit?
## qqnorm(residuals(lww_HA.lmm9, type="pearson"), main="Q-Q Normal-Residuals")
## qqline(residuals(lww_HA.lmm9, type="pearson"), col="red")

## The Q-Q plot for residuals follows approximately in the line, that means two distributions agree
## after linearly transforming the values. However, Q-Q plot for for random effect slightly differ
## from straight line but more reasonable than before and thus the model is appropriate.








# =======================================================================================================
#                         Non-linear Relationships
# =======================================================================================================


# For full model







#bs_lww_testLoc1.lmm1 <- lme(log(wetWeight)~bs(days1)+bs(days2)+midTrips+paintType+LocID+bs(days1):bs(days2)+bs(days1):midTrips+bs(days1):paintType+
	#bs(days2):midTrips+bs(days2):paintType+midTrips:paintType+LocID:bs(days1)+LocID:bs(days2)+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)


bs_lww_testLoc1.lmm1 <- lmer(log(wetWeight)~bs(days1)+bs(days2)+midTrips+paintType+LocID+bs(days1):bs(days2)+bs(days1):midTrips+bs(days1):paintType+
	 bs(days2):midTrips+bs(days2):paintType+midTrips:paintType+LocID:days1+LocID:days2+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)





####### Error Massages: Error in MEEM(object, conLin, control$niterEM) :
####### Singularity in backsolve at level 0, block 1



summary(bs_lww_testLoc1.lmm1)
anova(bs_lww_testLoc1.lmm1)






# For reduced model



#nonzerosamples1_testLoc1$boatID_id <- factor(nonzerosamples1_testLoc1$boatID)
#nonzerosamples1_testLoc1_new <- groupedData(wetWeight~days2|boatID_id, data=nonzerosamples1_testLoc1)


#bs_lww_testLoc1.lmm10 <- lme(log(wetWeight)~bs(days1)+bs(days2)+midTrips+paintType+LocID+
	#LocID:bs(days1)+LocID:bs(days2)+LocID:midTrips+LocID:paintType,
		#random=~1|boatID_id, data=nonzerosamples1_testLoc1_new, na.action=na.exclude)







bs_lww_testLoc1.lmm9 <- lmer(log(wetWeight)~bs(days1)+bs(days2)+midTrips+paintType+LocID+
	 	LocID:bs(days1)+LocID:bs(days2)+LocID:midTrips+LocID:paintType +
		 (1|boatID), REML=F, data=nonzerosamples1_testLoc1, na.action=na.exclude)






summary(bs_lww_testLoc1.lmm9)
anova(bs_lww_testLoc1.lmm9)


## Results are same as before







 print(bs_lww_testLoc1.lmm9, corr=FALSE)

### Table of estimates with 95% CI


bs_se9<-sqrt(diag(vcov(bs_lww_testLoc1.lmm9)))

bs_tab9<-cbind(Est=fixef(bs_lww_testLoc1.lmm9), LL=fixef(bs_lww_testLoc1.lmm9)-1.96*bs_se9, UL=fixef(bs_lww_testLoc1.lmm9)+1.96*bs_se9)

exp(bs_tab9)

#####################

bs_tvalue9<-fixef(bs_lww_testLoc1.lmm9)/bs_se9
bs_pvalue9<-2*(1-pnorm(abs(bs_tvalue9)))

#pvalue1<-sprintf("%.4f", 2*(1-pnorm(abs(tvalue1))))

bs_pvalue9<-sprintf("%.4f", bs_pvalue9)
bs_tab9.1<-cbind(bs_tab9, bs_tvalue9, bs_pvalue9)


bs_tab9.1




############### Residual plot, non-linear relationship ##########

scatter.smooth(fitted(bs_lww_testLoc1.lmm9), residuals(bs_lww_testLoc1.lmm9, type="pearson"), xlab="Fitted Values",
	ylab="Residuals", main="Residuals Plot")
abline(h=0, col="red")

## Residual plot looks good. Thus the model is reasonable.




############ Q-Q plot  ##############

## Q-Q plot for Random effect

#ref.boat<-ranef(bs_lww_testLoc1.lmm10)[[1]]
#ref.var.boat<-tapply(residuals(bs_lww_testLoc1.lmm10,type="pearson", level=1), nonzerosamples1_testLoc1_new$boatID_id, var)
#qqnorm(ref.boat, main="Q-Q Normal-Boat Random Effects")
#qqline(ref.boat, col="red")



ref.boat<-ranef(bs_lww_testLoc1.lmm9)[[1]]
qqnorm(ranef(bs_lww_testLoc1.lmm9)$boatID$"(Intercept)", main="Q-Q plot for the random effect" )
qqline(ref.boat, col="red")








##### Q-Q plot for Residuals ################

qqnorm(residuals(bs_lww_testLoc1.lmm9, type="pearson"), main="Q-Q Normal-Residuals")
qqline(residuals(bs_lww_testLoc1.lmm9, type="pearson"), col="red")




























# the analyses look sound to me.  Please construct a table for each one
# that shows the test order and conclusions.  It would be useful to
# report the classical and the bootstrap p-value, and the term being
# tested in the table.  That way the reader can easily see the
# development of the model.

# Cheers,

# Andrew


















































