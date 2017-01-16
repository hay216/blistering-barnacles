#This code analyses the EV II data set
#The code is going to read the excel spreadsheets directly using the RODBC library
#this library is required, along with MASS
library(RODBC)
library(stats)
library(ellipse)
library (MASS)
library(scatterplot3d)

#PART 1: READING THE DATA
#Set the working directory so that it reads to the correct path
setwd("H:/Projects/DEH Empirical Validation II/Reports")

#Make a connection to the excel workbook - note that the spreadsheet names must not include any spaces #otherwise the connection cannot read them
evIIData <-odbcConnectExcel("H:/Projects/DEH Empirical Validation II/Survey/Survey results - Oct 05.xls")

#Check what tables are held in the connected file
sqlTables(evIIData)

#Make four dataFrames using the sqlFetch function
samples <- sqlFetch(evIIData, "Samples")
activity <- sqlFetch(evIIData, "Activity")
vessel <- sqlFetch(evIIData, "Vessel")
plankton <- sqlFetch(evIIData, "Plankton")

#Beware that R may read empty rows so constrain the dataset!
samples <- samples[!is.na(samples$boatID), ]
activity <- activity[!is.na(activity$vesID), ]
vessel <- vessel[!is.na(vessel$boatID), ]
plankton <- plankton[!is.na(plankton$boatID), ]

#****************************************************************************************************************
#PART 2A: INITIAL DATA EXPLORATION - vessel
#Histogram and summary statistics of the vessel length and hull surface area
attach(vessel)
names(vessel)

bitmap("Vessel_length_hist.tif", type="tiffpack", width = 8, height = 8, res =600)
#windows()
par(mfrow=c(2,2))
par(xpd = NA)
hist(vessel$Length, prob = F, xlab = "Vessel length (m)", ylab = "Frequency", main = "")
text(5, 27, "a")
plot(density(vessel$Length), xlab = "Vessel length (m)", ylab = "Density", main = "" )
text(2, 0.13, "b")
hist(vessel$ApproxHullSA, prob = F, xlab = "Approx. hull surface area (m2)", ylab = "Frequency", main = "")
text(0, 43, "c")
plot(density(vessel$ApproxHullSA), xlab = "Approx. hull surface area (m2)", ylab = "Density", main = "" )
text(0, 0.022, "d")
dev.off()

#Aggregate the data by boat type (similar to tapply function) to see length and hull surface area statistics,
#use tapply for variance calculation, create different datagroups to investigate these statistics for similar
#vessels - note the exclusion of vessel #22 which was an outlier
#aggregate(vessel, list(VesselType = vessel$boatType), mean)
#tapply(vessel$ApproxHullSA, list(vessel$boatType), var)
g1 <- vessel[boatType %in% c("Motor cruiser", "Yacht"), ]
g2 <- vessel[boatType %in% c("Fishing vessel (Abalone mothership)", "Fishing vessel(Lobster/scallop)", "Fishing vessel (Long line)"), ]

attach(g2)
g2mod <- g2[boatID !=22, ]
g3 <- merge(g1, g2mod, all.x = T, all.y = T)
mean(g3$ApproxHullSA)
var(g3$ApproxHullSA)

#PART 2B: INITIAL DATA EXPLORATION - vessel activity patterns
#Read in the sample date (sDate), the last trip date (ltDate), the last slip date (lsDate) and the last anti-foul
#date (afDate).  Note use of the coercive "as" function because R read the excel dates as a factor
attach(activity)
names(activity)
sDate <- strptime(as(activity$sampleDate, "character"), "%Y-%m-%d")
ltDate <- strptime(as(activity$lastTripDate, "character"), "%Y-%m-%d")
lsDate <- strptime(as(activity$lastSlipDate, "character"), "%Y-%m-%d")
afDate <- strptime(as(activity$lastAFDate, "character"), "%Y-%m-%d")
iwDate <- strptime(as(activity$lastINWDate, "character"), "%Y-%m-%d")

#Remove the NA's from the last slip, in water clean and anti-foul date replace with an earlier date so long as
#there is no overlap
lsDate[is.na(lsDate)] <- strptime("1980-01-01", "%Y-%m-%d")
afDate[is.na(afDate)] <- strptime("1980-01-01", "%Y-%m-%d")
iwDate[is.na(iwDate)] <- strptime("1980-01-01", "%Y-%m-%d")

#Set up the covariate data frame which will also serve for the statistical modelling, enter the time between
#the sample date and the last trip date (days1) in days, the time between the sample date and the date the #vessel was last cleaned, slipped or antifouled (days2) in days, and the median number of trips per annum, #together with boat type and code
covars <- data.frame(array(NA, dim = c(nrow(activity), 8)))
names(covars) <- c("boatID", "days1", "days2", "midTrips", "boatType", "boatCode", "paintType", "paintRating")
covars$boatID <- activity$vesID
covars$boatType <- activity$boatType
covars$boatCode <- activity$boatCode
covars$paintType <- activity$Type
covars$paintRating <- activity$PerfRating
covars$days1 <- as.numeric(round((sDate - ltDate)/(60*60*24), 0))

#Assign days2 the same elements and class as afDate - note the same class is important for subsequent manipulations
covars$days2 <- afDate

#Set up a loops that find the latest date since the vessel was last cleaned, slipped or antifouled and thereby #calculate days2
 for (i in 1: nrow(covars)){
	if(lsDate[i] > afDate[i]) covars$days2[i] <- lsDate[i]
	if(iwDate[i] > covars$days2[i]) covars$days2[i] <- iwDate[i]
}
covars$days2 <- as.numeric(round((sDate - covars$days2), 0))

#Calculate median trip activity per annum
covars$midTrips <- activity$minTrips + (activity$maxTrips - activity$minTrips)/2
#table(covars$midTrips)

#Calculate the number of last trips less than or equal to 14 and 50 days and the number of vessels last
#cleaned within 1/2 and 1 years
attach(covars)
countDay14 <- covars[!is.na(days1) & days1 <= 14, 2]
countDay50 <- covars[!is.na(days1) & days1 <= 50, 2]
length(countDay14)
length(countDay50)
countYr1 <- covars[days2 <= 178, 3]
countYr2 <- covars[days2 <= 356, 3]
length(countYr1)
length(countYr2)

#Plot a histogram of the number of days between the sample date and the last trip date and last slip clean,
#in water clean or antifoul date, and the median trip activity
bitmap("Activity_hist.tif", type="tiffpack", width = 8, height = 8, res =600)
#windows()
par(mfrow=c(2,2))
par(xpd = NA)
hist(covars$days1, prob = F, xlab = "Days since last trip", ylab = "Frequency", main = "")
text(0, 40, "a")
hist(covars$days2, prob = F, xlab = "Days since last cleaned or antifouled", ylab = "Frequency", main = "")
text(0, 25, "b")
hist(covars$midTrips, prob = F, xlab = "Median number of trips per annum", ylab = "Frequency", main = "")
text(0, 40, "c")
dev.off()

#PART 2C: INITIAL DATA EXPLORATION - sample wet weight
#Name and select the data for graphical purposes - note the joining of the
#deck samples and fishing gear samples and the ordering of the samples
attach(samples)
names(samples)
wetHull <- samples[Location %in% c("Wet Hull"),  ]
wetHull <- wetHull[order(wetHull$wetWeight,decreasing = T), ]
intSam <- samples[Location %in% c("Internal spaces"), ]
intSam <- intSam[order(intSam$wetWeight, decreasing = T), ]
propSam <-samples[Location %in% c("Propeller"), ]
propSam <- propSam[order(propSam$wetWeight, decreasing = T), ]
rudSam <-samples[Location %in% c("Rudder"), ]
rudSam <- rudSam[order(rudSam$wetWeight, decreasing = T), ]
deckSam <- samples[Location %in% c("Fishing gear", "Deck", "Anchor"), ]
deckSam <- deckSam[order(deckSam$wetWeight, decreasing = T), ]

#Size frequency plot of the original data and then the ln transformed data
bitmap("Wet_weight_hist.tif", type="tiffpack", width = 8, height = 8, res =600)
#windows()
par(mfrow=c(2,2))
hist(wetHull$wetWeight, prob=F, xlab = "Wet weight (g)", ylab = "Count", main = " Wet Hull")
hist(intSam$wetWeight, prob = F, xlab = "Wet weight (g)", ylab = "Count", main = "Internal spaces")
hist(propSam$wetWeight, prob = F, xlab = "Wet weight (g)", ylab = "Count", main = "Propeller")
hist(rudSam$wetWeight, prob = F, xlab = "Wet weight (g)", ylab = "Count", main = "Rudder")
dev.off()

bitmap("Log_wet_weight_hist.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
par(mfrow=c(2,2))
hist(log(wetHull$wetWeight + 1), prob=F, xlab = "Ln(Wet weight + 1) (g)", ylab = "Count", main = " Wet Hull")
hist(log(intSam$wetWeight + 1), prob = F, xlab = "Ln(Wet weight + 1) (g)", ylab = "Count", main = "Internal spaces")
hist(log(propSam$wetWeight + 1), prob = F,xlab = "Ln(Wet weight + 1) (g)", ylab = "Count", main = "Propeller")
hist(log(rudSam$wetWeight + 1), prob = F, xlab = "Ln(Wet weight + 1) (g)", ylab = "Count", main = "Rudder")
dev.off()

#Create three data frames for the zero, non-zero data and all the data, note the use of aggregate functions
#rather than tapply since the former creates a data frame.  Note also the use of the factor command to
#exclude factors that only had zeroes, (it is otherwise redundant because R already interprets LocID as a factor),
#the use of the unique command to identify the individual factors that the Shapiro test is eventually
#applied to using the == function, and the storage of the results in a data frame (rather than a matrix)
#because this allows different classes of information to be stored (i.e. factor LocID and numeric Shapiro test statistics)
data <- data.frame(LocID = samples$LocID, wetWeight = samples$wetWeight)
data1 <- data[data$wetWeight > 0.5, ]
data0 <- data[data$wetWeight <= 0.5, ]
#data$LocID <- factor(data$LocID, exclude=NULL)
data1$LocID <- factor(data1$LocID, exclude=NULL)
data0$LocID <- factor(data0$LocID, exclude=NULL)
length.data <- aggregate(data$wetWeight, by = list(LocID = data$LocID), length)
length.data1 <- aggregate(data1$wetWeight, by = list(LocID = data1$LocID), length)
length.data0 <- aggregate(data0$wetWeight, by = list(LocID = data0$LocID), length)
unique.loc1 <- unique(data1$LocID)
unique.loc0 <- unique(data0$LocID)
#unique.loc <- unique(data$LocID)
names(length.data0)[[2]] <- "n0"
names(length.data1)[[2]] <- "n1"
names(length.data)[[2]] <- "n"

#Size frequency plot by location the data, note the use of the browser command to step through the loop
#and the need to specify LocID as a factor and the use of the unique command for the data0 and data1 datesets
#(because R holds over the LocID factors that were taken out to create these datasets)
bitmap("Figure2D4.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 1:12) {
	if(length.data$n[i] > 1) {
#browser()
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]) , xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

bitmap("Figure2D5.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 13:24) {
	if(length.data$n[i] > 1) {
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]), xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

bitmap("Figure2D6.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 25:36) {
	if(length.data$n[i] > 1) {
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]), xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

bitmap("Figure2D7.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 37:48) {
	if(length.data$n[i] > 1) {
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]), xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

bitmap("Figure2D8.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 49:50) {
	if(length.data$n[i] > 1) {
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]), xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

bitmap("Figure2D9.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 51:nrow(length.data)) {
	if(length.data$n[i] > 1) {
hist(log(data$wetWeight[data$LocID==length.data$LocID[i]]+1), main = as.character(length.data$LocID[i]), xlab = "Ln (Wet weight + 1) (g)")
	}
}
dev.off()

#Size frequency plot by location the non-zero data, note the use of the browser command to step through the loop
bitmap("Figure2D10.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 1:12) {
	if(length.data1$n1[i] > 1) {
#browser()
hist(log(data1$wetWeight[data1$LocID==length.data1$LocID[i]]), main = as.character(length.data1$LocID[i]), xlab = "Ln (Wet weight) (g)")
	}
}
dev.off()

bitmap("Figure2D11.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 13:24) {
	if(length.data1$n1[i] > 1) {
hist(log(data1$wetWeight[data1$LocID==length.data1$LocID[i]]), main = as.character(length.data1$LocID[i]), xlab = "Ln (Wet weight) (g)")
	}
}
dev.off()

bitmap("Figure2D12.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 25:nrow(length.data1)) {
	if(length.data1$n1[i] > 1) {
hist(log(data1$wetWeight[data1$LocID==length.data1$LocID[i]]), main = as.character(length.data1$LocID[i]), xlab = "Ln (Wet weight) (g)")
	}
}
dev.off()

#Box plot the non-zero data, rotating and setting the axis font to 80% of the default (see #?par)
bitmap("Figure2D13.tif", type = "tiffpack", width = 8, height = 8, res=600)
par(mfrow=c(2,2))
wetHull <- wetHull[wetHull$wetWeight > 0.5, ]
attach(wetHull)
plot(LocID, log(wetWeight), las = 2, xlab = "Location ID", ylab = "Ln (Wet weight) (g)", main = "Wet Hull")
intSam <- intSam[intSam$wetWeight > 0.5, ]
attach(intSam)
plot(LocID, log(wetWeight), las = 2, xlab = "Location ID", ylab = "Ln (Wet weight) (g)", main = "Internal spaces")
propSam <- propSam[propSam$wetWeight > 0.5, ]
attach(propSam)
plot(LocID, log(wetWeight), las = 2, xlab = "Location ID", ylab = "Ln (Wet weight) (g)", main = "Propeller")
rudSam <- rudSam[rudSam$wetWeight > 0.5, ]
attach(rudSam)
plot(LocID, log(wetWeight), las = 2, xlab = "Location ID", ylab = "Ln (Wet weight) (g)", main = "Rudder")
dev.off()

#Q-Q plot the non-zero data by Location ID
bitmap("Figure2D14.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 1:12) {
	if(length.data1$n[i] > 3) {
		qqnorm(log(data1$wetWeight[data1$LocID==unique.loc1[i]]), main = as.character(unique.loc1[i]))
		qqline(log(data1$wetWeight[data1$LocID==unique.loc1[i]]))
	}
}
dev.off()

bitmap("Figure2D15.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 13:24) {
	if(length.data1$n[i] > 3) {
		qqnorm(log(data1$wetWeight[data1$LocID==unique.loc1[i]]), main = as.character(unique.loc1[i]))
		qqline(log(data1$wetWeight[data1$LocID==unique.loc1[i]]))
	}
}
dev.off()

bitmap("Figure2D16.tif", type = "tiffpack", width = 8, height = 10, res=600)
par(mfrow=c(4,3))
for (i in 25:nrow(length.data1)) {
	if(length.data1$n[i] > 3) {
		qqnorm(log(data1$wetWeight[data1$LocID==unique.loc1[i]]), main = as.character(unique.loc1[i]))
		qqline(log(data1$wetWeight[data1$LocID==unique.loc1[i]]))
	}
}
dev.off()

#Perform a Shapiro-Wilks test on the non-zero data (the aggregate function rather than tapply above allows the use of merge functions by column name to construct the output table) and construct an output  table
#showing the sample size (n) number of zero's (n0), the number of non-zeros (n1), and the Shapiro test result
shapTest <- data.frame(array(NA, dim = c(nrow(length.data1),3)))
names(shapTest) <- c("LocID","w","p")
for (i in 1:nrow(length.data1)) {
	shapTest[i,1] <- as.character(unique.loc1[i])
	if(length.data1$n1[i] > 3) {
		shapTest[i,2] <- shapiro.test(log(data1$wetWeight[data1$LocID==length.data1$LocID[i]]))[[1]]
		shapTest[i,3] <- shapiro.test(log(data1$wetWeight[data1$LocID==length.data1$LocID[i]]))[[2]]
	}
}
test <- merge(length.data, length.data1, by.x = "LocID", by.y = "LocID", all.x = T, all.y = T)
test <- merge(test, length.data0, by.x = "LocID", by.y = "LocID", all.x = T, all.y = T)
test <- merge(test, shapTest, by.x = "LocID", by.y = "LocID", all.x = T, all.y = T)
write.table(test, file = "TableYY.csv", sep = ",")

#Construct a table showing the mean and variance of the non-zero data and merge it with the test table, note the use of tapply here because we are reading the results directly into a data frame - beware when reading
#this information out of the data frame because the factor information is still retained with it!
deltaData <- data.frame(array(NA, dim = c(nrow(length.data1), 5)))
names(deltaData) <- c("LocID", "muHat", "varHat", "deltaMu", "deltaVar")
deltaData$LocID <- unique.loc1
deltaData$muHat <- tapply(log(data1$wetWeight), data1$LocID, mean)
deltaData$varHat <- tapply(log(data1$wetWeight), data1$LocID, var)
delta <- merge(test, deltaData, by.x = "LocID", by.y = "LocID", all.x = T, all.y = T)

#Calculate the mean and variance of the delta distribution
for (i in 1:nrow(delta)) {
	if(!is.na(delta$w[i])) {
		dHat <- delta$n0[i]/delta$n[i]
		#bFunc <- besselJ(0.5*delta$varHat[[i]], delta$n1[i])
		#delta$deltaMu[i] <- (1-dHat)*exp(delta$muHat[[i]])*bFunc
		delta$deltaMu[i] <- (1-dHat)*exp(delta$muHat[[i]] + delta$varHat[[i]])
		delta$deltaVar[i] <- (1/delta$n[i])*exp(2*delta$muHat[[i]] + delta$varHat[[i]])*(dHat*(1-dHat) + 0.5*(1-dHat)*(2* delta$varHat[[i]] + (delta$varHat[[i]]^2)))
	}
}
write.table(delta, file = "TableD.csv", sep = ",")

#*****************************************************************************************************************
#PART 3A: STATISTICAL MODELLING ZERO DATA
#Examine the relationships between proportion of clean samples and time since the vessel was last used (days1), time since the vessel was last cleaned (days2), number of trips per annum (midTrips) and paint #type
attach(samples)
samSize1 <- aggregate(samples$wetWeight, by = list(boatID = samples$boatID, LocID = samples$LocID), length)
small.numbers<-function(x){
	return(length(x[x<=0.5]))
}
samSize2 <- aggregate(samples$wetWeight, by = list(boatID = samples$boatID, LocID = samples$LocID), small.numbers)
samSize1 <- cbind(samSize1, samSize2[ , 3])
dimnames(samSize1)[[2]][3] <- "n"
dimnames(samSize1)[[2]][4] <- "n0"

attach(samSize1)
testLoc <- c("IB", "HA", "PJ", "PB", "HD", "HAH", "HP", "HH", "HK")

bitmap("Figure3A1.tif", type = "tiffpack", width = 8, height = 10, res=600)
#windows( )
par(mfrow=c(4,3))
for (i in 1:length(testLoc)){
	dum <- merge(covars, samSize1[samSize1$LocID==testLoc[i], ], by.x = "boatID", by.y = "boatID")
	par(xpd = NA)
	plot(dum$days1, dum$n0/dum$n, main = testLoc[i], xlab = "Days since vessel last used", ylab = "Proportion of clean samples", cex = 1, col = "red", pch = 20)
	text(dum$days1, dum$n0/dum$n, labels = dum$boatID, pos = 3, cex = 0.6)
}
dev.off()

bitmap("Figure3A2.tif", type = "tiffpack", width = 8, height = 10, res=600)
#windows( )
par(mfrow=c(4,3))
for (i in 1:length(testLoc)){
	dum <- merge(covars, samSize1[samSize1$LocID==testLoc[i], ], by.x = "boatID", by.y = "boatID")
	par(xpd = NA)
	plot(dum$days2, dum$n0/dum$n, main = testLoc[i], xlab = "Days since vessel last cleaned", ylab = "Proportion of clean samples", cex = 1, col = "red", pch = 20)
	text(dum$days2, dum$n0/dum$n, labels = dum$boatID, pos = 3, cex = 0.6)
}
dev.off()

bitmap("Figure3A3.tif", type = "tiffpack", width = 8, height = 10, res=600)
#windows( )
par(mfrow=c(4,3))
for (i in 1:length(testLoc)){
	dum <- merge(covars, samSize1[samSize1$LocID==testLoc[i], ], by.x = "boatID", by.y = "boatID")
	par(xpd = NA)
	plot(dum$midTrips, dum$n0/dum$n, main = testLoc[i], xlab = "Median trips per annum", ylab = "Proportion of clean samples", cex = 1, col = "red", pch = 20)
	text(dum$midTrips, dum$n0/dum$n, labels = dum$boatID, pos = 3, cex = 0.6)
}
dev.off()

#Quasi-bionomial GLM model for the proportion of zero samples, in IB compare the two models - one weighted by the number of samples one without
dumIB <- merge(covars, samSize1[samSize1$LocID==testLoc[1], ], by.x = "boatID", by.y = "boatID")
modIB1 <- glm(n0/n ~ days2, data = dumIB, family = quasibinomial, method = "glm.fit", trace = F)
summary(modIB1)
bitmap("Figure3A4.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
par(mfrow=c(2,2))
plot(modIB1, main = testLoc[1])
dev.off()

modIB2 <- glm(n0/n ~ days2, data = dumIB, weight = n, family = quasibinomial, method = "glm.fit", trace = F)
summary(modIB2)
windows()
par(mfrow=c(2,2))
plot(modIB2, main = testLoc[1])

#GLM model for the proportion of zero samples, in the paddle wheel and booth (HK) compare the two models -
#one weighted by the number of samples one without
dumHK <- merge(covars, samSize1[samSize1$LocID==testLoc[9], ], by.x = "boatID", by.y = "boatID")
modHK1 <- glm(n0/n ~ days2, data = dumHK, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHK1)
bitmap("Figure3A5.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
par(mfrow=c(2,2))
plot(modHK1, main = testLoc[9])
dev.off()

modHK2 <- glm(n0/n ~ days2, data = dumHK, weight = n, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHK2)
windows()
par(mfrow=c(2,2))
plot(modHK2, main = testLoc[9])

#GLM model for the proportion of zero samples, in the cover plates (HD) compare the two models - one weighted by the number of samples one without
dumHD <- merge(covars, samSize1[samSize1$LocID==testLoc[5], ], by.x = "boatID", by.y = "boatID")
modHD1 <- glm(n0/n ~ days2, data = dumHD, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHD1)
bitmap("Figure3A6.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
par(mfrow=c(2,2))
plot(modHD1, main = testLoc[5])
dev.off()

modHD2 <- glm(n0/n ~ days2, data = dumHD, weight = n, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHD2)
windows()
par(mfrow=c(2,2))
plot(modHD2, main = testLoc[5])

#GLM model for the proportion of zero samples, in the depth sounder booth(HH) compare the two models - one weighted by the number of samples one without
dumHH <- merge(covars, samSize1[samSize1$LocID==testLoc[8], ], by.x = "boatID", by.y = "boatID")
modHH1 <- glm(n0/n ~ days2, data = dumHH, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHH1)
bitmap("Figure3A7.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
par(mfrow=c(2,2))
plot(modHH1, main = testLoc[5])
dev.off()

modHH2 <- glm(n0/n ~ days2, data = dumHH, weight = n, family = quasibinomial, method = "glm.fit", trace = F)
summary(modHH2)
windows()
par(mfrow=c(2,2))
plot(modHH2, main = testLoc[5])

#Plot up the predicted values from these GLM models
#windows()
bitmap("Zero_samples_GLM.tif", type="tiffpack", width = 8, height = 8, res = 600)
par(mfrow = c(2,2))
newx <- seq(0, 1000, len = 500)
newy <- predict(modIB1, newdata = data.frame(days2 = newx), type = "link", se.fit = T, conf.level = 0.95)
plot(dumIB$days2, dumIB$n0/dumIB$n, xlab="Days since vessel last cleaned", ylab = "Proportion of clean samples", main = testLoc[1])
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]] + 2*newy$se.fit), lty = 2)
lines(newx, quasibinomial()$linkinv(newy[[1]] - 2*newy$se.fit), lty = 2)

newy <- predict(modHK1, newdata = data.frame(days2 = newx), type = "link", se.fit = T, conf.level = 0.95)
plot(dumHK$days2, dumHK$n0/dumHK$n, xlab="Days since vessel last cleaned", ylab = "Proportion of clean samples", main = testLoc[9])
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]] + 2*newy$se.fit), lty = 2)
lines(newx, quasibinomial()$linkinv(newy[[1]] - 2*newy$se.fit), lty = 2)

newy <- predict(modHD1, newdata = data.frame(days2 = newx), type = "link", se.fit = T, conf.level = 0.95)
plot(dumHD$days2, dumHD$n0/dumHD$n, xlab="Days since vessel last cleaned", ylab = "Proportion of clean samples", main = testLoc[5])
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]] + 2*newy$se.fit), lty = 2)
lines(newx, quasibinomial()$linkinv(newy[[1]] - 2*newy$se.fit), lty = 2)

newy <- predict(modHH1, newdata = data.frame(days2 = newx), type = "link", se.fit = T, conf.level = 0.95)
plot(dumHH$days2, dumHH$n0/dumHH$n, xlab="Days since vessel last cleaned", ylab = "Proportion of clean samples", main = testLoc[8])
lines(newx, quasibinomial()$linkinv(newy[[1]]))
lines(newx, quasibinomial()$linkinv(newy[[1]] + 2*newy$se.fit), lty = 2)
lines(newx, quasibinomial()$linkinv(newy[[1]] - 2*newy$se.fit), lty = 2)
dev.off()

#************************************************************************************************************************
#PART 3B: ANALYSIS OF COVARIANCE for the hull quadrats using days2 and paintType
testLoc <- c("IB", "HA", "PJ", "PB", "HD", "HAH", "HP", "HH", "HK")
tempData <- samples[samples$wetWeight > 0.5, ]
tempData <- tempData[ , c("boatID", "LocID", "wetWeight")]
muData <- aggregate(log(tempData$wetWeight), by = list(boatID = tempData$boatID, LocID = tempData$LocID), mean)
modData<- merge(covars, muData[muData$LocID==testLoc[2], ], by.x = "boatID", by.y = "boatID")
modData$paintRating <- relevel(modData$paintRating, ref = "Very Poor")

#Output the raw data table for the report
write.table(modData, file= "Covariance_Data.csv", append = F, sep = ",")

#Plot the raw data
#windows()
bitmap("Covariance_data_plot1.tif", type="tiffpack", width = 8, height = 8, res = 600)
par(mfrow = c(2,1))
plot(modData$midTrips, modData$x, type = "n", xlab = "Median trips per annum", ylab = "Log10(mean hull quadrat wet weight)", cex = 0.6)
text(modData$midTrips, modData$x, substring(as.character(modData$boatID),1,2), cex = 0.6)
plot(modData$days2, modData$x, type = "n", xlab = "Days since last cleaned or antifouled", ylab = "Log10(mean hull quadrat wet weight)", cex = 0.6)
text(modData$days2, modData$x, substring(as.character(modData$boatID),1,2), cex = 0.6)
dev.off()

#windows()
bitmap("Covariance_data_plot2.tif", type="tiffpack", width = 8, height = 8, res = 600)
par(mfrow = c(2,1))
plot(modData$midTrips, modData$x, type = "n", xlab = "Median trips per annum", ylab = "Log10(mean hull quadrat wet weight)", cex = 0.6)
text(modData$midTrips, modData$x, substring(as.character(modData$paintType),1,1), cex = 0.6)
plot(modData$days2, modData$x, type = "n", xlab = "Days since last cleaned or antifouled", ylab = "Log10(mean hull quadrat wet weight)", cex = 0.6)
text(modData$days2, modData$x, substring(as.character(modData$paintType),1,1), cex = 0.6)
dev.off()

#Final model which is used with days2 and midTrips*paintType interaction
coMod6 <- lm(x ~ days2 + midTrips*paintType, modData)
summary(coMod6)
anova(coMod6)

#Basic model diagnostics
bitmap("Covariance_diagnostics.tif", type ="tiffpack", width = 8, height = 8, res =600)
#windows()
par(mfrow=c(2,2))
plot(coMod6)
dev.off()

#Create a table with the model diagnostics
temp <- modData[!is.na(modData$days2),]
temp <- temp[!is.na(temp$midTrips),]
temp <- temp[!is.na(temp$paintType),]
diag_tab <- temp$boatID
diag_tab <- cbind(diag_tab, round(coMod6$res, 2))
diag_tab <- cbind(diag_tab, round(rstandard(coMod6),2))
diag_tab <- cbind(diag_tab, round(rstudent(coMod6),2))
diag_tab <- cbind(diag_tab, round(hatvalues(coMod6),3))
diag_tab <- cbind(diag_tab, round(cooks.distance(coMod6),5))
write.table(diag_tab, file = "Covariance_diagnostics.csv", append = F, sep = ",")

#Model without apparent outliers
outliers <- c(21, 46)
temp <- subset(modData,boatID!= 24 & boatID!=51)
coMod6b <- lm(x ~ days2 + midTrips*paintType, data = temp)
summary(coMod6b)
anova(coMod6b)

#****************************************************************************************************
#Part 4A: IMEA analysis
#Read in the IMEA results and make a dataFrame using the sqlFetch function
imea <-odbcConnectExcel("H:/Projects/DEH Empirical Validation II/Reports/IMEA1.xls")
sqlTables(imea)
imea.data<-sqlFetch(imea, "Data")
imea.data <- imea.data[!is.na(imea.data$LocID), ]

#Create a data frame by location id with the mean of the data for comparison with the IMEA analysis
muData1 <- aggregate(log(samples$wetWeight +1), by = list(LocID = samples$LocID), mean)
varData1 <- aggregate(log(samples$wetWeight +1), by = list(LocID = samples$LocID), var)
countData <- aggregate(samples$wetWeight, by = list(LocID = samples$LocID), length)
tempData <- merge(muData1, varData1, by.x = "LocID", by.y="LocID")
tempData <- merge(tempData, countData, by.y = "LocID")
names(tempData)[[2]] <- "mu"
names(tempData)[[3]] <- "var"
names(tempData)[[4]] <- "n"
tempData <- merge(tempData, imea.data, by.x = "LocID", by.y = "LocID", all.x = T, all.y = T)
tempData$mu <- round(tempData$mu, 2)
tempData$var <- round(tempData$var, 2)
tempData$AvSevRat <- round(tempData$AvSevRat, 2)
tempData$AvOccRat <- round(tempData$AvOccRat, 2)
tempData$SOR <- round(tempData$SOR, 2)
tempData <- tempData[!is.na(tempData$mu), ]
tempData <- tempData[!is.na(tempData$SOR), ]
write.table(tempData, file= "H:/Projects/DEH Empirical Validation II/Reports/TableZZ.csv", append = F, sep = ",")

#Plot the IMEA score by mean wet weight for the different locations
attach(tempData)
bitmap("IMEA_scatterplot.tif", type = "tiffpack", width = 8, height = 8, res=600)
#windows()
plot(tempData$SOR, tempData$mu, ylab = "Ln (Mean wet weight) (g)", xlab = "IMEA score", pch = 19, cex = 0.5)
text(tempData$SOR, tempData$mu, labels = tempData$LocID, pos = 4, cex = 0.7)
dev.off()

#Perform a linear model on the IMEA score
mod8 <- lm(mu ~ SOR, data = tempData)
summary(mod8)
bitmap("IMEA_reg_diag.tif", type = "tiffpack", width = 8, height = 8, res=600)
par(mfrow=c(2,2))
plot(mod8)
dev.off()

#****************************************************************************************
#PART 5: Undaria fouling analysis
#Construct two new data frames showing the proportion of infected samples by location #category
und1.sum <- aggregate(samples$Und, by = list(boatID = samples$boatID, LocCode = samples$LocCode), sum, na.rm = T)
und1.length <- aggregate(samples$Und, by = list(boatID = samples$boatID, LocCode = samples$LocCode), length)
und1.prop <- und1.sum$x/und1.length$x
und1.sum <- data.frame(und1.sum, und1.prop)
und1.sum$LocCode <- as.numeric(und1.sum$LocCode)

und2.sum <- aggregate(plankton$Und, by = list(boatID = plankton$boatID), sum, na.rm = T)
und2.length <- aggregate(plankton$Und, by = list(boatID = plankton$boatID), length)
und2.prop <- und2.sum$x/und2.length$x

attach(vessel)
sam.loc <- vessel[c("samCode", "boatID")]
und1.sum <- merge(sam.loc, und1.sum, by.x = "boatID")

attach(und1.sum)
u1 <- und1.sum[LocCode ==1, ]
#a1$samCode <- as.numeric(a1$samCode)
u2 <- und1.sum[LocCode ==2, ]
#a2$samCode <- as.numeric(a2$samCode)
u3 <- und1.sum[LocCode ==3, ]
#a3$samCode <- as.numeric(a3$samCode)

bitmap("Undaria_positive_samples.tif", type = "tiffpack", width = 10, height = 10, res=600)
#windows()
par(mfrow=c(2,2))
plot(u1$boatID, u1$und1.prop, pch = u1$samCode, main = "Hull", xlim = c(30, 54), ylim = c(0,1), xlab = "", ylab = "")
plot(u2$boatID, u2$und1.prop, pch = u2$samCode, main = "Propeller rudder & anchor", xlim = c(30, 54), ylim = c(0,1), xlab = "", ylab = "")
plot(u3$boatID, u3$und1.prop, pch = u3$samCode, main = "Internal spaces", xlim = c(30, 54), ylim = c(0,1), xlab = "Boat reference", ylab = "Proportion of Undaria positive samples")
plot(u1$boatID, und2.prop, pch = u2$samCode, main = "Plankton samples", xlim = c(30, 54), ylim = c(0,1), xlab = "", ylab = "")
par(xpd = NA)
legend(18,1.5,  c("Sandringham yacht club", "Royal Yacht Club of Victoria", "Hobsons Bay Yacht Club"), pch=c(3,4,5), cex=0.8)
dev.off()

#Make a new data frame showing the proportion of positive samples by location and print this out as a table
Und3.sum <- aggregate(samples$Und, by = list(boatID = samples$boatID, LocID = samples$LocID), sum, na.rm = T)
Und3.length <- aggregate(samples$Und, by = list(boatID = samples$boatID, LocID = samples$LocID), length)

Und3.prop <- Und3.sum$x/Und3.length$x
Und3.sum <- data.frame(Und3.sum, Und3.length, Und3.prop)
Und3.sum <- Und3.sum[Und3.prop !=0, ]
write.table(Und3.sum, file= "Undaria_table.csv", append = F, sep = ",")

#Save the R workspace
save.image("H:/Projects/DEH Empirical Validation II/Reports/Data analysis - Nov 06.RData")

#************************************************************
#Models which were tried but not used in the analysis
#coMod1 <- lm(x ~ days2*paintType, modData)
#summary(coMod1)
#anova(coMod1)

#coMod2 <- lm(x ~ days2*paintRating, modData)
#summary(coMod2)
#anova(coMod2)

#coMod3 <- lm(x ~ days2 + paintRating, modData)
#summary(coMod3)

#coMod4 <- lm(x ~ midTrips + paintType, modData)
#summary(coMod4)

#coMod5 <- lm(x ~ midTrips*paintType, modData)
#summary(coMod5)

#coMod7 <- lm(x ~ days2*midTrips*paintType, modData)
#summary(coMod7)
#anova(coMod7)
