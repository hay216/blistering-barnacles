################################################################################
################################################################################
## Title: Imputation functions
## Author: Steve Lane
## Date: Wednesday, 29 March 2017
## Synopsis: Functions to run the censored regression imputation models
## Time-stamp: <2017-04-24 10:18:40 (slane)>
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Perform level 2 imputation.
################################################################################
################################################################################
## Function takes in the whole dataset, summarises within-boat wet weights by
## the median, performs one iteration of random forest imputation using mice,
## then passes back the imputed data.
lvl2Imp <- function(data){
    data <- data %>%
        mutate(wwImp = ifelse(wetWeight < 1.5, runif(n(), 0, 1.5), wetWeight))
    lvl2 <- data %>%
        group_by(boatID, days1, days2, midTrips, paintType, boatType,
                 ApproxHullSA) %>%
        summarise(m = median(wwImp)) %>% ungroup()
    pred.mat <- matrix(1 - diag(ncol(lvl2)), ncol(lvl2))
    pred.mat[1, ] <- pred.mat[, 1] <- 0
    miLvl2 <- mice(lvl2, m = 1, method = "rf", predictorMatrix = pred.mat,
                   maxit = 20, printFlag = FALSE)
    lvl2 <- complete(miLvl2) %>% select(-m)
    list(lvl2 = lvl2)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Create stan data function
################################################################################
################################################################################
createStanData <- function(lvl1, lvl2){
    ## Use model matrix to get indicators of what I need. These will be created
    ## alphabetically/numerically, just need to be careful feeding them in to
    ## stan.
    lvl1Mat <- model.matrix(
        ~ wetWeight + LocID + boatID + cens, data = lvl1) %>% as.data.frame()
    lvl2Mat <- model.matrix(
        ~ days1 + days2 + midTrips + ApproxHullSA + paintType + boatType +
            days1:paintType + days1:boatType + days2:paintType +
            days2:boatType + midTrips:paintType + midTrips:boatType +
            ApproxHullSA:paintType + ApproxHullSA:boatType, data = lvl2) %>%
        as.data.frame() %>% mutate(boatID = lvl2$boatID)
    ## Can pass list with too many variables to stan, so just add all as it
    ## makes it easier in the long run.
    obsData <- lvl1Mat %>% filter(cens == 0)
    censData <- lvl1Mat %>% filter(cens == 1)
    stanData <- c(with(obsData,
                       list(Y = wetWeight, N = nrow(obsData), numLoc = 3,
                            locID = cbind(LocIDKeel, LocIDRudder),
                            numBoat = max(boatID), boatID = boatID)),
                  with(censData,
                       list(nCens = nrow(censData), U = 1.5, 
                            locIDCens = cbind(LocIDKeel, LocIDRudder),
                            boatIDCens = boatID)),
                  with(lvl2Mat,
                       list(days1 = days1, days2 = days2, midTrips = midTrips,
                            hullSA = ApproxHullSA, numPaint = 3,
                            paintType = cbind(paintType2, paintType3),
                            numType = 3,
                            boatType = cbind(boatType2, boatType3),
                            days1paint = cbind(`days1:paintType2`,
                                               `days1:paintType3`),
                            days1boat = cbind(`days1:boatType2`,
                                              `days1:boatType3`),
                            days2paint = cbind(`days2:paintType2`,
                                               `days2:paintType3`),
                            days2boat = cbind(`days2:boatType2`,
                                              `days2:boatType3`),
                            midTripspaint = cbind(`midTrips:paintType2`,
                                                  `midTrips:paintType3`),
                            midTripsboat = cbind(`midTrips:boatType2`,
                                                 `midTrips:boatType3`),
                            ApproxHullSApaint =
                                cbind(`ApproxHullSA:paintType2`,
                                      `ApproxHullSA:paintType3`),
                            ApproxHullSAboat =
                                cbind(`ApproxHullSA:boatType2`,
                                      `ApproxHullSA:boatType3`)))
                  )
    ## Do a scaled version as well (centre and divide by sd)
    lvl2 <- lvl2 %>%
        mutate(
            days1 = as.numeric(scale(days1)),
            days2 = as.numeric(scale(days2)),
            midTrips = as.numeric(scale(midTrips)),
            ApproxHullSA = as.numeric(scale(ApproxHullSA))
        )
    lvl2Mat <- model.matrix(
        ~ days1 + days2 + midTrips + ApproxHullSA + paintType + boatType +
            days1:paintType + days1:boatType + days2:paintType +
            days2:boatType + midTrips:paintType + midTrips:boatType +
            ApproxHullSA:paintType + ApproxHullSA:boatType, data = lvl2) %>%
        as.data.frame() %>% mutate(boatID = lvl2$boatID)
    stanDataSc <- c(with(obsData,
                         list(Y = wetWeight, N = nrow(obsData), numLoc = 3,
                              locID = cbind(LocIDKeel, LocIDRudder),
                              numBoat = max(boatID), boatID = boatID)),
                    with(censData,
                         list(nCens = nrow(censData), U = 1.5, 
                              locIDCens = cbind(LocIDKeel, LocIDRudder),
                              boatIDCens = boatID)),
                    with(lvl2Mat,
                         list(days1 = days1, days2 = days2, midTrips = midTrips,
                              hullSA = ApproxHullSA, numPaint = 3,
                              paintType = cbind(paintType2, paintType3),
                              numType = 3,
                              boatType = cbind(boatType2, boatType3),
                              days1paint = cbind(`days1:paintType2`,
                                                 `days1:paintType3`),
                              days1boat = cbind(`days1:boatType2`,
                                                `days1:boatType3`),
                              days2paint = cbind(`days2:paintType2`,
                                                 `days2:paintType3`),
                              days2boat = cbind(`days2:boatType2`,
                                                `days2:boatType3`),
                              midTripspaint = cbind(`midTrips:paintType2`,
                                                    `midTrips:paintType3`),
                              midTripsboat = cbind(`midTrips:boatType2`,
                                                   `midTrips:boatType3`),
                              ApproxHullSApaint =
                                  cbind(`ApproxHullSA:paintType2`,
                                        `ApproxHullSA:paintType3`),
                              ApproxHullSAboat =
                                  cbind(`ApproxHullSA:boatType2`,
                                        `ApproxHullSA:boatType3`)))
                    )
    list(stanData = stanData, stanDataSc = stanDataSc)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Iterate between censored model and level 2 imputation.
################################################################################
################################################################################
## Function takes in the original dataset, an initial level 2 imputed/completed
## dataset, the regression model to be run (name of stan file), and the number
## of imputation iterations. Further arguments to rstan can be passed in via ...
iterImp <- function(origData, initData, modelname, impIter, ...){
    ## First, figure out for one iteration.
    ## Create local copies of datsets (that get overwritten)
    data <- origData
    impData <- initData
    ## Compile model
    model <- stan_model(paste0("../stan/", modelname, ".stan"))
    ## Create stan data (this may be a bit tricky when it comes to
    ## interactions... get it working first) But I can add all the possible
    ## interaction terms into the stanData list, and the program will only use
    ## what it needs to.
    for(i in seq_len(impIter)){
        obsData <- impData %>% filter(cens == 0)
        censData <- impData %>% filter(cens == 1)
        stanData <- c(with(obsData,
                           list(N = nrow(obsData), days1 = days1S,
                                days2 = days2S, midTrips = midTripsS,
                                numPaint = max(as.integer(paintTypeF)),
                                paintType = as.integer(paintTypeF),
                                numLoc = max(locIDInt),
                                boatID = boatIDInt, numBoat = max(boatIDInt),
                                locID = locIDInt, Y = wetWeight)),
                      with(censData,
                           list(nCens = nrow(censData), days1Cens = days1S,
                                days2Cens = days2S, midTripsCens = midTripsS,
                                paintTypeCens = as.integer(paintTypeF),
                                locIDCens = locIDInt, boatIDCens = boatIDInt,
                                U = 1.5)))
        ## Only need to save the yCens output for imputation
        hmcOut <- sampling(model, data = stanData, pars = "yCens",
                           iter = 500, chains = 4,
                           control = list(adapt_delta = 0.9))
        ## Extract the censored posterior predictions. For imputation, we want
        ## to make a draw from the posterior, so sample a draw from
        ## observation.
        yCens <- extract(hmcOut, "yCens")$yCens
        sampInd <- sample(seq_len(nrow(yCens)), ncol(yCens))
        yCens <- diag(yCens[sampInd, ])
        censData <- censData %>%
            mutate(wwImp = yCens,
                   wwLog = log(yCens))
        ## This dataset is not really needed, except to merge back onto the
        ## original dat so that we can iterate through.
        impCensData <- bind_rows(obsData, censData)
        ## Now merge censored predictions back on to the original data to feed
        ## into the level 2 imputation function
        data <- left_join(
            data %>% select(-wwImp, -wwLog),
            impCensData %>% select(boatIDInt, obsNum, wwImp, wwLog)
        )
        ## Now, pass into the level 2 imputation
        lvl2 <- lvl2Imp(data)
        ## And then store in the temp data frame to go back to the start...
        impData <- lvl2$impData
    }
    ## Return the fully completed/imputed data (as a list) to fit a final
    ## version of the model outside of this function (in stan). Also return the
    ## level 2 data to monitor convergence of the imputations (if required).
    obsData <- impData %>% filter(cens == 0)
    censData <- impData %>% filter(cens == 1)
    stanData <- c(with(obsData,
                       list(N = nrow(obsData), days1 = days1S,
                            days2 = days2S, midTrips = midTripsS,
                            numPaint = max(as.integer(paintTypeF)),
                            paintType = as.integer(paintTypeF),
                            numLoc = max(locIDInt),
                            boatID = boatIDInt, numBoat = max(boatIDInt),
                            locID = locIDInt, Y = wetWeight)),
                  with(censData,
                       list(nCens = nrow(censData), days1Cens = days1S,
                            days2Cens = days2S, midTripsCens = midTripsS,
                            paintTypeCens = as.integer(paintTypeF),
                            locIDCens = locIDInt, boatIDCens = boatIDInt,
                            U = 1.5)))
    list(stanData = stanData, lvl2Data = lvl2$lvl2)
}
################################################################################
################################################################################
