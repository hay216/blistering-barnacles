################################################################################
################################################################################
## Title: Imputation functions
## Author: Steve Lane
## Date: Wednesday, 29 March 2017
## Synopsis: Functions to run the censored regression imputation models
## Time-stamp: <2017-03-29 12:23:01 (slane)>
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
    lvl2 <- data %>% group_by(days1S, days2S, midTripsS, paintTypeF) %>%
        summarise(m = median(wwLog)) %>% ungroup()
    miLvl2 <- mice(lvl2, m = 1, method = "rf", maxit = 1)
    lvl2 <- complete(miLvl2) %>% mutate(boatIDInt = 1:n())
    impData <- left_join(
        data %>% select(-days1S, -days2S, -midTripsS, -paintTypeF),
        lvl2 %>% select(-m)
    )
    list(impData = impData, lvl2 = lvl2)
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
