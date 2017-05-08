################################################################################
################################################################################
## Title: Imputation functions
## Author: Steve Lane
## Date: Wednesday, 29 March 2017
## Synopsis: Functions to run the censored regression imputation models
## Time-stamp: <2017-05-08 08:16:32 (slane)>
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Function to load packages (and install if required).
################################################################################
################################################################################
ipak <- function(pkg){
    ## Check for github packages (throw away github username)
    chk.git <- gsub(".*/", "", pkg)    
    new.pkg <- pkg[!(chk.git %in% installed.packages()[, "Package"])]
    if(!(length(new.pkg) == 0)){
        git.ind <- grep("/", new.pkg)
        if(length(git.ind) == 0){
            install.packages(new.pkg, dependencies = TRUE,
                             repos = "https://cran.csiro.au/")
        } else {
            devtools::install_github(new.pkg[git.ind])
        }
    }
    sapply(chk.git, require, character.only = TRUE)
}
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
        summarise(m = median(wwImp)) %>% ungroup() %>%
        mutate(
            days1 = as.numeric(scale(days1)),
            days2 = as.numeric(scale(days2)),
            midTrips = as.numeric(scale(midTrips)),
            ApproxHullSA = as.numeric(scale(ApproxHullSA))
        )
    pred.mat <- matrix(1 - diag(ncol(lvl2)), ncol(lvl2))
    pred.mat[1, ] <- pred.mat[, 1] <- 0
    miLvl2 <- mice(lvl2, m = 1, method = "rf", predictorMatrix = pred.mat,
                   maxit = 20, printFlag = FALSE)
    lvl2 <- complete(miLvl2) %>% select(-m)
    return(lvl2)
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
    ## I've pre-scaled now, so removed what was here (see previous commits).
    return(stanData)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Function to create coefficient summary data
################################################################################
################################################################################
sumMC <- function(draws, qnts){
    if(is.na(ncol(draws))){
        summ <- quantile(draws, qnts, names = FALSE)
        summ <- tibble::as_tibble(t(summ))
        names(summ) <- names(qnts)
    } else {
        summ <- apply(draws, 2, quantile, probs = qnts, names = FALSE)
        summ <- tibble::as_tibble(t(summ))
        names(summ) <- names(qnts)
    }
    summ
}
summRename <- function(summList){
    summ <- lapply(names(summList), function(nm){
        dat <- summList[[nm]]
        nr <- nrow(dat)
        if(nr > 1){
            dat <- dat %>%
                mutate(coef = paste0(nm, seq_len(nr)))
        } else {
            dat <- dat %>%
                mutate(coef = nm)
        }
    })
    bind_rows(summ)
}
################################################################################
################################################################################
