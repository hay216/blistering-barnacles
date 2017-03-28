################################################################################
################################################################################
## Title: Censored MLE
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Tests fitting a censored model to the biofouling data
## Time-stamp: <2017-03-29 07:45:46 (slane)>
################################################################################
################################################################################
ipak <- function(pkg){
    ## Check for github packages (throw away github username)
    chk.git <- gsub(".*/", "", pkg)    
    new.pkg <- pkg[!(chk.git %in% installed.packages()[, "Package"])]
    if(!(length(new.pkg) == 0)){
        git.ind <- grep("/", new.pkg)
        if(length(git.ind) == 0){
            install.packages(new.pkg, dependencies = TRUE)
        } else {
            devtools::install_github(new.pkg[git.ind])
        }
    }
    sapply(chk.git, require, character.only = TRUE)
}
## Add github packages using gitname/reponame format
packages <- c("dplyr", "rstan", "mice", "randomForest", "forcats")
ipak(packages)
samplesdata <- read.csv("../data/samples.csv")
## Filter out crap, and transform stuff. Scale the logged variables (just cause
## it's easy at the moment, do something better later).
## Double-checked these, and 51 doesn't seem bad? Don't know why it was chosen
## as an outlier?
data1 <- samplesdata %>% filter(!(boatID %in% c(24, 51)),
                                LocID %in% c("HA", "PJ", "HP")) %>%
    select(-boatType, -boatCode, -paintRating, -Location, -LocCode, -samLab,
           -wetWeight1, -wetWeight2) %>% na.omit() %>%
    mutate(cens = ifelse(wetWeight < 1.5, 1, 0),
           paintTypeInt = as.integer(factor(paintType)),
           locIDInt = as.integer(factor(LocID)),
           days1S = as.numeric(scale(log(days1 + 0.1))),
           days2S = as.numeric(scale(log(days2 + 0.1))),
           midTripsS = as.numeric(scale(log(midTrips + 0.1))),
           boatIDInt = as.numeric(as.factor(boatID)))
obsData <- data1 %>% filter(cens == 0)
censData <- data1 %>% filter(cens == 1)
stanData <- c(with(obsData, list(N = nrow(obsData), days1 = days1S,
                                 days2 = days2S, midTrips = midTripsS,
                                 numPaint = max(paintTypeInt),
                                 paintType = paintTypeInt,
                                 numLoc = max(locIDInt),
                                 locID = locIDInt, Y = wetWeight)),
              with(censData, list(nCens = nrow(censData), days1Cens = days1S,
                                   days2Cens = days2S, midTripsCens = midTripsS,
                                   paintTypeCens = paintTypeInt,
                                   locIDCens = locIDInt,
                                   U = 1.5)))

################################################################################
################################################################################
## Begin Section: Compile stan model
################################################################################
################################################################################
## At present, I'm ignoring the clustering...
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()/2)
model <- stan_model("../stan/censored-mle.stan")
## Try fitting the MLE
mleOutput <- optimizing(model, data = stanData)
## HMC sampling
hmcOutput <- sampling(model, data = stanData,
                      pars = c("mu", "beta1", "beta2", "beta3", "alpha1",
                               "alpha2", "sigma"), iter = 500, chains = 4)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Now include clustering at the boat level.
################################################################################
################################################################################
stanData2 <- c(stanData, list(boatID = obsData$boatIDInt),
               list(boatIDCens = censData$boatIDInt),
               list(numBoat = max(obsData$boatIDInt)))
model2 <- stan_model("../stan/censored-mle-boat.stan")
## Try fitting the MLE
mleOutput2 <- optimizing(model2, data = stanData2)
## HMC sampling
system.time(
hmcOutput2 <- sampling(model2, data = stanData2,
                       pars = c("mu", "beta1", "beta2", "beta3", "alpha1",
                                "alpha2", "alphaBoat", "sigma_alpha1",
                                "sigma_alpha2", "sigma_alphaBoat",
                                "sigma", "muHatCens"), iter = 500, chains = 4)
)
## That works, but unfortunately predictions are too large (this is due to the
## generated data section not obeying the censoring).
## Function for predictions
predFun <- function(model, U){
    muHat <- extract(model, "muHatCens")$muHatCens
    sigma <- extract(model, "sigma")$sigma
    predY <- sapply(seq_len(nrow(muHat)), function(i){
        upper <- plnorm(U, muHat[i,], sigma[i])
        u <- sapply(upper, function(u) runif(1, 0, u))
        qlnorm(u, muHat[i,], sigma[i])
    })
    t(predY)
}
yCens2 <- predFun(hmcOutput2, 1.5)

model3 <- stan_model("../stan/censored-mle-boat2.stan")
## HMC sampling
system.time(
hmcOutput3 <- sampling(model3, data = stanData2,
                       pars = c("mu", "beta1", "beta2", "beta3", "alpha1",
                                "alpha2", "alphaBoat", "sigma_alpha1",
                                "sigma_alpha2", "sigma_alphaBoat",
                                "sigma", "yCens"), iter = 1000, chains = 4,
                       control = list(adapt_delta = 0.9))
)
## This works, but adaptation is slower and I get more divergent transitions.
yCens3 <- extract(hmcOutput3, "yCens")$yCens

################################################################################
################################################################################
## Begin Section: Manual imputation, using MICE
################################################################################
################################################################################
## Take output from hmcOutput3, average within boat, then use boat-level data
## for imputation.
## First prepare data
data <- samplesdata %>% filter(!(boatID %in% c(24, 51)),
                               LocID %in% c("HA", "PJ", "HP")) %>%
    select(-boatType, -boatCode, -paintRating, -Location, -LocCode, -samLab,
           -wetWeight1, -wetWeight2) %>%
    mutate(cens = ifelse(wetWeight < 1.5, 1, 0),
           paintTypeInt = as.integer(factor(paintType)),
           locIDInt = as.integer(factor(LocID)),
           days1S = as.numeric(scale(log(days1 + 0.1))),
           days2S = as.numeric(scale(log(days2 + 0.1))),
           midTripsS = as.numeric(scale(log(midTrips + 0.1))),
           boatIDInt = as.numeric(as.factor(boatID)),
           paintTypeF = factor(paintTypeInt))
lvl2 <- data %>% group_by(days1S, days2S, midTripsS, paintTypeF) %>%
    summarise(m = mean(wetWeight)) %>% ungroup()
miLvl2 <- mice(lvl2, m = 1, method = c(rep("pmm", 3), "polyreg", "pmm"),
               maxit = 10)
miLvl2 <- mice(lvl2, m = 1, method = "rf", maxit = 10)
## Not gettting anything for paint type... that's because it needs to be a
## factor :) fixed now.
## days1, midTrips and paintType are the variables with missing data. What about
## using randomForests?
## First, put in some random values.
lvl2 <- lvl2 %>%
    mutate(days1Miss = ifelse(is.na(days1S), 1, 0),
           midTripsMiss = ifelse(is.na(midTripsS), 1, 0),
           paintTypeMiss = ifelse(is.na(paintType), 1, 0),
           days1S = ifelse(is.na(days1S), mean(days1S, na.rm = TRUE), days1S),
           midTripsS = ifelse(is.na(midTripsS), mean(midTripsS, na.rm = TRUE),
                              midTripsS),
           paintType = factor(ifelse(is.na(paintType), "Hard", paintType)))
## Now fit models
## First days1
miDays1 <- randomForest(days1S ~ days2S + midTripsS + paintType + m, data = lvl2,
                        ntree = 2000)
## Predict where days1 is actually missing
lvl2 <- bind_rows(
    lvl2 %>% filter(days1Miss == 0),
    lvl2 %>% filter(days1Miss == 1) %>%
    mutate(days1 = predict(miDays1, newdata = .))
)
## midTrips
miMidTrips <- randomForest(midTripsS ~ days1S + days2S + paintType + m,
                           data = lvl2, ntree = 2000)
## Predict where days1 is actually missing
lvl2 <- bind_rows(
    lvl2 %>% filter(midTripsMiss == 0),
    lvl2 %>% filter(midTripsMiss == 1) %>%
    mutate(midTrips = predict(miMidTrips, newdata = .))
)
## paintType
miPaintType <- randomForest(paintType ~ days1S + days2S + midTripsS + m,
                            data = lvl2, ntree = 2000)
## Predict where days1 is actually missing
lvl2 <- bind_rows(
    lvl2 %>% filter(paintTypeMiss == 0),
    lvl2 %>% filter(paintTypeMiss == 1) %>%
    mutate(paintType = predict(miPaintType, newdata = .))
)

## OK, now join this back on to the full data and fit in Stan
impData <- left_join(
    data %>% select(boatIDInt, wetWeight, cens, locIDInt),
    lvl2 %>% mutate(boatIDInt = seq_len(nrow(lvl2)),
                    paintTypeInt = as.integer(factor(paintType)))
)
obsDataImp <- impData %>% filter(cens == 0)
censDataImp <- impData %>% filter(cens == 1)
stanDataImp <- c(with(obsDataImp, list(N = nrow(obsDataImp), days1 = days1S,
                                       days2 = days2S, midTrips = midTripsS,
                                       numPaint = max(paintTypeInt),
                                       paintType = paintTypeInt,
                                       numLoc = max(locIDInt),
                                       locID = locIDInt, Y = wetWeight,
                                       boatID = boatIDInt,
                                       numBoat = max(boatIDInt))),
                 with(censDataImp, list(nCens = nrow(censDataImp),
                                        days1Cens = days1S,
                                        days2Cens = days2S,
                                        midTripsCens = midTripsS,
                                        paintTypeCens = paintTypeInt,
                                        locIDCens = locIDInt,
                                        boatIDCens = boatIDInt,
                                        U = 1.5)))

system.time(
hmcImp1 <- sampling(model3, data = stanDataImp,
                    pars = c("mu", "beta1", "beta2", "beta3", "alpha1",
                             "alpha2", "alphaBoat", "sigma_alpha1",
                             "sigma_alpha2", "sigma_alphaBoat",
                             "sigma", "yCens"), iter = 1000, chains = 4,
                    control = list(adapt_delta = 0.9))
)
## This works, so I now need to automate it...
