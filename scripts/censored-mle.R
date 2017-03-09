################################################################################
################################################################################
## Title: Censored MLE
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Tests fitting a censored model to the biofouling data
## Time-stamp: <2017-03-09 16:43:27 (slane)>
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
packages <- c("dplyr", "rstan")
ipak(packages)
samplesdata <- read.csv("../data/samples.csv")
## Filter out crap, and transform stuff. Scale the logged variables (just cause
## it's easy at the moment, do something better later).
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
                                   U = rep(1.5, nrow(censData)))))

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
hmcOutput2 <- sampling(model2, data = stanData2,
                       pars = c("mu", "beta1", "beta2", "beta3", "alpha1",
                                "alpha2", "alphaBoat", "sigma_alpha1",
                                "sigma_alpha2", "sigma_alphaBoat",
                                "sigma", "yCens"), iter = 500, chains = 4)
## That works, but unfortunately
