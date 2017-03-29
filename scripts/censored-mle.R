################################################################################
################################################################################
## Title: Censored MLE
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Tests fitting a censored model to the biofouling data
## Time-stamp: <2017-03-29 12:33:04 (slane)>
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
packages <- c("dplyr", "rstan", "mice", "ggplot2")
ipak(packages)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()/2)
samplesdata <- read.csv("../data/samples.csv")
source("../scripts/imputation-functions.R")
## Filter out crap, and transform stuff. Scale the logged variables (just cause
## it's easy at the moment, do something better later).
## Double-checked these, and 51 doesn't seem bad? Don't know why it was chosen
## as an outlier?

################################################################################
################################################################################
## Begin Section: Imputation, using MICE + stan in a custom-rolled iterative
## function.
################################################################################
################################################################################
## First time around, use uniform imputation for censored outcome.
## First prepare data.
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
           paintTypeF = factor(paintTypeInt),
           wwImp = wetWeight,
           wwLog = log(ifelse(wetWeight < 1.5, runif(1, 0, 1.5), wetWeight))) %>%
    group_by(boatIDInt) %>% mutate(obsNum = seq_len(n())) %>% ungroup()
## Do the initial imputation
initialImp <- lvl2Imp(data)
## Iterate between stan and mice
imp1 <- iterImp(data, initialImp$impData, "censored-mle-boat2", 2)
## Yay, that works :)
################################################################################
################################################################################

## Have a look at boxplots of imputed data
test <- bind_rows(
    data %>% group_by(boatIDInt) %>%
    summarise(days1S = mean(days1S), days2S = mean(days2S),
              midTripsS = mean(midTripsS)) %>%
    mutate(impNo = 0),
    imp1$lvl2Data %>% select(days1S, days2S, midTripsS) %>%
    mutate(impNo = 1)
)
pl1 <- ggplot(test, aes(x = factor(impNo), y = days1S)) + geom_boxplot()
pl2 <- ggplot(test, aes(x = factor(impNo), y = midTripsS)) + geom_boxplot()
