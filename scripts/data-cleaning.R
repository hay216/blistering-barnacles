################################################################################
################################################################################
## Title: Data cleaning
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Cleans data for manuscript and model fitting.
## Time-stamp: <2017-04-12 14:39:15 (slane)>
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
packages <- c("dplyr")
ipak(packages)
samplesdata <- read.csv("../data-raw/samples.csv")
## Bring in vessel data as well.
vessels <- read.csv("../data-raw/vessel.csv")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Filter out locations not used and mutate data as necessary.
################################################################################
################################################################################
## First time around, use uniform imputation for censored outcome.
## First prepare data.
data <- samplesdata %>% filter(LocID %in% c("HA", "PJ", "HP")) %>%
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
## Merge on vessel data
data <- left_join(
    data,
    vessels
) %>%
    mutate(
        samLoc = recode(samLoc,
                        "SaNAringham Yacht Club" = "Sandringham Yacht Club",
                        "Hobsons Bay Yacht Club " = "Hobsons Bay Yacht Club"),
        boatType = recode(boatType,
                          "Fishing vessel (Lobster/scallop)" = "Fishing",
                          "Fishing vessel (Abalone mothership)" = "Fishing",
                          "Fishing vessel (Long line)" = "Fishing",
                          "Ferry" = "Other",
                          "Tug" = "Other")
    )
## Save as rds for further use.
if(!dir.exists("../data/")) dir.create("../data/")
saveRDS(data, file = "../data/biofouling.rds")
################################################################################
################################################################################
