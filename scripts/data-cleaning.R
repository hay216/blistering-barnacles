################################################################################
################################################################################
## Title: Data cleaning and imputation
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Cleans data for manuscript and model fitting, and performs
## imputation on the vessel level.
## Time-stamp: <2017-04-28 16:21:48 (slane)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("dplyr", "mice", "parallel")
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
## First prepare data.
data <- samplesdata %>% filter(LocID %in% c("HA", "PJ", "HP")) %>%
    select(-boatType, -boatCode, -paintRating, -Location, -LocCode, -samLab,
           -wetWeight1, -wetWeight2)
vessels <- vessels %>%
    distinct(boatID, .keep_all = TRUE) %>%
    select(boatID, samLoc, boatType, ApproxHullSA) %>%
    filter(!is.na(boatID))
## Merge on vessel data
data <- left_join(
    data,
    vessels,
    by = "boatID"
) %>%
    mutate(
        samLoc = recode(samLoc,
                        "SaNAringham Yacht Club" = "Sandringham Yacht Club",
                        "Hobsons Bay Yacht Club " = "Hobsons Bay Yacht Club"),
        boatType = recode(boatType,
                          "Fishing vessel (Lobster/scallop)" = "Fishing",
                          "Fishing vessel (Abalone mothership)" = "Fishing",
                          "Fishing vessel (Long line)" = "Fishing",
                          "Motor cruiser" = "Motor cruiser/Other",
                          "Ferry" = "Motor cruiser/Other",
                          "Tug" = "Motor cruiser/Other"),
        boatType = as.factor(boatType),
        LocID = recode(LocID,
                       "HA" = "Hull", "HP" = "Keel", "PJ" = "Rudder"),
        LocID = as.factor(LocID),
        cens = ifelse(wetWeight < 1.5, 1, 0),
        paintType = as.factor(paintType)
    )
## Data for imputations/modelling
impData <- data %>% select(-samLoc, -cens, -LocID)
## Level 1 data
lvl1Data <- data %>% select(boatID, wetWeight, LocID, cens)
## Loop to create multiple imputations - give a loop as I want to start each
## imputation off with a random draw from a U(0, 1.5) for the censored data just
## to inject a little randomness into it.
## Create 20 imputations, join to full data, and also create stan data
set.seed(787, "L'Ecuyer")
impList <- mclapply(1:20, function(i){
    imp <- lvl2Imp(impData)
    stanData <- createStanData(lvl1Data, imp)
    list(lvl2 = imp, stanData = stanData)
})
## Save as rds for further use.
if(!dir.exists("../data/")) dir.create("../data/")
saveRDS(data, file = "../data/biofouling.rds")
saveRDS(impList, file = "../data/imputations.rds")
################################################################################
################################################################################
