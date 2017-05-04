################################################################################
################################################################################
## Title: Post processing
## Author: Steve Lane
## Date: Thursday, 04 May 2017
## Synopsis: Post process the output from the regression models
## Time-stamp: <2017-05-04 20:49:04 (slane)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "rstan", "loo", "ggplot2")
ipak(packages)
theme_set(theme_bw())
m0 <- readRDS("../data/censored-mle-m0.rds")
m1 <- readRDS("../data/censored-mle-m1.rds")
m2 <- readRDS("../data/censored-mle-m2.rds")
m3 <- readRDS("../data/censored-mle-m3.rds")
m4 <- readRDS("../data/censored-mle-m4.rds")
imps <- readRDS("../data/imputations.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: M1 figures
################################################################################
################################################################################
a1 <- extract(m1, "alphaBoat")$alphaBoat
a1Sum <- t(apply(a1, 2, quantile, probs = c(0.1, 0.5, 0.9)))
a1Dat <- tibble(low = a1Sum[,1], mid = a1Sum[,2], high = a1Sum[,3],
                boatID = 1:nrow(a1Sum))
set.seed(13)
lvl2 <- imps[[sample(seq_along(imps), 1)]]$lvl2
a1Dat <- left_join(a1Dat, lvl2) %>%
    gather(type, value, -low, -mid, -high, -boatID, -paintType, -boatType) %>%
    mutate(type = recode(type,
                         days1 = "Days since last used",
                         days2 = "Days since last cleaned",
                         midTrips = "Median number of trips",
                         ApproxHullSA = "Hull surface area")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m1, c("betaDays1", "betaDays2", "betaMidTrips",
                        "betaHullSA"))
slopes <- sapply(slopes, quantile, probs = 0.5)
slopes <- tibble(intercept = 0, slope = slopes,
                 type = unique(a1Dat$type))
plM1boat <- ggplot(a1Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                              colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept), data = slopes) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1boat.pdf", plM1boat)
plM1paint <- ggplot(a1Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                               colour = `Paint type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept), data = slopes) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1paint.pdf", plM1paint)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: M3 figures
################################################################################
################################################################################
a3 <- extract(m3, "alphaBoat")$alphaBoat
a3Sum <- t(apply(a3, 2, quantile, probs = c(0.1, 0.5, 0.9)))
a3Dat <- tibble(low = a3Sum[,1], mid = a3Sum[,2], high = a3Sum[,3],
                boatID = 1:nrow(a3Sum))
a3Dat <- left_join(a3Dat, lvl2) %>%
    select(-ApproxHullSA, -days2) %>%
    gather(type, value, -low, -mid, -high, -boatID, -paintType, -boatType) %>%
    mutate(type = recode(type,
                         days1 = "Days since last used",
                         midTrips = "Median number of trips")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m3, c("betaDays1", "betaDays2", "betaMidTrips",
                        "betaPaint", "betaType", "betaDaysType",
                        "betaTripsType", "betaTripsPaint"))
intType <- c(0, median(slopes$betaType[,1]), median(slopes$betaType[,2]))
slopeDaysType <- c(median(slopes$betaDays1),
                   median(slopes$betaDays1 + slopes$betaDaysType[,1]),
                   median(slopes$betaDays1 + slopes$betaDaysType[,2]))
slopeTripsType <- c(median(slopes$betaMidTrips),
                    median(slopes$betaMidTrips + slopes$betaTripsType[,1]),
                    median(slopes$betaMidTrips + slopes$betaTripsType[,2]))
parsType <- tibble(intercept = rep(intType, 2),
                   slope = c(slopeDaysType, slopeTripsType),
                   `Vessel type` = factor(rep(levels(a3Dat$`Vessel type`), 2)),
                   type = rep(c("Days since last used",
                                "Median number of trips"), each = 3))
plM3boat <- ggplot(a3Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                              colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type, ncol = 2) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Vessel type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3boat.pdf", plM3boat)
intPaint <- c(0, median(slopes$betaPaint[,1]), median(slopes$betaPaint[,2]))
slopeTripsPaint <- c(median(slopes$betaMidTrips),
                     median(slopes$betaMidTrips + slopes$betaTripsPaint[,1]),
                     median(slopes$betaMidTrips + slopes$betaTripsPaint[,2]))
parsType <- tibble(intercept = intPaint,
                   slope = slopeTripsPaint,
                   `Paint type` = factor(levels(a3Dat$`Paint type`)),
                   type = rep("Median number of trips", each = 3))
plM3paint <- ggplot(a3Dat %>% filter(type == "Median number of trips"),
                    aes(x = value, y = mid, ymin = low, ymax = high,
                        colour = `Paint type`)) +
    geom_pointrange(fatten = 1) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Paint type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3paint.pdf", plM3paint)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: M4 figures
################################################################################
################################################################################
a4 <- extract(m4, "alphaBoat")$alphaBoat
a4Sum <- t(apply(a4, 2, quantile, probs = c(0.1, 0.5, 0.9)))
a4Dat <- tibble(low = a4Sum[,1], mid = a4Sum[,2], high = a4Sum[,3],
                boatID = 1:nrow(a4Sum))
a4Dat <- left_join(a4Dat, lvl2) %>%
    select(-ApproxHullSA, -midTrips) %>%
    gather(type, value, -low, -mid, -high, -boatID, -paintType, -boatType) %>%
    mutate(type = recode(type,
                         days1 = "Days since last used",
                         days2 = "Days since last cleaned",
                         midTrips = "Median number of trips")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m4, c("betaDays1", "betaDays2", "betaType",
                        "betaDaysType"))
intType <- c(0, median(slopes$betaType[,1]), median(slopes$betaType[,2]))
slopeDaysType <- c(median(slopes$betaDays1),
                   median(slopes$betaDays1 + slopes$betaDaysType[,1]),
                   median(slopes$betaDays1 + slopes$betaDaysType[,2]))
parsType <- tibble(intercept = intType, slope = slopeDaysType,
                   `Vessel type` = factor(levels(a4Dat$`Vessel type`)),
                   type = rep("Days since last used", 3))
plM4type <- ggplot(a4Dat %>% filter(type == "Days since last used"),
                   aes(x = value, y = mid, ymin = low, ymax = high,
                       colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Vessel type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM4type.pdf", plM4type)
plM4days2 <- ggplot(a4Dat %>% filter(type == "Days since last cleaned"),
                    aes(x = value, y = mid, ymin = low, ymax = high,
                        colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    geom_abline(aes(slope = median(slopes$betaDays2), intercept = intercept,
                    colour = `Vessel type`), data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM4days2.pdf", plM4days2)
################################################################################
################################################################################
