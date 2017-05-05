################################################################################
################################################################################
## Title: Post processing
## Author: Steve Lane
## Date: Thursday, 04 May 2017
## Synopsis: Post process the output from the regression models
## Time-stamp: <2017-05-05 14:23:13 (slane)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "rstan", "loo", "ggplot2", "RColorBrewer")
ipak(packages)
## This is set for readable text when included at half page width.
theme_set(theme_bw())
m0 <- readRDS("../data/censored-mle-m0.rds")
m1 <- readRDS("../data/censored-mle-m1.rds")
m2 <- readRDS("../data/censored-mle-m2.rds")
m3 <- readRDS("../data/censored-mle-m3.rds")
m4 <- readRDS("../data/censored-mle-m4.rds")
imps <- readRDS("../data/imputations.rds")
biofoul <- readRDS("../data/biofouling.rds")
vessels <- biofoul %>% distinct(boatID, .keep_all = TRUE)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Observed weight figures
################################################################################
################################################################################
histData <- biofoul %>%
    filter(wetWeight >= 1.5) %>%
    mutate(wwLog = log(wetWeight)) %>%
    select(LocID, `Weight (gm)` = wetWeight,
           `Weight (gm), log-scale` = wwLog) %>%
    gather("logged", "ww", 2:3)
plHist <- ggplot(histData, aes(x = ww)) + geom_histogram(bins = 11) +
    facet_grid(LocID ~ logged, scales = "free_x") +
    xlab("") +
    ylab("Count") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/obs-hist.pdf", plHist)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Imputation figures
################################################################################
################################################################################
lvl2Imp <- lapply(imps, function(x) x$lvl2) %>% bind_rows() %>%
    mutate(nummi = rep(seq_along(imps), each = nrow(vessels)))
vessImps <- vessels %>%
    mutate(
        days1 = as.numeric(scale(days1)),
        days2 = as.numeric(scale(days2)),
        midTrips = as.numeric(scale(midTrips)),
        ApproxHullSA = as.numeric(scale(ApproxHullSA)),
        nummi = 0
    )
vessImps <- bind_rows(vessImps, lvl2Imp)
plDays1 <- ggplot(vessImps, aes(x = factor(nummi), y = days1)) +
    geom_boxplot(outlier.size = 0.5) +
    xlab("Imputation number") +
    ylab("") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/imp-days1.pdf", plDays1, width = 4.9, height = 4.9)
plMidTrips <- ggplot(vessImps, aes(x = factor(nummi), y = midTrips)) +
    geom_boxplot(outlier.size = 0.5) +
    xlab("Imputation number") +
    ylab("") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/imp-trips.pdf", plMidTrips, width = 4.9, height = 4.9)
set.seed(76)
impSelect <- sample(seq_along(imps), 15)
plPaint <- ggplot(vessImps %>% filter(nummi %in% c(0, impSelect),
                                      !is.na(paintType)),
                  aes(x = paintType)) +
    geom_bar() +
    facet_wrap(~ factor(nummi)) +
    xlab("Anti-fouling paint type") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../graphics/imp-paint.pdf", plPaint, width = 4.9, height = 4.9)
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
    scale_colour_brewer(palette = "Dark2") +
    theme_bw(base_size = 7.7) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1boat.pdf", plM1boat, width = 4.9, height = 4.9)
plM1paint <- ggplot(a1Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                               colour = `Paint type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept), data = slopes) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw(base_size = 7.7) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1paint.pdf", plM1paint, width = 4.9, height = 4.9)
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
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3boat.pdf", plM3boat, height = 3.5)
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
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Paint type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3paint.pdf", plM3paint, width = 3.5, height = 3.5)
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
                         days2 = "Days since last cleaned")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m4, c("betaDays1", "betaDays2", "betaType",
                        "betaDaysType"))
intType <- c(0, median(slopes$betaType[,1]), median(slopes$betaType[,2]))
slopeDaysType <- c(median(slopes$betaDays1),
                   median(slopes$betaDays1 + slopes$betaDaysType[,1]),
                   median(slopes$betaDays1 + slopes$betaDaysType[,2]))
parsType <- tibble(intercept = rep(intType, 2),
                   slope = c(slopeDaysType, rep(median(slopes$betaDays2), 3)),
                   `Vessel type` = factor(rep(levels(a4Dat$`Vessel type`), 2)),
                   type = c(rep("Days since last used", 3),
                            rep("Days since last cleaned", 3)))
parsPaint <- tibble(intercept = intType,
                    slope = rep(median(slopes$betaDays2), 3),
                    `Vessel type` = factor(levels(a4Dat$`Vessel type`)),
                    type = rep("Days since last cleaned"))
plM4type <- ggplot(a4Dat,
                   aes(x = value, y = mid, ymin = low, ymax = high,
                       colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Vessel type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM4type.pdf", plM4type, height = 3.5)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Regression coefficients
################################################################################
################################################################################

## ## Something like the following:
## plotData <- coefs %>%
##         filter(lvl == lookupData$lvl[i]) %>%
##         mutate(label = factor(label, levels = ords)) %>%
##         group_by(model) %>%
##         arrange(label) %>%
##         filter(est < 10, est > -10, lower > -10)
##     plComp <- ggplot(
##         plotData, aes(x = label, y = est, ymin = lower, ymax = upper, colour =
##                           model)
##     ) +
##         geom_pointrange(position = position_dodge(width = 0.75)) +
##         geom_hline(yintercept = 0, linetype = 2) +
##         coord_flip() +
##         ylab("log(Coefficient)") +
##         xlab("Variable")
