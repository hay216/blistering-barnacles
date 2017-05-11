################################################################################
################################################################################
## Title: Post-processing comparison
## Author: Steve Lane
## Date: Thursday, 11 May 2017
## Synopsis: Performs comparison between outcome models.
## Time-stamp: <2017-05-11 17:32:56 (slane)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "tibble", "rstan", "loo", "ggplot2",
              "RColorBrewer")
ipak(packages)
## This is set for readable text when included at half page width.
theme_set(theme_bw())
m0 <- readRDS("../data/censored-mle-m0-robust.rds")
m1 <- readRDS("../data/censored-mle-m1-robust.rds")
m2 <- readRDS("../data/censored-mle-m2-robust.rds")
m3 <- readRDS("../data/censored-mle-m3-robust.rds")
m4 <- readRDS("../data/censored-mle-m4-robust.rds")
m0N <- readRDS("../data/censored-mle-m0.rds")
m1N <- readRDS("../data/censored-mle-m1.rds")
m2N <- readRDS("../data/censored-mle-m2.rds")
m3N <- readRDS("../data/censored-mle-m3.rds")
m4N <- readRDS("../data/censored-mle-m4.rds")
biofoul <- readRDS("../data/biofoul.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: looic table for robust model
################################################################################
################################################################################
m0ll <- extract_log_lik(m0)
m0loo <- loo(m0ll)
m1ll <- extract_log_lik(m1)
m1loo <- loo(m1ll)
m2ll <- extract_log_lik(m2)
m2loo <- loo(m2ll)
m3ll <- extract_log_lik(m3)
m3loo <- loo(m3ll)
m4ll <- extract_log_lik(m4)
m4loo <- loo(m4ll)
looTab <- as.data.frame(compare(m0loo, m1loo, m2loo, m3loo, m4loo))
names(looTab) <- c("LOOIC", "se(LOOIC)", "ELPD", "se(ELPD)", "Eff. P",
                   "se(Eff. P)")
looTab <- looTab %>%
    mutate(Model = c("M4", "M3", "M2", "M1", "M0"))
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: looic for normal model
################################################################################
################################################################################
m0Nll <- extract_log_lik(m0N)
m0Nloo <- loo(m0Nll)
m1Nll <- extract_log_lik(m1N)
m1Nloo <- loo(m1Nll)
m2Nll <- extract_log_lik(m2N)
m2Nloo <- loo(m2Nll)
m3Nll <- extract_log_lik(m3N)
m3Nloo <- loo(m3Nll)
m4Nll <- extract_log_lik(m4N)
m4Nloo <- loo(m4Nll)
looTabN <- as.data.frame(compare(m0Nloo, m1Nloo, m2Nloo, m3Nloo, m4Nloo))
names(looTabN) <- c("LOOIC", "se(LOOIC)", "ELPD", "se(ELPD)", "Eff. P",
                    "se(Eff. P)")
looTabN <- looTabN %>%
    mutate(Model = c("M4", "M3", "M2", "M1", "M0"))
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Join up for output.
################################################################################
################################################################################
compTab <- left_join(looTab, looTabN, by = "Model")
diffs <- rbind(
    compare(m4Nloo, m4loo),
    compare(m3Nloo, m3loo),
    compare(m2Nloo, m2loo),
    compare(m1Nloo, m1loo),
    compare(m0Nloo, m0loo)
    )
compTab <- cbind(compTab, diffs)
compTab$elpd_diff <- 2*compTab$elpd_diff
compTab$se <- sqrt(2)*compTab$se
saveRDS(compTab, "../data/looic-compare.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: PPCs for M0: prop < cens; median; iqr
################################################################################
################################################################################
yPPC <- extract(m0, "y_ppc")$y_ppc
ppc <- apply(yPPC, 1, function(x){
    pr <- mean(x < log(1.5))
    med <- median(exp(x))
    iqr <- IQR(exp(x))
    tibble(`Prop(hat(Y)<1.5)` = pr, `Median(hat(Y))` = med,
           `IQR(hat(Y))` = iqr)
}) %>% bind_rows() %>%
    mutate(model = "O2")
yPPCN <- extract(m0N, "y_ppc")$y_ppc
ppcN <- apply(yPPCN, 1, function(x){
    pr <- mean(x < 1.5)
    med <- median(x)
    iqr <- IQR(x)
    tibble(`Prop(hat(Y)<1.5)` = pr, `Median(hat(Y))` = med,
           `IQR(hat(Y))` = iqr)
}) %>% bind_rows() %>%
    mutate(model = "O1")
allPPC <- bind_rows(ppc, ppcN) %>%
    gather(ppc, value, -model)
## Observed data
obs <- tibble(ppc = c("Prop(hat(Y)<1.5)", "Median(hat(Y))", "IQR(hat(Y))"),
              value = with(biofoul,
                           c(mean(wetWeight < 1.5), median(wetWeight),
                             IQR(wetWeight))))
plPPC <- ggplot(allPPC, aes(x = value)) +
    geom_histogram() +
    facet_grid(model ~ ppc, scales = "free", labeller = label_parsed) +
    geom_vline(aes(xintercept = value), data = obs) +
    xlab("") +
    ylab("Count") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/ppc-compare.pdf", plPPC)
################################################################################
################################################################################
