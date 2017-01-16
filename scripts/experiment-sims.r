################################################################################
# 
# Simulation script
# 
# Andrew Robinson 15-April-2015
#
################################################################################


the.seed <- 10000000

reps <- 1000
cores <- 25

MAIL <- TRUE

mail <- function(address, subject, message, attach = NULL) {
  if (MAIL)
    if (is.null(attach)) {
      system(paste("echo '", message,
                   "' | mutt -s '", subject,
                   "' ", address, sep=""))
    } else {
      system(paste("echo '", message,
                   "' | mutt -s '", subject, "'",
                   " ", address,
                   paste(" -a ", attach, collapse=" "), sep=""))
    }
}

mail.and.stop <- function() {
  mail("mensurationist@gmail.com",
       "Sims on server: ERROR.",
       "Stopped.")
  stop()
}

options(error = mail.and.stop)


library(lme4)    # For mixed-effects models

library(Amelia)  # For imputation

library(equivalence) # For data that I know

library(parallel) # Speed things up

data(ufc)

str(ufc)

ufc <- subset(ufc, !is.na(Dbh))
ufc$Species[ufc$Species %in% c("F","FG")] <- "GF"
ufc$Species <- factor(ufc$Species)
ufc$ht.measured <- is.na(ufc$Height)

ufc <- subset(ufc,
              Species %in% names(sort(-table(Species)))[1:6],
              select = c("ht.measured","Plot","Dbh.in","Species"))

## Fit a pair of models to compare using LRT

ufc.whole <- ufc

isna.1.w <- glmer(ht.measured ~ Dbh.in + Species + (1 | Plot),
                  data = ufc.whole, family = binomial)
isna.0.w <- glmer(ht.measured ~ Dbh.in + (1 | Plot),
                  data = ufc.whole, family = binomial)

anova(isna.0.w, isna.1.w)

################################################################################

imputeLRT <- function(h0, h1, imputed.data.list) {
  m <- length(imputed.data.list)
  h0.models <- lapply(imputed.data.list,
                      function (dataset) {
                        update(h0, data = dataset)
                      })
  h0.dev.funs <- lapply(imputed.data.list,
                        function (dataset) {
                          update(h0, data = dataset,
                                 devFunOnly = TRUE)
                        })
  Q.bar.0 <-
    colMeans(do.call(rbind,
                     lapply(h0.models,
                            function (fit) {
                              c(getME(fit, "theta"), fixef(fit))
                            })))  
  h1.models <- lapply(imputed.data.list,
                      function (dataset) {
                        update(h1, data = dataset)
                      })
  h1.dev.funs <- lapply(imputed.data.list,
                        function (dataset) {
                          update(h1, data = dataset,
                                 devFunOnly = TRUE)
                        })
  Q.bar.1 <-
    colMeans(do.call(rbind,
                     lapply(h1.models,
                            function (fit) {
                              c(getME(fit, "theta"), fixef(fit))
                            })))
  d.prime.m.bar <- mean(unlist(lapply(1:m,
                                      function(i) {
                                        anova(h0.models[[i]],
                                              h1.models[[i]])$Chisq[2]
                                      })))
  d.L.bar <- mean(unlist(lapply(1:m,
                                function(i) {
                                  h0.dev.funs[[i]](Q.bar.0) - 
                                    h1.dev.funs[[i]](Q.bar.1) 
                                })))
  p0 <- length(Q.bar.0)
  p1 <- length(Q.bar.1)
  k <- p1 - p0
  rL <- (m + 1) / ((m - 1) * k) * (d.prime.m.bar - d.L.bar) # 3.8
  D.L <- d.L.bar / ((1 + rL) * k)                           # 3.7
  v <- k * (m - 1) 
  w <- ifelse(v > 4,                                        # 2.7 
              4 + (v - 4) * (1 + (1 - v/2) / rL)^2,
              v / 2 * (1 + 1/k) * (1 + 1/rL)^2)
  Pval <- 1 - pf(D.L, k, w)
  return(list(
    D.L = D.L,
    Pval = Pval,
    rL = rL,
    d.L.bar = d.L.bar,
    d.prime.m.bar = d.prime.m.bar,
    Q.bar.0 = Q.bar.0,
    Q.bar.1 = Q.bar.1,
    k = k,
    w = w,
    p0 = p0,
    p1 = p1,
    m = m))
}


isna.1 <- glmer(ht.measured ~ Dbh.in + Species + (1 | Plot),
                data = ufc, family = binomial)
isna.0 <- glmer(ht.measured ~ Dbh.in + (1 | Plot),
                data = ufc, family = binomial)  

imputeLRT(isna.0, isna.1,
          list(ufc, ufc, ufc, ufc, ufc))

anova(isna.0, isna.1)



simulate <- function(m = 5, n = 30) {
  make.missing <- sample(1:nrow(ufc), size = n * 3, replace = FALSE)
  is.na(ufc$Species[make.missing[1:n]]) <- TRUE
  is.na(ufc$Dbh.in[make.missing[(n+1):(2*n)]]) <- TRUE
  is.na(ufc$ht.measured[make.missing[(2*n+1):(3*n)]]) <- TRUE
  isna.1 <- glmer(ht.measured ~ Dbh.in + Species + (1 | Plot),
                  data = ufc, family = binomial)
  isna.0 <- glmer(ht.measured ~ Dbh.in + (1 | Plot),
                  data = ufc, family = binomial)  
  naive.p <- anova(isna.0, isna.1)[["Pr(>Chisq)"]][2]
  ## Impute datasets using Amelia
  ufc.imputes <- amelia(ufc,
                        m = m,
                        cs = "Plot",
                        noms = "Species",
                        p2s = 0)
  imputed <- imputeLRT(isna.0, isna.1, ufc.imputes$imputations)
  return(list(naive.p = naive.p, imputed.p = imputed$Pval))
}

simulate()

frame <- expand.grid(ms = c(2, 5, 10, 20, 50, 100),
                     ns = c(5, 10, 20, 30, 50, 100),
                     reps = 1:reps)

results <- mclapply(1:nrow(frame),
                    function(i) 
                    simulate(frame$ms[i], frame$ns[i]),
                    mc.cores = cores)

results <- do.call(rbind, results)

experiment <- cbind(frame, results)

save(experiment, file = "mi.RData")

mail(address = "andrewpr",
     subject = "MI simulations",
     message = "All done.",
     attach = "mi.RData")


