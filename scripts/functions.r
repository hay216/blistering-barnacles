# ==============================================================================
#                           FUNCTIONS
# ==============================================================================
## Function to produce histrogram in a pairs plot
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col='blue', cex.main = 0.8,...)
}

## Function to produce kernel density estimate in a pairs plot
panel.density <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  den <- density(x, na.rm = T)
  rug(x, ticksize = 0.06, lwd = 1, col = 'red')
  y <- den$y; y <- y/max(y)
  x <- den$x
  lines(x, y)
}

## Function to produce correlation coefficients in a pairs plot
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1 <- abs(cor(x, y, use ='pairwise.complete.obs', method = 'pearson'))
  r2 <- abs(cor(x, y, use ='pairwise.complete.obs', method = 'kendall'))
  txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
  txt2 <- format(c(r2, 0.123456789), digits=digits)[1]
  txt1 <- paste('Pearson = ', txt1, sep='')
  txt2 <- paste('Kendall = ',txt2, sep='')
  txt <- paste(prefix, txt1, txt2, sep=" ")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = 1)
  text(0.5, 0.6, txt1)
  text(0.5, 0.4, txt2)
}

## Function to calculate the deviance (2*likelihood ratio) of two mixed models
pboot<-function(m0,m1){
  s <- simulate(m0)
  L0<-logLik(refit(m0,s))
  L1<-logLik(refit(m1,s))
  2*(L1-L0)
}

imputeLRT <- function(h0, h1, imputed.data.object) {
    if (!(length(fixef(h0)) <  length(fixef(h1))))
        stop("The first argument needs the smaller model.")
    imputed.data.list <- imputed.data.object$imputations
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
                                        #  browser()
    return(list(
        D.L = D.L,
        Pval = Pval,
        rL = rL,
        d.L.bar = d.L.bar,
        d.prime.m.bar = d.prime.m.bar,
#    Q.bar.0 = Q.bar.0,
#    Q.bar.1 = Q.bar.1,
        k = k,
        w = w,
        p0 = p0,
        p1 = p1,
        m = m))
}

imputeLRT.c <- cmpfun(imputeLRT)
