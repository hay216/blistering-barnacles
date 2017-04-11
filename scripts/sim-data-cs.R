################################################################################
################################################################################
## Title: Simulate Data (Cross-Sectional)
## Author: Steve Lane
## Date: Tuesday, 11 April 2017
## Synopsis: Script to generate data for simulation study. Generates a single
## data set, so that it can be called from a master function. This means it can
## be made reproducible either in serial or parallel (e.g. using
## RNGkind("L'Ecuyer-CMRG") from the parallel package).
## Time-stamp: <2017-04-11 17:15:14 (slane)>
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

## Two independent normals and a categorical (three levels)
n <- 100
## Random normals
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)
## Categorical
x3 <- sample(c("a", "b", "c"), n, replace = TRUE, prob = c(0.2, 0.3, 0.5))
## Categorical effects (a = 1, b = 0, c = -1)
e3 <- ifelse(x3 == "a", 1, ifelse(x3 == "b", 0, -1))
## Substantive model
y <- exp(3 + x1 - x2 + e3 + rnorm(n, 0, 2))
## Test lm
ly <- log(y)
lm1 <- lm(ly ~ x1 + x2 + x3)
