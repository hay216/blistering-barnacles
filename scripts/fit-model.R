#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
################################################################################
################################################################################
## Title: Fit model
## Author: Steve Lane
## Date: Friday, 21 April 2017
## Time-stamp: <2017-04-21 11:32:43 (slane)>
## Synopsis: Script that drives the censored regression model. Designed to be
## called from the Makefile, it requires the model name, a seed for rng, and
## number of iterations to be set on the command line, or prior to sourcing the
## script.
################################################################################
################################################################################
if(!(length(args) %in% 2:4)){
    stop("Four arguments must be supplied: model name, myseed, scaled and iter.\n scaled (use scaled inputs) is optional and defaults to TRUE.\n iter (for HMC) is optional, and defaults to 4000.\nRscript fit-model.R mname=model-name myseed=my-seed scaled=1 iter=num-iter",
         call. = FALSE)
} else {
    if(length(args) == 2){
        ## Default if option not specified
        iter <- 4000
        scaled <- 1
    }
    hasOpt <- grepl("=", args)
    argLocal <- strsplit(args[hasOpt], "=")
    for(i in seq_along(argLocal)){
        value <- NA
        tryCatch(value <- as.double(argLocal[[i]][2]), warning = function(e){})
        if(!is.na(value)){
            ## Assume int/double
            assign(argLocal[[i]][1], value, inherits = TRUE)
        } else {
            assign(argLocal[[i]][1], argLocal[[i]][2], inherits = TRUE)
        }
    }
}
ipak <- function(pkg){
    ## Check for github packages (throw away github username)
    chk.git <- gsub(".*/", "", pkg)    
    new.pkg <- pkg[!(chk.git %in% installed.packages()[, "Package"])]
    if(!(length(new.pkg) == 0)){
        git.ind <- grep("/", new.pkg)
        if(length(git.ind) == 0){
            install.packages(new.pkg, dependencies = TRUE,
                             repos = "https://cran.csiro.au")
        } else {
            devtools::install_github(new.pkg[git.ind])
        }
    }
    sapply(chk.git, require, character.only = TRUE)
}
## Add github packages using gitname/reponame format
packages <- c("rstan", "parallel")
ipak(packages)
## Load stan model (it should already be compiled from the compile-model.R
## script)
rstan_options(auto_write = TRUE)
## Want cores to be one, we're only running one chain, then combining. Each
## imputation will be sent out via mclapply.
options(mc.cores = parallel::detectCores()/2)
model <- stan_model(paste0("../stan/", mname, ".stan"))
## Load data
impList <- readRDS("../data/imputations.rds")
set.seed(myseed)
if(scaled == 1){
    out <- mclapply(impList, function(dat){
        locMod <- sampling(model, data = dat$stanDataSc, iter = iter,
                           chains = 1, cores = 1)
        locMod
    })
    outname <- paste0("../data/", mname, "-scaled.rds")
} else {
    out <- mclapply(impList, function(dat){
        locMod <- sampling(model, data = dat$stanData, iter = iter,
                           chains = 1, cores = 1)
        locMod
    })
    outname <- paste0("../data/", mname, ".rds")
}
output <- sflist2stanfit(out)
saveRDS(output, file = outname)
################################################################################
################################################################################
