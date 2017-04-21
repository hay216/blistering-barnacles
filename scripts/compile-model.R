#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
################################################################################
################################################################################
## Title: Compile model
## Author: Steve Lane
## Date: Friday, 21 April 2017
## Time-stamp: <2017-04-21 09:42:24 (slane)>
## Synopsis: Script that compiles the stan model. Designed to be called
## from the Makefile, it requires the model name to be set on the command line,
## or prior to sourcing the script.
################################################################################
################################################################################
if(length(args) != 1){
    stop("One argument must be supplied: stan model name.\nRscript compile-model.R mname=model-name", call. = FALSE)
} else {
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
packages <- "rstan"
ipak(packages)
## Compile stan model
rstan_options(auto_write = TRUE)
model <- stan_model(paste0("../stan/", mname, ".stan"))
################################################################################
################################################################################
