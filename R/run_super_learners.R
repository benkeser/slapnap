#! /usr/bin/env Rscript

# load libraries
library(SuperLearner)

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"

# load data
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
dat <- dat[complete.cases(dat),]

# may need at some point to source in code to do super learner
# library preparation etc...
# source()

# if short run, use a simple library
if(reduce_library){
  SL.library <- c("SL.mean", "SL.glm")
}else{
  SL.library <- c("SL.mean", "SL.glm")  
}

# get names of predictors
non_pred_names <- c("pc.ic50", "pc.ic80", "iip",
                    "dichotomous.1", "dichotomous.2",
                    "seq.id.lanl","seq.id.catnap")
pred_names <- colnames(dat)[!(colnames(dat) %in% non_pred_names)]

set.seed(123125)

#' function to run super learner and cv super learner on a single outcome
#' @param outcome_name String name of outcome
#' @param pred_names Vector of string names of predictor variables
sl_one_outcome <- function(outcome_name, 
                           pred_names,                           
                           save_dir = "/home/slfits/",
                           fit_name = paste0("fit_", outcome_name, ".RData"),
                           cv_fit_name = paste0("cvfit_", outcome_name, ".RData"),
                           reduce_covs = FALSE,
                           ...){
        pred <- dat[ , pred_names]

        if(reduce_covs){
          pred <- pred[ , 1:3]
        }

        fit <- SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
        save(fit, file = paste0(save_dir, fit_name))

        cv_fit <- CV.SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
        save(cv_fit, file = paste0(save_dir, cv_fit_name))
        return(invisible(NULL))
}

sl_ic50 <- sl_one_outcome(outcome_name = "pc.ic50",
                          pred_names = pred_names,
                          family = "gaussian",
                          SL.library = SL.library, 
                          cvControl = list(V = 10),
                          method = "method.CC_LS",
                          reduce_covs = reduce_covs)

if(!reduce_outcomes){
  sl_ic80 <- sl_one_outcome(outcome_name = "pc.ic80",
                            pred_names = pred_names,
                           family = "gaussian",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_LS")

  sl_iip <- sl_one_outcome(outcome_name = "iip",
                            pred_names = pred_names,
                           family = "gaussian",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_LS")

  sl_dichotomous1 <- sl_one_outcome(outcome_name = "dichotomous.1",
                            pred_names = pred_names,
                           family = "binomial",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_nloglik")

  sl_dichotomous2 <- sl_one_outcome(outcome_name = "dichotomous.2",
                            pred_names = pred_names,
                           family = "binomial",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                         method = "method.CC_nloglik")
}