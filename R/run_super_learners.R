#! /usr/bin/env Rscript

## run the super learners for all outcomes, pre-defined variable groups

# load libraries
library(SuperLearner)
source("variable_groups.R")

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"

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
#' @param save_dir name of directory to save results to
#' @param fit_name name of fits (defaults to fit_<outcome_name>.rds)
#' @param cv_fit_name name of CV fits (defaults to cvfit_<outcome_name>.rds)
#' @param reduce_covs Flag to reduce the number of covariates under consideration
#' @param save_full_object Flag for whether or not to save the full fitted object, or just the fitted values
sl_one_outcome <- function(outcome_name, 
                           pred_names,                           
                           save_dir = "/home/slfits/",
                           fit_name = paste0("fit_", outcome_name, ".rds"),
                           cv_fit_name = paste0("cvfit_", outcome_name, ".rds"),
                           reduce_covs = FALSE,
                           save_full_object = TRUE,
                           ...){
        pred <- dat[ , pred_names]

        if(reduce_covs){
          pred <- pred[ , 1:3]
        }

        fit <- SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
        if (save_full_object) {
            saveRDS(fit, file = paste0(save_dir, fit_name))    
        } 
        saveRDS(fit$SL.predict, file = paste0(save_dir, gsub(".RData", ".rds", gsub("fit_", "fitted_", fit_name))))

        cv_fit <- CV.SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
        if (save_full_object) {
            saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
        } 
        saveRDS(cv_fit$SL.predict, file = paste0(save_dir, gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name))))
        saveRDS(cv_fit$folds, file = paste0(save_dir, gsub("cvfitted_", "cvfolds_", gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name)))))
        return(invisible(NULL))
}

## run full super learner
sl_ic50 <- sl_one_outcome(outcome_name = "pc.ic50",
                          pred_names = pred_names,
                          family = "gaussian",
                          SL.library = SL.library, 
                          cvControl = list(V = 10),
                          method = "method.CC_LS",
                          reduce_covs = reduce_covs)
## run super learners on pre-defined groups
all_var_groups <- get_variable_groups(dat, pred_names)
if (reduce_groups) {
    this_name <- names(all_var_groups)[1]
    sl_ic50_i <- sl_one_outcome(outcome_name = "pc.ic50",
                                    pred_names = pred_names[!(pred_names %in% all_var_groups[[1]])],
                                    fit_name = paste0("fitted_pc.ic50_minus_", this_name, ".rds"),
                                    cv_fit_name = paste0("cvfitted_pc.ic50_minus_", this_name, ".rds"),
                                    family = "gaussian",
                                    SL.library = SL.library,
                                    cvControl = list(V = 10),
                                    method = "method.CC_LS",
                                    reduce_covs = reduce_covs,
                                    save_full_object = FALSE)
} else {
    for (i in 1:length(all_var_groups)) {
        if (length(all_var_groups[i]) == 0) {

        } else {
            this_name <- names(all_var_groups)[i]
            sl_ic50_i <- sl_one_outcome(outcome_name = "pc.ic50",
                                        pred_names = pred_names[!(pred_names %in% all_var_groups[[i]])],
                                        fit_name = paste0("fitted_pc.ic50_minus_", this_name, ".rds"),
                                        cv_fit_name = paste0("cvfitted_pc.ic50_minus_", this_name, ".rds"),
                                        family = "gaussian",
                                        SL.library = SL.library,
                                        cvControl = list(V = 10),
                                        method = "method.CC_LS",
                                        reduce_covs = reduce_covs,
                                        save_full_object = FALSE)
        }
    }
}

if(!reduce_outcomes){
  sl_ic80 <- sl_one_outcome(outcome_name = "pc.ic80",
                            pred_names = pred_names,
                           family = "gaussian",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_LS")
  for (i in 1:length(all_var_groups)) {
        if (length(all_var_groups[i]) == 0) {

        } else {
            this_name <- names(all_var_groups)[i]
            sl_ic80_i <- sl_one_outcome(outcome_name = "pc.ic80",
                                        pred_names = pred_names[!(pred_names %in% all_var_groups[[i]])],
                                        fit_name = paste0("fitted_pc.ic80_minus_", this_name, ".rds"),
                                        cv_fit_name = paste0("cvfitted_pc.ic80_minus_", this_name, ".rds"),
                                        family = "gaussian",
                                        SL.library = SL.library,
                                        cvControl = list(V = 10),
                                        method = "method.CC_LS",
                                        reduce_covs = reduce_covs,
                                        save_full_object = FALSE)
        }
    }
  sl_iip <- sl_one_outcome(outcome_name = "iip",
                            pred_names = pred_names,
                           family = "gaussian",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_LS")
  for (i in 1:length(all_var_groups)) {
        if (length(all_var_groups[i]) == 0) {

        } else {
            this_name <- names(all_var_groups)[i]
            sl_iip_i <- sl_one_outcome(outcome_name = "iip",
                                        pred_names = pred_names[!(pred_names %in% all_var_groups[[i]])],
                                        fit_name = paste0("fitted_iip_minus_", this_name, ".rds"),
                                        cv_fit_name = paste0("cvfitted_iip_minus_", this_name, ".rds"),
                                        family = "gaussian",
                                        SL.library = SL.library,
                                        cvControl = list(V = 10),
                                        method = "method.CC_LS",
                                        reduce_covs = reduce_covs,
                                        save_full_object = FALSE)
        }
    }
  sl_dichotomous1 <- sl_one_outcome(outcome_name = "dichotomous.1",
                            pred_names = pred_names,
                           family = "binomial",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                           method = "method.CC_nloglik")
  for (i in 1:length(all_var_groups)) {
        if (length(all_var_groups[i]) == 0) {

        } else {
            this_name <- names(all_var_groups)[i]
            sl_dichotomous1_i <- sl_one_outcome(outcome_name = "dichotomous.1",
                                        pred_names = pred_names[!(pred_names %in% all_var_groups[[i]])],
                                        fit_name = paste0("fitted_dichotomous.1_minus_", this_name, ".rds"),
                                        cv_fit_name = paste0("cvfitted_dichotomous.1_minus_", this_name, ".rds"),
                                        family = "binomial",
                                        SL.library = SL.library,
                                        cvControl = list(V = 10),
                                        method = "method.CC_nloglik",
                                        reduce_covs = reduce_covs,
                                        save_full_object = FALSE)
        }
    }
  sl_dichotomous2 <- sl_one_outcome(outcome_name = "dichotomous.2",
                            pred_names = pred_names,
                           family = "binomial",
                           SL.library = SL.library, 
                           cvControl = list(V = 10),
                         method = "method.CC_nloglik")
  for (i in 1:length(all_var_groups)) {
        if (length(all_var_groups[i]) == 0) {

        } else {
            this_name <- names(all_var_groups)[i]
            sl_dichotomous2_i <- sl_one_outcome(outcome_name = "dichotomous.2",
                                        pred_names = pred_names[!(pred_names %in% all_var_groups[[i]])],
                                        fit_name = paste0("fitted_dichotomous.2_minus_", this_name, ".rds"),
                                        cv_fit_name = paste0("cvfitted_dichotomous.2_minus_", this_name, ".rds"),
                                        family = "binomial",
                                        SL.library = SL.library,
                                        cvControl = list(V = 10),
                                        method = "method.CC_nloglik",
                                        reduce_covs = reduce_covs,
                                        save_full_object = FALSE)
        }
    }
}