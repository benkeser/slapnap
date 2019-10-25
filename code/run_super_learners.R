#! /usr/bin/env Rscript

## run the super learners for all outcomes, pre-defined variable groups, and individual variables
## this involves:
##  (1) a regression of outcome on all features: for use as the "best possible outcome predictor"
##  (2) a regression of outcome on reduced set of features (created by removing the pre-defined group of interest): for use in group variable importance, conditional on all other features being in the model
##  (3) a regression of outcome on only pre-defined group + confounders: for use in group marginal feature importance
##  (4) a regression of outcome on each individual feature + confounders: for use in individual-level marginal feature importance
##  (5) a regression of outcome on only confounders: for use in individual-level marginal feature importance

## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------
# load libraries
library(SuperLearner)
source("/home/lib/variable_groups.R")
source("/home/lib/super_learner_libraries.R")

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"

# ~~~!~~~! add option to not do cross-val sl fitting !~~~!~~~ #
no_cv <- Sys.getenv("no_cv") == "TRUE"

# load data
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
dat <- dat[complete.cases(dat),]

# may need at some point to source in code to do super learner
# library preparation etc...
# source()

# if short run, use a simple library
if(reduce_library){
  SL.library <- default_library_reduced # defined in super_learner_libraries.R
}else{
  SL.library <- default_library # defined in super_learner_libraries.R
}

# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]
# get names of outcomes
all_outcome_names <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")
# if reduce_outcomes, only run on ic50
if (reduce_outcomes) {
    outcome_names <- "log10.pc.ic50"
} else {
    outcome_names <- all_outcome_names
}
# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
# if reduce_groups, only run on the cd4 binding site
if (reduce_groups) {
    var_groups <- all_var_groups[1]
} else {
    var_groups <- all_var_groups
}
V <- 10
# determine SL options based on outcome name
get_sl_options <- function(outcome_name) {
    if (grepl("dichot", outcome_name)) {
        sl_fam <- "binomial"
        cv_ctrl_lst <- list(V = V, stratifyCV = TRUE)
        sl_method <- "tmp_method.CC_nloglik"
    } else {
        sl_fam <- "gaussian"
        cv_ctrl_lst <- list(V = V)
        sl_method <- "tmp_method.CC_LS"
    }
    return(list(fam = sl_fam, ctrl = cv_ctrl_lst, method = sl_method))
}
# if reduce_covs, only do individual imp on that number
if (reduce_covs) {
    num_covs <- 10
    var_inds <- pred_names[!grepl("geog", pred_names)][1:(num_covs - length(all_geog_vars))]
} else {
    num_covs <- length(pred_names) - length(all_geog_vars)
    var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]
}

set.seed(123125)

#' function to run super learner and cv super learner on a single outcome
#' @param outcome_name String name of outcome
#' @param pred_names Vector of string names of predictor variables
#' @param save_dir name of directory to save results to
#' @param fit_name name of fits (defaults to fit_<outcome_name>.rds)
#' @param cv_fit_name name of CV fits (defaults to cvfit_<outcome_name>.rds)
#' @param reduce_covs Flag to reduce the number of covariates under consideration
#' @param save_full_object Flag for whether or not to save the full fitted object, or just the fitted values
#' @param run_cv Whether or not to run the cv super learner
sl_one_outcome <- function(outcome_name,
                           pred_names,
                           save_dir = "/home/slfits/",
                           fit_name = paste0("fit_", outcome_name, ".rds"),
                           cv_fit_name = paste0("cvfit_", outcome_name, ".rds"),
                           reduce_covs = FALSE,
                           save_full_object = TRUE,
                           run_cv = TRUE,
                           ...){
        pred <- dat[ , pred_names]

        if(reduce_covs){
          pred <- pred[ , 1:num_covs]
        }

        fit <- SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
        if (save_full_object) {
            saveRDS(fit, file = paste0(save_dir, fit_name))
        }
        # save super learner predictions
        saveRDS(fit$SL.predict, file = paste0(save_dir, gsub(".RData", ".rds", gsub("fit_", "fitted_", fit_name))))
        # save super learner weights
        saveRDS(fit$coef, file = paste0(save_dir, "slweights_", fit_name))

        if(run_cv){
          cv_fit <- CV.SuperLearner(Y = dat[ , outcome_name], X = pred, ...)
          if (save_full_object) {
              saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
          }
          saveRDS(cv_fit$SL.predict, file = paste0(save_dir, gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name))))
          saveRDS(cv_fit$folds, file = paste0(save_dir, gsub("cvfitted_", "cvfolds_", gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name)))))
        }
        return(invisible(NULL))
}

## ----------------------------------------------------------------------------
## (1) run full super learners for each outcome (unless reduce_outcomes = TRUE)
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name)
    sl_fit_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names, family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, reduce_covs = reduce_covs, run_cv = !no_cv)
}
## ----------------------------------------------------------------------------
## (2)+(3) run super learners for each outcome (unless reduce_outcomes = TRUE)
##         on (2) reduced set of features defined by removing group of interest
##         and (3) group of interest + geographic confounders
## ----------------------------------------------------------------------------
## run super learners on pre-defined groups
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name)
    for (j in 1:length(var_groups)) {
        if (length(var_groups[j]) != 0) {
            this_group_name <- names(var_groups)[j]
            ## fit based on removing group of interest
            sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[!(pred_names %in% var_groups[[j]])], fit_name = paste0("fitted_", this_outcome_name, "_minus_", this_group_name, ".rds"), cv_fit_name = paste0("cvfitted_", this_outcome_name, "_minus_", this_group_name, ".rds"), family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, reduce_covs = reduce_covs, run_cv = !no_cv, save_full_object = FALSE)
            ## fit based on only group of interest + geographic confounders
            sl_fit_marginal_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% var_groups[[j]]) | (pred_names %in% all_geog_vars)], fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"), cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"), family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, reduce_covs = reduce_covs, run_cv = !no_cv, save_full_object = FALSE)
        }
    }
}
## ----------------------------------------------------------------------------
## (4) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on each individual feature + confounders
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name)
    for (j in 1:length(var_inds)) {
        this_var_name <- var_inds[j]
        ## fit SL of only this variable plus geographic confounders
        sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% var_inds[j]) | (pred_names %in% all_geog_vars)], fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"), cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"), family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, reduce_covs = FALSE, run_cv = !no_cv, save_full_object = FALSE)
    }
}

## ----------------------------------------------------------------------------
## (5) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on only confounders
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name)
    sl_geog_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)], fit_name = paste0("fitted_", this_outcome_name, "_geog.rds"), cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog.rds"), family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, reduce_covs = FALSE, run_cv = !no_cv, save_full_object = FALSE)
}
