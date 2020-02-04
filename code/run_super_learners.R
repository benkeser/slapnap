#! /usr/bin/env Rscript

# ~DB2BW: we should change this description when all is said and done. 
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
source("/home/lib/utils.R")

#----------------- 
# Temp options
#----------------- 
#~DB~ delete eventually
# read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# ~DB2BW: delete?~
## option to run individual-level vimp; defaults to FALSE
# run_indi_vimp <- Sys.getenv("run_indi_vimp") == "TRUE"

# load data and subset to complete cases
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
dat <- dat[complete.cases(dat),]

# make super learner library
SL.library <- make_sl_library_vector(opts = opts)


# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]


# get names of outcomes
outcome_names <- c(
  ifelse("ic50" %in% opts$outcomes, "log10.pc.ic50", NULL),
  ifelse("ic80" %in% opts$outcomes, "log10.pc.ic80", NULL),
  ifelse("iip" %in% opts$outcomes, "iip", NULL),
  ifelse("sens1" %in% opts$outcomes, "dichotomous.1", NULL),
  ifelse("sens2" %in% opts$outcomes, "dichotomous.2", NULL)
)

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]

# set number of CV folds
V <- 5

set.seed(123125)
## ----------------------------------------------------------------------------
## (1) run full super learners for each outcome (unless reduce_outcomes = TRUE)
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name, V = V)
    ## generate a set of outer folds for sample splitting for VIM hypothesis testing
    outer_folds <- make_folds(dat[, this_outcome_name], V = 2, stratified = grepl("dichot", this_outcome_name))
    saveRDS(outer_folds, file = paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
    ## do the fitting
    sl_fit_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names, 
                               family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, 
                               method = sl_opts$method, reduce_covs = reduce_covs, opts = opts)
    #~DB2BW: Not sure where this gets used to know whether to flag it based on e.g., importance_grp~
    sl_split_fit_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names, 
                                     family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, 
                                     method = sl_opts$method, reduce_covs = reduce_covs, 
                                     outer_folds = outer_folds, full_fit = TRUE, opts = opts)
}

#~DB2BW: I'm going to let you add appropriate flags here. You can use the opts object 
#        above, which is a named list. See Dockerfile for valid values for the flags pertaining
#        to variable importance. It might be worth double checking the new sl_one_outcome function
#        to make sure everything that you need is saved. 

## ----------------------------------------------------------------------------
## (2)+(3) run super learners for each outcome (unless reduce_outcomes = TRUE)
##         on (2) reduced set of features defined by removing group of interest
##         and (3) group of interest + geographic confounders
## ----------------------------------------------------------------------------
# here are conditional variable importance fits
## run super learners on pre-defined groups
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name, V = V)
    ## read in the full folds, for this group's cv folds
    outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
    for (j in 1:length(var_groups)) {
        if (length(var_groups[j]) != 0) {
            this_group_name <- names(var_groups)[j]
            ## fit based on removing group of interest
            sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, 
                                        pred_names = pred_names[!(pred_names %in% var_groups[[j]])], 
                                        fit_name = paste0("fitted_", this_outcome_name, "_minus_", this_group_name, ".rds"), 
                                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_minus_", this_group_name, ".rds"), 
                                        family = sl_opts$fam, SL.library = SL.library, 
                                        cvControl = sl_opts$ctrl, method = sl_opts$method, 
                                        save_full_object = FALSE, outer_folds = outer_folds, 
                                        full_fit = FALSE, opts = opts)
            ## fit based on only group of interest + geographic confounders
            sl_fit_marginal_ij <- sl_one_outcome(outcome_name = this_outcome_name, 
                                                 pred_names = pred_names[(pred_names %in% var_groups[[j]]) | (pred_names %in% all_geog_vars)], 
                                                 fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"), 
                                                 cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"), 
                                                 family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, 
                                                 method = sl_opts$method, save_full_object = FALSE, 
                                                 outer_folds = outer_folds, full_fit = FALSE, opts = opts)
        }
    }
}
## ----------------------------------------------------------------------------
## (4) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on each individual feature + confounders
## ----------------------------------------------------------------------------
# here are marginal fits
if (run_indi_vimp) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        sl_opts <- get_sl_options(this_outcome_name)
        outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
        for (j in 1:length(var_inds)) {
            this_var_name <- var_inds[j]
            ## fit SL of only this variable plus geographic confounders
            sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, 
                                        pred_names = pred_names[(pred_names %in% var_inds[j]) | (pred_names %in% all_geog_vars)], 
                                        fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"), 
                                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"), 
                                        family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, 
                                        method = sl_opts$method, reduce_covs = FALSE, run_cv = !no_cv, save_full_object = FALSE, 
                                        outer_folds = outer_folds, full_fit = FALSE, opts = opts)
        }
    }
}
## ----------------------------------------------------------------------------
## (5) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on only confounders
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name)
    outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
    sl_geog_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)], 
                                fit_name = paste0("fitted_", this_outcome_name, "_geog.rds"), 
                                cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog.rds"), 
                                family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, 
                                method = sl_opts$method, reduce_covs = FALSE, run_cv = !no_cv, 
                                save_full_object = FALSE, outer_folds = outer_folds, full_fit = TRUE, 
                                opts = opts)
}
