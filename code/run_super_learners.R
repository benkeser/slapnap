#! /usr/bin/env Rscript

#  (1) Run regression for the user-defined set of outcomes, in opts$outcomes: regression of outcome on all features, for use as the "best possible outcome predictor"
#  (2) If "cond" is in opts$importance_grp, run regression of each outcome in opts$outcomes on each reduced set of features (created by removing the pre-defined group of interest): for use in group variable importance, conditional on all other features being in the model
#  (3) If "marg" is in opts$importance_grp, run regression of each outcome in opts$outcomes on only pre-defined group + confounders: for use in group marginal feature importance
#  (4) If "marg" is in opts$importance_grp OR opts$importance_ind, run a regression of outcome on only confounders: for use in individual-level and group-level marginal feature importance
#  (5) If "cond" is in opts$importance_ind, run regression of outcome on all features except that one, for use in individual-level conditional importance
# (6) If "marg" is in opts$importance_ind, run regression of outcome on only feature of interest + geographic confounders (note that this is a glm always)
## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------
# load libraries
library("SuperLearner")
library("dplyr")
source("/home/lib/variable_groups.R")
source("/home/lib/super_learner_libraries.R")
source("/home/lib/utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# load data and subset to complete cases
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
nprevious <- length(dat[,1])
saveRDS(nprevious, "/home/slfits/nprevious.rds")
dat <- dat[complete.cases(dat),]

# save for report compilation later
ncomplete <- length(dat[,1])
saveRDS(ncomplete, "/home/slfits/ncomplete.rds")

# make super learner library
SL.library <- make_sl_library_vector(opts = opts)

# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]

# get names of outcomes
outcome_names <- c(
    switch("ic50" %in% opts$outcomes, "log10.pc.ic50", NULL),
    switch("ic80" %in% opts$outcomes, "log10.pc.ic80", NULL),
    switch("iip" %in% opts$outcomes, "iip", NULL),
    switch("sens1" %in% opts$outcomes, "dichotomous.1", NULL),
    switch("sens2" %in% opts$outcomes, "dichotomous.2", NULL)
)

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]

# set number of CV folds
V <- 5 # move to user control?

set.seed(123125)
## ----------------------------------------------------------------------------
## (1) run full super learners for each outcome specified in outcome_names
## ----------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name, V = V)
    ## do the fitting
    sl_fit_i <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name, pred_names = pred_names, family = sl_opts$fam, SL.library = SL.library,
        cvControl = sl_opts$ctrl, method = sl_opts$method, opts = opts)
    ## if we need any type of importance, generate splits for VIM hypothesis testing
    if (!((length(opts$importance_grp) == 0) & (length(opts$importance_ind) == 0))) {
        outer_folds <- make_folds(dat[, this_outcome_name], V = 2, stratified = grepl("dichot", this_outcome_name))
        saveRDS(outer_folds, file = paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
    }
    ## if conditional importance is desired, fit the full regression
    if (("cond" %in% opts$importance_grp) | ("cond" %in% opts$importance_ind)) {
        sl_split_fit_i <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name, pred_names = pred_names,
            fit_name = paste0("fitted_", this_outcome_name, "_for_vimp.rds"), cv_fit_name = paste0("cvfitted_", this_outcome_name, "_for_vimp.rds"),
            family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl, method = sl_opts$method, outer_folds = outer_folds, full_fit = TRUE, opts = opts)
    }
}

## ----------------------------------------------------------------------------
## Only run the following code if neither opts$importance_grp nor opts$importance_ind is empty
## ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# (2) If "cond" is in opts$importance_grp, run regression of each outcome in outcome_names on the reduced set of features defined by removing the group of interest
# (3) If "marg" is in opts$importance_grp, run regression of each outcome in outcome_names on the set of features defined by the group of interest + confounders
# (4) If "marg" is in opts$importance_grp, run regression of each outcome in outcome_names on geographic confounders only
# ----------------------------------------------------------------------------
if (("cond" %in% opts$importance_grp) | ("marg" %in% opts$importance_grp)) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        ## read in the full folds, for this group's cv folds
        outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
        for (j in 1:length(all_var_groups)) {
            if (length(all_var_groups[j]) != 0) {
                this_group_name <- names(all_var_groups)[j]
                ## fit based on removing group of interest
                if ("cond" %in% opts$importance_grp) {
                    sl_fit_ij <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name,
                        pred_names = pred_names[!(pred_names %in% all_var_groups[[j]])],
                        fit_name = paste0("fitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"),
                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"),
                        family = sl_opts$fam, SL.library = SL.library,
                        cvControl = sl_opts$ctrl, method = sl_opts$method,
                        save_full_object = FALSE, outer_folds = outer_folds,
                        full_fit = FALSE, opts = opts)
                }
                ## fit based on only group of interest + geographic confounders
                if ("marg" %in% opts$importance_grp) {
                    sl_fit_marginal_ij <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name,
                        pred_names = pred_names[(pred_names %in% all_var_groups[[j]]) | (pred_names %in% all_geog_vars)],
                        fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"),
                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"),
                        family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                        method = sl_opts$method, save_full_object = FALSE,
                        outer_folds = outer_folds, full_fit = FALSE, opts = opts)
                }
            }
        }
        ## if "marg" is in opts$importance_grp, fit a regression of outcome on geographic confounders only
        if ("marg" %in% opts$importance_grp) {
            sl_geog_i <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)],
                fit_name = paste0("fitted_", this_outcome_name, "_geog.rds"),
                cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog.rds"),
                family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                method = sl_opts$method,
                save_full_object = FALSE, outer_folds = outer_folds, full_fit = TRUE,
                opts = opts)
        }
        ## if "marg" is in opts$importance_ind, fit a glm regression of outcome on geographic confounders only
        if ("marg" %in% opts$importance_ind) {
            sl_geog_glm_i <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)],
                fit_name = paste0("fitted_", this_outcome_name, "_geog_glm.rds"),
                cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog_glm.rds"),
                family = sl_opts$fam, SL.library = "SL.glm", cvControl = sl_opts$ctrl,
                method = sl_opts$method,
                save_full_object = FALSE, outer_folds = outer_folds, full_fit = TRUE,
                opts = opts)
        }
    }
}
# ----------------------------------------------------------------------------
# (5) If "cond" is in opts$importance_ind, run regressions dropping each individual feature
# (6) If "marg" is in opts$importance_ind, run regressions with each feature + geographic confounders (note that for all variables, this is a glm)
# ----------------------------------------------------------------------------
if (("cond" %in% opts$importance_ind) | ("marg" %in% opts$importance_ind)) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
        for (j in 1:length(var_inds)) {
            this_var_name <- var_inds[j]
            # if conditional, do regression of everything but this one
            if ("cond" %in% opts$importance_ind) {
                sl_fit_ij <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name,
                    pred_names = pred_names[!(pred_names %in% var_inds[j]) | (pred_names %in% all_geog_vars)],
                    fit_name = paste0("fitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"),
                    cv_fit_name = paste0("cvfitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"),
                    family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                    method = sl_opts$method,  save_full_object = FALSE,
                    outer_folds = outer_folds, full_fit = FALSE, opts = opts)
            }
            # if marginal, do glm of this + confounders
            if ("marg" %in% opts$importance_ind) {
                sl_fit_ij <- sl_one_outcome(dat = dat, outcome_name = this_outcome_name,
                    pred_names = pred_names[(pred_names %in% var_inds[j]) | (pred_names %in% all_geog_vars)],
                    fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"),
                    cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"),
                    family = sl_opts$fam, SL.library = "SL.glm", cvControl = sl_opts$ctrl,
                    method = sl_opts$method,  save_full_object = FALSE,
                    outer_folds = outer_folds, full_fit = FALSE, opts = opts)
            }
        }
    }
}
