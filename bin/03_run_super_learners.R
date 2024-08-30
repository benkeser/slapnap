#! /usr/bin/env -S Rscript --vanilla

#  (1) Run regression for the user-defined set of outcomes, in opts$outcomes: regression of outcome on all features, for use as the "best possible outcome predictor"
#  (2) If "cond" is in opts$importance_grp, run regression of each outcome in opts$outcomes on each reduced set of features (created by removing the pre-defined group of interest): for use in group variable importance, conditional on all other features being in the model
#  (3) If "marg" is in opts$importance_grp, run regression of each outcome in opts$outcomes on only pre-defined group + confounders: for use in group marginal feature importance
#  (4) If "marg" is in opts$importance_grp OR opts$importance_ind, run a regression of outcome on only confounders: for use in individual-level and group-level marginal feature importance
#  (5) If "cond" is in opts$importance_ind, run regression of outcome on all features except that one, for use in individual-level conditional importance
# (6) If "marg" is in opts$importance_ind, run regression of outcome on only feature of interest + geographic confounders (note that this is a glm always)
# For (5) and (6), if opts$ind_importance_type == "sitewise", aggregate all residues at a site for individual importance; otherwise, use individual residues
# ---------------------------------------------------------------------------
# Set up args, variables, functions
# ---------------------------------------------------------------------------
# load libraries
library("SuperLearner")
library("dplyr")
library("slapnap") # source our function library

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()
# path.home <- opts$output

# TODO: Use opts
path.data.analysis <- Sys.getenv("analysis") #file.path(path.data, "analysis")
path.ml.slfits <- Sys.getenv("slfits") # file.path(path.data, "analysis")

# If h2o is listed as a learner, initiate h2o cluster
h2o_here <- !(all(grepl("h2oboost", opts$learners) == FALSE))
if (h2o_here) {
    library("h2o")
    # initiate h2o cluster
    h2o.init(max_mem_size = "32G")
    # don't print the progress bar
    h2o::h2o.no_progress()
}

# load data and subset to complete cases
analysis_data_names <- list.files(path.data.analysis)
analysis_data_names <- get_analysis_dataset_name(analysis_data_names, opts)
dat <- readRDS(paste0(path.data.analysis, analysis_data_names))

# make super learner library
SL.library <- make_sl_library_vector(opts = opts)
print(SL.library)
# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]

# get names of outcomes
outcome_names <- get_outcome_names(opts)
one_nab <- length(opts$nab) == 1
all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
all_labels <- c("IC-50", "IC-80", "IIP", ifelse(one_nab, "Sensitivity", "Estimated sensitivity"), "Multiple sensitivity")
nice_outcomes <- opts$outcomes
for(i in seq_along(all_outcomes)){
    nice_outcomes <- gsub(all_outcomes[i], all_labels[i], nice_outcomes)
}

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
# get individual variables -- either sitewise or residuewise
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- get_individual_features(pred_names[!grepl("geog", pred_names)][1:num_covs], opts$ind_importance_type)

# set number of CV folds
V <- as.numeric(opts$nfolds)

# check the outcomes to see if we can run them or not
run_sl_vimp_bools <- check_outcomes(dat, outcome_names, V)
run_sl_vimp_bools2 <- lapply(run_sl_vimp_bools, function(x){
    x[c("ic50", "ic80", "iip", "sens1", "sens2") %in% opts$outcomes]
})

# print a message for any that are false (first SL, then vimp)
for (i in 1:length(outcome_names)) {
    this_outcome_name <- nice_outcomes[i]
    if (!run_sl_vimp_bools2$run_sl[i]) {
        this_class <- (0:1)[which.min(table(dat[, outcome_names[i]]))]
        print(paste0("The number of observations in class ", this_class," is less than or equal to the number of CV folds for requested outcome '", this_outcome_name, "'. The SuperLearner will not be run for this outcome."))
    }
    if (!run_sl_vimp_bools2$run_vimp[i] & any(opts$importance_grp != "" | (opts$importance_ind != "" | opts$importance_ind != "pred"))) {
        this_class <- (0:1)[which.min(table(dat[, outcome_names[i]]))]
        print(paste0("The number of observations in class ", this_class, " is too small to run population variable importance for requested outcome '", this_outcome_name, "'."))
    }
}

# ----------------------------------------------------------------------------
# (1) run full super learners for each outcome specified in outcome_names
# ----------------------------------------------------------------------------
# get the sample-splitting folds for variable importance
if (V > 1) {
    sample_splitting_folds <- vimp::make_folds(y = seq_len(V), V = 2)
}
for (i in 1:length(outcome_names)) {
    set.seed(123125)
    this_outcome_name <- outcome_names[i]
    sl_opts <- get_sl_options(this_outcome_name, V = V)
    if (V <= 1) {
        sample_splitting_folds <- vimp::make_folds(y = dat[, this_outcome_name], V = 2)
    }
    saveRDS(sample_splitting_folds, paste0(path.ml.slfits, "ss_folds_", this_outcome_name, ".rds"))
    # do the fitting, if there are enough outcomes
    if (run_sl_vimp_bools2$run_sl[i]) {
        print(paste0("Fitting ", nice_outcomes[i]))
        sl_fit_i <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name, pred_names = pred_names, family = sl_opts$fam, SL.library = SL.library,
            cvControl = sl_opts$ctrl, method = sl_opts$method, opts = opts, h2o_here = h2o_here, ss_folds = sample_splitting_folds)
    }
}

# ----------------------------------------------------------------------------
# Only run the following code if neither opts$importance_grp nor opts$importance_ind is empty
# ----------------------------------------------------------------------------
# get the individual AA SL library
ind_sl_lib <- switch(
    (grepl("residue", opts$ind_importance_type)) + 1,
    SL.library, "SL.glm"
)
# ----------------------------------------------------------------------------
# (2) If "cond" is in opts$importance_grp, run regression of each outcome in outcome_names on the reduced set of features defined by removing the group of interest
# (3) If "marg" is in opts$importance_grp, run regression of each outcome in outcome_names on the set of features defined by the group of interest + confounders
# (4) If "marg" is in opts$importance_grp, run regression of each outcome in outcome_names on geographic confounders only
# ----------------------------------------------------------------------------
if (("cond" %in% opts$importance_grp) | ("marg" %in% opts$importance_grp | "marg" %in% opts$importance_ind)) {
    for (i in 1:length(outcome_names)) {
        set.seed(4747)
        this_outcome_name <- outcome_names[i]
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        # set up validation rows for CV SuperLearner
        cross_fitting_folds <- readRDS(paste0(path.ml.slfits, "cvfolds_", this_outcome_name, ".rds"))
        sl_opts$ctrl$validRows <- cross_fitting_folds
        # only do this if we have enough obs to run it
        if (run_sl_vimp_bools2$run_vimp[i]) {
            cat("Fitting reduced learners for outcome", nice_outcomes[i], "\n")
            if ("marg" %in% opts$importance_grp | "cond" %in% opts$importance_grp) {
                for (j in 1:length(all_var_groups)) {
                    if (length(all_var_groups[j]) != 0) {
                        this_group_name <- names(all_var_groups)[j]
                        print(paste0("Fitting reduced learners for group variable importance of ", this_group_name))
                        # fit based on removing group of interest
                        if ("cond" %in% opts$importance_grp) {
                            sl_fit_ij <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name,
                                pred_names = pred_names[!(pred_names %in% all_var_groups[[j]])],
                                fit_name = paste0("fitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"),
                                cv_fit_name = paste0("cvfitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"),
                                family = sl_opts$fam, SL.library = SL.library,
                                cvControl = sl_opts$ctrl, method = sl_opts$method,
                                save_full_object = FALSE, ss_folds = sample_splitting_folds,
                                full_fit = FALSE, opts = opts, h2o_here = h2o_here)
                        }
                        # fit based on only group of interest + geographic confounders
                        if ("marg" %in% opts$importance_grp) {
                            sl_fit_marginal_ij <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name,
                                pred_names = pred_names[(pred_names %in% all_var_groups[[j]]) | (pred_names %in% all_geog_vars)],
                                fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"),
                                cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"),
                                family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                                method = sl_opts$method, save_full_object = FALSE, ss_folds = sample_splitting_folds,
                                full_fit = TRUE, opts = opts, h2o_here = h2o_here)
                        }
                    }
                }
            }
            # if "marg" is in opts$importance_grp or in opts$importance_ind and opts$ind_importance_type is sitewise, fit a regression of outcome on geographic confounders only
            if ("marg" %in% opts$importance_grp | ("marg" %in% opts$importance_ind & grepl("site", opts$ind_importance_type))) {
                sl_geog_i <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)],
                    fit_name = paste0("fitted_", this_outcome_name, "_geog.rds"),
                    cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog.rds"),
                    family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                    method = sl_opts$method,
                    save_full_object = FALSE, ss_folds = sample_splitting_folds, full_fit = FALSE,
                    opts = opts, h2o_here = h2o_here)
            }
            # if "marg" is in opts$importance_ind, fit a simple regression of outcome on geographic confounders only
            if ("marg" %in% opts$importance_ind & grepl("residue", opts$ind_importance_type)) {
                sl_geog_glm_i <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)],
                    fit_name = paste0("fitted_", this_outcome_name, "_geog_glm.rds"),
                    cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog_glm.rds"),
                    family = sl_opts$fam, SL.library = "SL.glm", cvControl = sl_opts$ctrl,
                    method = sl_opts$method,
                    save_full_object = FALSE, ss_folds = sample_splitting_folds, full_fit = FALSE,
                    opts = opts, h2o_here = h2o_here)
            }
        }
    }
}
# ----------------------------------------------------------------------------
# (5) If "cond" is in opts$importance_ind, run regressions dropping each individual feature (residuewise) or each site (a group of residues)
# (6) If "marg" is in opts$importance_ind, run regressions with each feature/site + geographic confounders (if residue-wise, this is a glm)
# ----------------------------------------------------------------------------
if (("cond" %in% opts$importance_ind) | ("marg" %in% opts$importance_ind)) {
    for (i in 1:length(outcome_names)) {
        set.seed(1234)
        this_outcome_name <- outcome_names[i]
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        # set up validation rows for CV SuperLearner
        cross_fitting_folds <- readRDS(paste0(path.ml.slfits, "cvfolds_", this_outcome_name, ".rds"))
        sl_opts$ctrl$validRows <- cross_fitting_folds
        if (run_sl_vimp_bools2$run_vimp[i]) {
            print(paste0("Fitting reduced learners for individual variable importance for outcome ", nice_outcomes[i]))
            for (j in 1:length(var_inds)) {
                this_var_name <- names(var_inds)[j]
                # if conditional, do regression of everything but this one
                if ("cond" %in% opts$importance_ind) {
                    sl_fit_ij <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name,
                        pred_names = pred_names[!(pred_names %in% var_inds[[j]]) | (pred_names %in% all_geog_vars)],
                        fit_name = paste0("fitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"),
                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"),
                        family = sl_opts$fam, SL.library = SL.library, cvControl = sl_opts$ctrl,
                        method = sl_opts$method, save_full_object = FALSE, ss_folds = sample_splitting_folds,
                        full_fit = FALSE, opts = opts, h2o_here = h2o_here)
                }
                # if marginal, do regression of this + confounders
                if ("marg" %in% opts$importance_ind) {
                    sl_fit_ij <- sl_one_outcome(complete_dat = dat, outcome_name = this_outcome_name,
                        pred_names = pred_names[(pred_names %in% var_inds[[j]]) | (pred_names %in% all_geog_vars)],
                        fit_name = paste0("fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"),
                        cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"),
                        family = sl_opts$fam, SL.library = ind_sl_lib, cvControl = sl_opts$ctrl,
                        method = sl_opts$method, save_full_object = FALSE, ss_folds = sample_splitting_folds,
                        full_fit = TRUE, opts = opts, h2o_here = h2o_here)
                }
            }
        }
    }
}

# shutdown h2o
if (h2o_here) {
    h2o.shutdown(prompt = FALSE)
}
