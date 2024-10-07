#! /usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Set up args, variables, functions
# ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("vimp")
library("dplyr")
source("/home/lib/04_variable_groups.R")
source("/home/lib/03_super_learner_libraries.R")
source("/home/lib/00_utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# load data and subset to complete cases
analysis_data_names <- list.files("/home/dat/analysis")
analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
complete_dat <- readRDS(paste0("/home/dat/analysis/", analysis_data_name))

# make super learner library
SL.library <- make_sl_library_vector(opts = opts)

# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(complete_dat))) # geography seems to be first column of relevant data
pred_names <- colnames(complete_dat)[geog_idx:ncol(complete_dat)]

# get names of outcomes
outcome_names <- get_outcome_names(opts)

# get variable groups
all_var_groups <- get_variable_groups(complete_dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
# get individual variables -- either sitewise or residuewise
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- get_individual_features(pred_names[!grepl("geog", pred_names)][1:num_covs], opts$ind_importance_type)

# set number of CV folds
V <- as.numeric(opts$nfolds)

# check the outcomes to see if we can run them or not
run_sl_vimp_bools <- check_outcomes(complete_dat, outcome_names, V)
run_sl_vimp_bools2 <- lapply(run_sl_vimp_bools, function(x){
    x[c("ic50", "ic80", "iip", "sens1", "sens2") %in% opts$outcomes]
})
# ---------------------------------------------------------------------------
# get variable importance! but only run if one of opts$importance_grp or opts$importance_ind is not empty
# ---------------------------------------------------------------------------
# logical: should we estimated VIM using cross-fitting?
cf_vim <- (length(opts$learners) == 1 & opts$cvtune & opts$cvperf) | (length(opts$learners) > 1 & opts$cvperf)
set.seed(474747)
# if none of them, then don't run variable importance
if (((length(opts$importance_grp) == 0) & (length(opts$importance_ind) == 0))) {
    print("Variable importance was not requested during this run. If you desire variable importance, please use the environment variables 'importance_grp' (for group importance) or 'importance_ind' (for individual-variable importance).")
} else { # otherwise, do run variable importance
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        vimp_opts <- get_vimp_options(this_outcome_name)
        ss_folds <- readRDS(paste0("/home/slfits/ss_folds_", this_outcome_name, ".rds"))
        # subset to the complete data for this outcome (or all outcomes, if specified)
        if(!opts$same_subset | !(("ic80" %in% opts$outcomes | "iip" %in% opts$outcomes) & length(opts$outcomes) > 1)){
          complete_cases_idx <- complete.cases(complete_dat[, c(this_outcome_name, pred_names)])
        } else {
          complete_cases_idx <- complete.cases(complete_dat)
        }
        dat <- complete_dat[complete_cases_idx, ]
        y <- dat[, this_outcome_name]
        if (run_sl_vimp_bools2$run_vimp[i]) {
            # create output list
            eval(parse(text = paste0(this_outcome_name, '_vimp_lst <- make_vimp_list(all_var_groups, var_inds)')))
            # read in the full fit
            if (cf_vim) {
                full_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, ".rds"))
                full_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, ".rds"))
                cross_fitting_folds <- vimp::get_cv_sl_folds(full_folds)
            } else {
                full_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, ".rds"))
            }
            # if "marg" is in opts$importance_grp AND/OR using sitewise individual importance, read in the fit with geographic confounders only
            if (("marg" %in% opts$importance_grp) | 
                ("marg" %in% opts$importance_ind & grepl("site", opts$ind_importance_type))) {
                # load geog-only fit corresponding to this outcome
                if (cf_vim) {
                    geog_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_geog.rds"))
                } else {
                    geog_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_geog.rds"))
                }
            }
            # if "marg" is in opts$importance_ind AND using residuewise individual importance, read in the simple fit with geographic confounders only
            if (("marg" %in% opts$importance_ind) & grepl("residue", opts$ind_importance_type)) {
                # load simple geog-only fit corresponding to this outcome
                if (cf_vim) {
                    geog_glm_preds <-readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_geog_glm.rds"))
                } else {
                    geog_glm_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_geog_glm.rds"))
                }
            }
            # -----------------------------------
            # group variable importance
            # -----------------------------------
            # if "cond" is in opts$importance_grp, run this loop
            if ("cond" %in% opts$importance_grp) {
                for (j in 1:length(all_var_groups)) {
                    this_group_name <- names(all_var_groups)[j]
                    if (cf_vim) {
                        cond_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, 
                                                     "_conditional_", this_group_name, ".rds"))
                        # get conditional, CV-VIM
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_group_name, 
                                                                  " <- vimp::cv_vim(Y = y, ",
                                                                  "cross_fitted_f1 = full_preds, ",
                                                                  "cross_fitted_f2 = cond_preds, ",
                                                                  "indx = which(pred_names %in% all_var_groups[[j]]), ",
                                                                  "run_regression = FALSE, sample_splitting = TRUE, ",
                                                                  "alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, ",
                                                                  "cross_fitting_folds = cross_fitting_folds, ",
                                                                  "sample_splitting_folds = ss_folds, V = V / 2, ",
                                                                  "na.rm = TRUE, scale = 'identity')"))))
                    } else {
                        cond_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_conditional_", 
                                                     this_group_name, ".rds"))
                        # get conditional, non-CV VIM
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_group_name, 
                                                                  " <- vimp::vim(Y = y, ",
                                                                  "f1 = full_preds, ",
                                                                  "f2 = cond_preds, ",
                                                                  "indx = which(pred_names %in% all_var_groups[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type = vimp_opts$vimp_measure, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                # merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$grp_conditional <- merge_vim(", paste(paste0(this_outcome_name, 
                                         "_cond_", names(all_var_groups)), collapse = ", "), ")")))
            }
            # if "marg" is in opts$importance_grp, run this loop
            if ("marg" %in% opts$importance_grp) {
                for (j in 1:length(all_var_groups)) {
                    this_group_name <- names(all_var_groups)[j]
                    if (cf_vim) {
                        # get marginal, CV-VIM
                        marg_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_marginal_", 
                                                     this_group_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_group_name, 
                                                                  " <- vimp::cv_vim(Y = y, ",
                                                                  "cross_fitted_f1 = marg_preds, ",
                                                                  "cross_fitted_f2 = geog_preds, ",
                                                                  "indx = which(pred_names %in% all_var_groups[[j]]), ",
                                                                  "run_regression = FALSE, sample_splitting = TRUE, ",
                                                                  "alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, ",
                                                                  "cross_fitting_folds = cross_fitting_folds, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "V = V / 2, na.rm = TRUE, scale = 'identity')"))))
                    } else {
                        # get marginal, non-CV VIM
                        marg_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_marginal_", 
                                                     this_group_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_group_name, 
                                                                  " <- vimp::vim(Y = y, ",
                                                                  "f1 = marg_preds, ",
                                                                  "f2 = geog_preds, ",
                                                                  "indx = which(pred_names %in% all_var_groups[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type = vimp_opts$vimp_measure, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                # merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$grp_marginal <- merge_vim(", 
                                         paste(paste0(this_outcome_name, "_marg_", names(all_var_groups)), collapse = ", "), ")")))
            }
            # -----------------------------------
            # individual variable importance
            # -----------------------------------
            # if "cond" is in opts$importance_ind, run this loop
            if ("cond" %in% opts$importance_ind) {
                for (j in 1:length(var_inds)) {
                    this_var_name <- names(var_inds)[j]
                    if (cf_vim) {
                        # get individual, CV-VIM (conditional)
                        indi_cond_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, 
                                                          "_conditional_", this_var_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_var_name, 
                                                                  " <- vimp::cv_vim(Y = y, ",
                                                                  "cross_fitted_f1 = full_preds, ",
                                                                  "cross_fitted_f2 = indi_cond_preds, ",
                                                                  "indx = which(pred_names %in% var_inds[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type = vimp_opts$vimp_measure, ",
                                                                  "cross_fitting_folds = cross_fitting_folds, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "V = V / 2, na.rm = TRUE, scale = 'identity')"))))
                    } else {
                        # get individual, non-CV VIM (conditional)
                        indi_cond_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, 
                                                          "_conditional_", this_var_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_var_name, 
                                                                  " <- vimp::vim(Y = y, ",
                                                                  "f1 = full_preds, ",
                                                                  "f2 = indi_cond_preds, ",
                                                                  "indx = which(pred_names %in% var_inds[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type = vimp_opts$vimp_measure, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                # merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$ind_conditional <- merge_vim(", 
                                         paste(paste0(this_outcome_name, "_cond_", names(var_inds)), collapse = ", "), ")")))
            }
            # if "marg" is in opts$importance_ind, run this loop
            if ("marg" %in% opts$importance_ind) {
                if (grepl("residue", opts$ind_importance_type)) {
                    reduced_preds <- geog_glm_preds
                } else {
                    reduced_preds <- geog_preds
                }
                for (j in 1:length(var_inds)) {
                    this_var_name <- names(var_inds)[j]
                    if (cf_vim) {
                        # get individual, CV-VIM (marginal)
                        indi_marg_preds <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, 
                                                          "_marginal_", this_var_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_var_name, 
                                                                  " <- vimp::cv_vim(Y = y, ",
                                                                  "cross_fitted_f1 = indi_marg_preds, ", 
                                                                  "cross_fitted_f2 = reduced_preds, ",
                                                                  "indx = which(pred_names %in% var_inds[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type =  vimp_opts$vimp_measure, ",
                                                                  "cross_fitting_folds = cross_fitting_folds, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "V = V / 2, na.rm = TRUE, scale = 'identity')"))))
                    } else {
                        # get individual, non-CV VIM (marginal)
                        indi_marg_preds <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, 
                                                          "_marginal_", this_var_name, ".rds"))
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_var_name, 
                                                                  " <- vimp::vim(Y = y, ",
                                                                  "f1 = indi_marg_preds, ",
                                                                  "f2 = reduced_preds, ",
                                                                  "indx = which(pred_names %in% var_inds[[j]]), ",
                                                                  "run_regression = FALSE, alpha = 0.05, delta = 0, ",
                                                                  "type = vimp_opts$vimp_measure, ",
                                                                  "sample_splitting_folds = ss_folds, ",
                                                                  "na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                # merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$ind_marginal <- merge_vim(", 
                                         paste(paste0(this_outcome_name, "_marg_", names(var_inds)), collapse = ", "), ")")))
            }
            # save them off
            eval(parse(text = paste0("saveRDS(", this_outcome_name, "_vimp_lst, file = '/home/slfits/", 
                                     paste0(this_outcome_name, "_vimp"), ".rds')")))
        }
    }
}
