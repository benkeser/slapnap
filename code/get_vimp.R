#! /usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("vimp")
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
analysis_data_names <- list.files("/home/dat/analysis")
analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
complete_dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)

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
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]

# set number of CV folds
V <- as.numeric(opts$nfolds)

# check the outcomes to see if we can run them or not
run_sl_vimp_bools <- check_outcomes(complete_dat, outcome_names, V)
run_sl_vimp_bools2 <- lapply(run_sl_vimp_bools, function(x){
    x[c("ic50", "ic80", "iip", "sens1", "sens2") %in% opts$outcomes]
})
## ---------------------------------------------------------------------------
## get variable importance! but only run if one of opts$importance_grp or opts$importance_ind is not empty
## ---------------------------------------------------------------------------
set.seed(474747)
# if none of them, then don't run variable importance
if (((length(opts$importance_grp) == 0) & (length(opts$importance_ind) == 0))) {
    print("Variable importance was not requested during this run. If you desire variable importance, please use the environment variables 'importance_grp' (for group importance) or 'importance_ind' (for individual-variable importance).")
} else { # otherwise, do run variable importance
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        vimp_opts <- get_vimp_options(this_outcome_name)
        # subset to the complete data for this outcome (or all outcomes, if specified)
        if(!opts$same_subset | !(("ic80" %in% opts$outcomes | "iip" %in% opts$outcomes) & length(opts$outcomes) > 1)){
          complete_cases_idx <- complete.cases(complete_dat[, c(this_outcome_name, pred_names)])
        } else {
          complete_cases_idx <- complete.cases(complete_dat)
        }
        dat <- complete_dat[complete_cases_idx, ]
        if (run_sl_vimp_bools2$run_vimp[i]) {
            ## create output list
            eval(parse(text = paste0(this_outcome_name, '_vimp_lst <- make_vimp_list(all_var_groups, var_inds)')))
            eval(parse(text = paste0(this_outcome_name, '_cv_vimp_lst <- make_vimp_list(all_var_groups, var_inds)')))
            ## load in outer folds for VIM
            outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
            ## if "cond" is in either opts$importance_grp or opts$importance_ind, read it in
            if (("cond" %in% opts$importance_grp) | ("cond" %in% opts$importance_ind)) {
                ## load full fit corresponding to this outcome
                full_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_for_vimp.rds"))
                if (opts$cvperf) {
                    full_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_for_vimp.rds"))
                    full_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_for_vimp.rds"))
                    full_cv_folds_vec <- get_cv_folds(full_cv_folds)
                    full_cv_fit_lst <- lapply(as.list(1:length(unique(full_cv_folds_vec))), function(x) full_cv_fit[full_cv_folds_vec == x])
                }
            }
            # if "marg" is in opts$importance_grp, read in the "sl" fit with geographic confounders only
            if (("marg" %in% opts$importance_grp)) {
                ## load geog-only fit corresponding to this outcome
                geog_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_geog.rds"))
                if (opts$cvperf) {
                    geog_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_geog.rds"))
                    geog_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_geog.rds"))
                    geog_cv_folds_vec <- get_cv_folds(geog_cv_folds)
                    geog_cv_fit_lst <- lapply(as.list(1:length(unique(geog_cv_folds_vec))), function(x) geog_cv_fit[geog_cv_folds_vec == x])
                }
            }
            # if "marg" is in opts$importance_ind, read in the glm fit with geographic confounders only
            if (("marg" %in% opts$importance_ind)) {
                ## load geog-only fit corresponding to this outcome
                geog_glm_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_geog_glm.rds"))
                if (opts$cvperf) {
                    geog_glm_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_geog_glm.rds"))
                    geog_glm_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_geog_glm.rds"))
                    geog_glm_cv_folds_vec <- get_cv_folds(geog_cv_folds)
                    geog_glm_cv_fit_lst <- lapply(as.list(1:length(unique(geog_glm_cv_folds_vec))), function(x) geog_glm_cv_fit[geog_glm_cv_folds_vec == x])
                }
            }
            ## -----------------------------------
            ## group variable importance
            ## -----------------------------------
            ## if "cond" is in opts$importance_grp, run this loop
            if ("cond" %in% opts$importance_grp) {
                for (j in 1:length(all_var_groups)) {
                    this_group_name <- names(all_var_groups)[j]
                    cond_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"))
                    cond_folds <- list(outer_folds = outer_folds)
                    ## get conditional, non-cv vimp
                    suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_group_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = full_fit, f2 = cond_fit, indx = which(pred_names %in% all_var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = cond_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
                    if (opts$cvperf) {
                        cond_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_conditional_", this_group_name, ".rds"))
                        cond_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_conditional_", this_group_name, ".rds"))
                        cond_cv_folds_vec <- get_cv_folds(cond_cv_folds)
                        cond_cv_fit_lst <- lapply(as.list(1:length(unique(cond_cv_folds_vec))), function(x) cond_cv_fit[cond_cv_folds_vec == x])
                        cond_folds <- list(outer_folds = outer_folds, inner_folds = list(inner_folds_1 = full_cv_folds_vec, inner_folds_2 = cond_cv_folds_vec))
                        ## get conditional, cv vimp
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_cond_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = full_cv_fit_lst, f2 = cond_cv_fit_lst, indx = which(pred_names %in% all_var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = cond_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                ## merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$grp_conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cond_", names(all_var_groups)), collapse = ", "), ")")))
                if (opts$cvperf) {
                    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$grp_conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cv_cond_", names(all_var_groups)), collapse = ", "), ")")))
                }
            }
            ## if "marg" is in opts$importance_grp, run this loop
            if ("marg" %in% opts$importance_grp) {
                for (j in 1:length(all_var_groups)) {
                    this_group_name <- names(all_var_groups)[j]
                    marg_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
                    marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2)
                    ## get marginal, non-cv vimp
                    suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_group_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = marg_fit, f2 = geog_fit, indx = which(pred_names %in% all_var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = marg_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
                    if (opts$cvperf) {
                        marg_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
                        marg_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
                        marg_cv_folds_vec <- get_cv_folds(marg_cv_folds)
                        marg_cv_fit_lst <- lapply(as.list(1:length(unique(marg_cv_folds_vec))), function(x) marg_cv_fit[marg_cv_folds_vec == x])
                        marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = marg_cv_folds_vec, inner_folds_2 = geog_cv_folds_vec))
                        ## get marginal, cv vimp
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = marg_cv_fit_lst, f2 = geog_cv_fit_lst, indx = which(pred_names %in% all_var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = marg_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                ## merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$grp_marginal <- merge_vim(", paste(paste0(this_outcome_name, "_marg_", names(all_var_groups)), collapse = ", "), ")")))
                if (opts$cvperf) {
                    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$grp_marginal <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", names(all_var_groups)), collapse = ", "), ")")))
                }
            }
            ## -----------------------------------
            ## individual variable importance
            ## -----------------------------------
            ## if "cond" is in opts$importance_ind, run this loop
            if ("cond" %in% opts$importance_ind) {
                for (j in 1:length(var_inds)) {
                    this_var_name <- var_inds[j]
                    indi_cond_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"))
                    indi_cond_folds <- list(outer_folds = outer_folds)
                    ## get individual, non-cv vimp
                    suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_var_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = full_fit, f2 = indi_cond_fit, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_cond_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
                    if (opts$cvperf) {
                        indi_cv_cond_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_conditional_", this_var_name, ".rds"))
                        indi_cv_cond_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_conditional_", this_var_name, ".rds"))
                        indi_cv_cond_folds_vec <- get_cv_folds(indi_cv_cond_folds)
                        indi_cv_cond_fit_lst <- lapply(as.list(1:length(unique(indi_cv_cond_folds_vec))), function(x) indi_cv_cond_fit[indi_cv_cond_folds_vec == x])
                        indi_cond_folds <- list(outer_folds = outer_folds, inner_folds = list(inner_folds_1 = full_cv_folds_vec, inner_folds_2 = indi_cv_cond_folds_vec))
                        ## get individual, cv vimp
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_cond_", this_var_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = full_cv_fit_lst, f2 = indi_cv_cond_fit_lst, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_cond_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                ## merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$ind_conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cond_", var_inds), collapse = ", "), ")")))
                if (opts$cvperf) {
                    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$ind_conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cv_cond_", var_inds), collapse = ", "), ")")))
                }
            }
            ## if "marg" is in opts$importance_ind, run this loop
            if ("marg" %in% opts$importance_ind) {
                for (j in 1:length(var_inds)) {
                    this_var_name <- var_inds[j]
                    indi_marg_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
                    indi_marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2)
                    ## get individual, non-cv vimp
                    suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_var_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = indi_marg_fit, f2 = geog_glm_fit, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_marg_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
                    if (opts$cvperf) {
                        ## cv
                        indi_cv_marg_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
                        indi_cv_marg_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
                        indi_cv_marg_folds_vec <- get_cv_folds(indi_cv_marg_folds)
                        indi_cv_marg_fit_lst <- lapply(as.list(1:length(unique(indi_cv_marg_folds_vec))), function(x) indi_cv_marg_fit[indi_cv_marg_folds_vec == x])
                        indi_marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = indi_cv_marg_folds_vec, inner_folds_2 = geog_glm_cv_folds_vec))
                        ## get individual, cv vimp
                        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_var_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = indi_cv_marg_fit_lst, f2 = geog_glm_cv_fit_lst, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_marg_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
                    }
                }
                ## merge together
                eval(parse(text = paste0(this_outcome_name, "_vimp_lst$ind_marginal <- merge_vim(", paste(paste0(this_outcome_name, "_marg_", var_inds), collapse = ", "), ")")))
                if (opts$cvperf) {
                    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$ind_marginal <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", var_inds), collapse = ", "), ")")))
                }
            }
            ## save them off
            eval(parse(text = paste0("saveRDS(", this_outcome_name, "_vimp_lst, file = '/home/slfits/", paste0(this_outcome_name, "_vimp"), ".rds')")))
            eval(parse(text = paste0("saveRDS(", this_outcome_name, "_cv_vimp_lst, file = '/home/slfits/", paste0(this_outcome_name, "_cv_vimp"), ".rds')")))
        }
    }
}
