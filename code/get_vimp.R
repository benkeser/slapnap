#! /usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("vimp")
library("dplyr")
source("/home/lib/variable_groups.R")
source("/home/lib/utils.R")

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"
no_cv <- Sys.getenv("no_cv") == "TRUE"

## option to run individual-level vimp; defaults to FALSE
run_indi_vimp <- Sys.getenv("run_indi_vimp") == "TRUE"

# load data
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
dat <- dat[complete.cases(dat),]

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
# if reduce_covs, only do individual imp on that number
if (reduce_covs) {
    num_covs <- 10
    var_inds <- pred_names[!grepl("geog", pred_names)][1:(num_covs - length(all_geog_vars))]
} else {
    num_covs <- length(pred_names) - length(all_geog_vars)
    var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]
}

## ---------------------------------------------------------------------------
## get variable importance!
## ---------------------------------------------------------------------------
set.seed(474747)
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    vimp_opts <- get_vimp_options(this_outcome_name)
    ## create output list
    eval(parse(text = paste0(this_outcome_name, '_vimp_lst <- make_vimp_list(var_groups, var_inds)')))
    eval(parse(text = paste0(this_outcome_name, '_cv_vimp_lst <- make_vimp_list(var_groups, var_inds)')))
    ## load full fit corresponding to this outcome
    full_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, ".rds"))
    full_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, ".rds"))
    full_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, ".rds"))
    full_cv_folds_vec <- get_cv_folds(full_cv_folds)
    full_cv_fit_lst <- lapply(as.list(1:length(unique(full_cv_folds_vec))), function(x) full_cv_fit[full_cv_folds_vec == x])
    ## load geog-only fit corresponding to this outcome
    geog_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_geog.rds"))
    geog_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_geog.rds"))
    geog_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_geog.rds"))
    geog_cv_folds_vec <- get_cv_folds(geog_cv_folds)
    geog_cv_fit_lst <- lapply(as.list(1:length(unique(geog_cv_folds_vec))), function(x) geog_cv_fit[geog_cv_folds_vec == x])
    ## load in outer folds for VIM
    outer_folds <- readRDS(paste0("/home/slfits/", this_outcome_name, "_outer_folds.rds"))
    ## group variable importance
    for (j in 1:length(var_groups)) {
        this_group_name <- names(var_groups)[j]
        ## non-cv
        cond_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_minus_", this_group_name, ".rds"))
        marg_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
        ## cv
        cond_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_minus_", this_group_name, ".rds"))
        marg_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
        cond_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_minus_", this_group_name, ".rds"))
        marg_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
        cond_cv_folds_vec <- get_cv_folds(cond_cv_folds)
        marg_cv_folds_vec <- get_cv_folds(marg_cv_folds)
        cond_cv_fit_lst <- lapply(as.list(1:length(unique(cond_cv_folds_vec))), function(x) cond_cv_fit[cond_cv_folds_vec == x])
        marg_cv_fit_lst <- lapply(as.list(1:length(unique(marg_cv_folds_vec))), function(x) marg_cv_fit[marg_cv_folds_vec == x])
        cond_folds <- list(outer_folds = outer_folds, inner_folds = list(inner_folds_1 = full_cv_folds_vec, inner_folds_2 = cond_cv_folds_vec))
        marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = marg_cv_folds_vec, inner_folds_2 = geog_cv_folds_vec))
        ## get conditional, non-cv vimp
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cond_", this_group_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = full_fit, f2 = cond_fit, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = cond_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
        ## get marginal, non-cv vimp
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_group_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = marg_fit, f2 = geog_fit, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = marg_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
        ## get conditional, cv vimp
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_cond_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = full_cv_fit_lst, f2 = cond_cv_fit_lst, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = cond_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
        ## get marginal, cv vimp
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = marg_cv_fit_lst, f2 = geog_cv_fit_lst, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = marg_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
    }
    ## merge together
    eval(parse(text = paste0(this_outcome_name, "_vimp_lst$conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cond_", names(var_groups)), collapse = ", "), ")")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cv_cond_", names(var_groups)), collapse = ", "), ")")))
    eval(parse(text = paste0(this_outcome_name, "_vimp_lst$marginal <- merge_vim(", paste(paste0(this_outcome_name, "_marg_", names(var_groups)), collapse = ", "), ")")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$marginal <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", names(var_groups)), collapse = ", "), ")")))
    ## individual variable importance
    if (run_indi_vimp) {
        for (j in 1:length(var_inds)) {
            this_var_name <- var_inds[j]
            ## non-cv
            indi_fit <- readRDS(paste0("/home/slfits/fitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
            ## cv
            indi_cv_fit <- readRDS(paste0("/home/slfits/cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
            indi_cv_folds <- readRDS(paste0("/home/slfits/cvfolds_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
            indi_cv_folds_vec <- get_cv_folds(indi_cv_folds)
            indi_cv_fit_lst <- lapply(as.list(1:length(unique(indi_cv_folds_vec))), function(x) indi_cv_fit[indi_cv_folds_vec == x])
            indi_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = indi_cv_folds_vec, inner_folds_2 = geog_cv_folds_vec))
            ## get individual, non-cv vimp
            suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_marg_", this_var_name, " <- vimp::vim(Y = dat[, this_outcome_name], f1 = indi_fit, f2 = geog_fit, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_folds$outer_folds, na.rm = TRUE, scale = 'identity')"))))
            ## get individual, cv vimp
            suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_var_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = indi_cv_fit_lst, f2 = geog_cv_fit_lst, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = vimp_opts$vimp_measure, folds = indi_folds, V = V, na.rm = TRUE, scale = 'identity')"))))
        }
        ## merge together
        eval(parse(text = paste0(this_outcome_name, "_vimp_lst$individual <- merge_vim(", paste(paste0(this_outcome_name, "_marg_", var_inds), collapse = ", "), ")")))
        eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$individual <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", var_inds), collapse = ", "), ")")))
    }
    ## save them off
    eval(parse(text = paste0("saveRDS(", this_outcome_name, "_vimp_lst, file = '/home/slfits/", paste0(this_outcome_name, "_vimp"), ".rds')")))
    eval(parse(text = paste0("saveRDS(", this_outcome_name, "_cv_vimp_lst, file = '/home/slfits/", paste0(this_outcome_name, "_cv_vimp"), ".rds')")))
}
