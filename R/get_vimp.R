#! /usr/bin/env Rscript

## get variable importance!

# load libraries
library("SuperLearner")
library("vimp")
source("/home/lib/variable_groups.R")

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"

# load data
analysis_data_name <- list.files("/home/dat/analysis")
dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
dat <- dat[complete.cases(dat),]

# get pre-defined groups
# get names of predictors
non_pred_names <- c("pc.ic50", "pc.ic80", "iip",
                    "dichotomous.1", "dichotomous.2",
                    "seq.id.lanl","seq.id.catnap")
pred_names <- colnames(dat)[!(colnames(dat) %in% non_pred_names)]
all_var_groups <- get_variable_groups(dat, pred_names)

# load all superlearner fits
sl_fit_names <- list.files("/home/slfits")
sl_fitted_names <- paste0("/home/slfits/", sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("cvfit_", sl_fit_names)])
all_sl_fits <- lapply(as.list(sl_fitted_names), readRDS)
names(all_sl_fits) <- sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("cvfit_", sl_fit_names)]
full_sl_fits <- all_sl_fits[grepl("fitted_", sl_fitted_names) & !grepl("cvfitted_", sl_fitted_names) & !grepl("minus", sl_fitted_names)]
full_cv_sl_fits <- all_sl_fits[grepl("cvfitted_", sl_fitted_names) & !grepl("minus", sl_fitted_names)]
reduced_sl_fits <- all_sl_fits[grepl("fitted_", sl_fitted_names) & !grepl("cvfitted_", sl_fitted_names) & grepl("minus", sl_fitted_names)]
reduced_cv_sl_fits <- all_sl_fits[grepl("cvfitted_", sl_fitted_names) & grepl("minus", sl_fitted_names)]
full_cv_sl_folds <- all_sl_fits[grepl("cvfolds_", sl_fit_names) & !grepl("minus", sl_fit_names)]

## split out into continuous, dichotomous
continuous_sl_fits <- full_sl_fits[grepl("ic50", names(full_sl_fits)) | grepl("ic80", names(full_sl_fits)) | grepl("iip", names(full_sl_fits))]
continuous_reduced_sl_fits <- reduced_sl_fits[grepl("ic50", names(reduced_sl_fits)) | grepl("ic80", names(reduced_sl_fits)) | grepl("iip", names(reduced_sl_fits))]
continuous_cv_sl_fits <- full_cv_sl_fits[grepl("ic50", names(full_cv_sl_fits)) | grepl("ic80", names(full_cv_sl_fits)) | grepl("iip", names(full_cv_sl_fits))]
continuous_reduced_cv_sl_fits <- reduced_cv_sl_fits[grepl("ic50", names(reduced_cv_sl_fits)) | grepl("ic80", names(reduced_cv_sl_fits)) | grepl("iip", names(reduced_cv_sl_fits))]
continuous_cv_sl_folds <- full_cv_sl_folds[grepl("ic50", names(full_cv_sl_folds)) | grepl("ic80", names(full_cv_sl_folds)) | grepl("iip", names(full_cv_sl_folds))]
binary_sl_fits <- full_sl_fits[grepl("dichotomous.1", names(full_sl_fits)) | grepl("dichotomous.2", names(full_sl_fits))]
binary_reduced_sl_fits <- reduced_sl_fits[grepl("dichotomous.1", names(reduced_sl_fits)) | grepl("dichotomous.2", names(reduced_sl_fits))]
binary_cv_sl_fits <- full_cv_sl_fits[grepl("dichotomous.1", names(full_cv_sl_fits)) | grepl("dichotomous.2", names(full_cv_sl_fits))]
binary_reduced_cv_sl_fits <- reduced_cv_sl_fits[grepl("dichotomous.1", names(reduced_cv_sl_fits)) | grepl("dichotomous.2", names(reduced_cv_sl_fits))]
binary_cv_sl_folds <- full_cv_sl_folds[grepl("dichotomous.1", names(full_cv_sl_folds)) | grepl("dichotomous.2", names(full_cv_sl_folds))]

## get the outcomes in the list
all_outcome_var_lst <- unique(gsub(".rds", "", gsub("cv", "", gsub("fitted_", "", sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("folds", sl_fit_names) & !grepl(".RData", sl_fit_names)]))))
all_outcome_lst <- all_outcome_var_lst[!grepl("minus", all_outcome_var_lst)]
all_var_grp_lst <- all_outcome_var_lst[grepl("minus", all_outcome_var_lst)]

# for continuous outcomes, do r-squared
continuous_outcomes <- all_outcome_lst[grepl("ic50", all_outcome_lst) | grepl("ic80", all_outcome_lst) | grepl("iip", all_outcome_lst)]
continuous_grps <- tail(unlist(strsplit(all_var_grp_lst[grepl("ic50", all_outcome_lst) | grepl("ic80", all_outcome_lst) | grepl("iip", all_outcome_lst)], "_", fixed = TRUE)), n = 1)
continuous_outcome_vimp <- vector("list", length = length(continuous_outcomes))
continuous_outcome_cv_vimp <- vector("list", length = length(continuous_outcomes))

set.seed(474747)
for (i in 1:length(continuous_outcome_vimp)) {
    ## make sub-folds for non-cv
    sub_folds <- sample(1:2, length(dat[, continuous_outcomes[i]]), replace = TRUE, prob = c(0.5, 0.5))

    continuous_outcome_vimp[[i]] <- vector("list", length = length(continuous_grps))
    continuous_outcome_cv_vimp[[i]] <- vector("list", length = length(continuous_grps))
    for (j in 1:length(continuous_reduced_sl_fits)) {
        continuous_outcome_vimp[[i]][[j]] <- vimp::vim(Y = dat[, continuous_outcomes[i]], 
                                                            f1 = continuous_sl_fits[[i]],
                                                            f2 = continuous_reduced_sl_fits[[j]],
                                                            indx = which(pred_names %in% unlist(all_var_groups[grepl(continuous_grps[j], names(all_var_groups))])),
                                                            run_regression = FALSE,
                                                            alpha = 0.05,
                                                            type = "r_squared",
                                                            folds = sub_folds)
    }
    continuous_outcome_cv_vimp[[i]] <- vector("list", length = length(continuous_reduced_cv_sl_fits))
    ## make a vector of folds
    max_in_fold <- max(unlist(lapply(continuous_cv_sl_folds[[i]], function(x) length(x))))
    folds_tmp <- lapply(continuous_cv_sl_folds[[i]], function(x) {
        if (length(x) < max_in_fold) {
            c(x, rep(NA, max_in_fold - length(x)))
        } else {
            x
        }
    })
    fold_row_nums <- as.vector(do.call(cbind, folds_tmp))
    folds_init <- rep(as.numeric(names(folds_tmp)), 
                      each = round(length(dat[, continuous_outcomes[i]])/length(folds_tmp)))
    folds_mat <- cbind(fold_row_nums, folds_init)
    folds <- folds_mat[order(folds_mat[, 1]), 2][!is.na(folds_mat[, 1])]
    ## make a list of the outcome, fitted values
    full_lst <- lapply(as.list(1:length(unique(folds))), function(x) continuous_cv_sl_fits[[i]][folds == x])
    for (j in 1:length(continuous_reduced_cv_sl_fits)) {
        redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) continuous_reduced_cv_sl_fits[[j]][folds == x])
    
        continuous_outcome_cv_vimp[[i]][[j]] <- vimp::cv_vim(Y = dat[, continuous_outcomes[i]],
                                                             f1 = full_lst,
                                                             f2 = redu_lst,
                                                             indx = which(pred_names %in% unlist(all_var_groups[grepl(continuous_grps[j], names(all_var_groups))])),
                                                             run_regression = FALSE,
                                                             alpha = 0.05,
                                                             type = "r_squared",
                                                             folds = folds)
    }
}

# for binary outcomes, do AUC (this only, for now)
binary_outcomes <- all_outcome_lst[grepl("dichotomous.1", all_outcome_lst) | grepl("dichotomous", all_outcome_lst) | grepl("iip", all_outcome_lst)]
binary_grps <- all_var_grp_lst[grepl("ic50", all_outcome_lst) | grepl("ic80", all_outcome_lst) | grepl("iip", all_outcome_lst)]
binary_outcome_vimp <- vector("list", length = length(binary_outcomes))
binary_outcome_cv_vimp <- vector("list", length = length(binary_outcomes))

set.seed(474747)
for (i in 1:length(binary_outcome_vimp)) {
    ## make sub-folds for non-cv
    sub_folds <- sample(1:2, length(dat[, binary_outcomes[i]]), replace = TRUE, prob = c(0.5, 0.5))

    binary_outcome_vimp[[i]] <- vector("list", length = length(binary_grps))
    binary_outcome_cv_vimp[[i]] <- vector("list", length = length(binary_grps))
    for (j in 1:length(binary_reduced_sl_fits)) {
        binary_outcome_vimp[[i]][[j]] <- vimp::vim(Y = dat[, binary_outcomes[i]], 
                                                            f1 = binary_sl_fits[[i]],
                                                            f2 = binary_reduced_sl_fits[[j]],
                                                            indx = which(pred_names %in% unlist(all_var_groups[grepl(binary_grps[j], names(all_var_groups))])),
                                                            run_regression = FALSE,
                                                            alpha = 0.05,
                                                            type = "auc",
                                                            folds = sub_folds)
    }
    binary_outcome_cv_vimp[[i]] <- vector("list", length = length(reduced_cv_sl_fits))
    max_in_fold <- max(unlist(lapply(binary_cv_sl_folds[[i]], function(x) length(x))))
    folds_tmp <- lapply(binary_cv_sl_folds[[i]], function(x) {
        if (length(x) < max_in_fold) {
            c(x, rep(NA, max_in_fold - length(x)))
        } else {
            x
        }
    })
    fold_row_nums <- as.vector(do.call(cbind, folds_tmp))
    folds_init <- rep(as.numeric(names(folds_tmp)), 
                      each = round(length(dat[, binary_outcomes[i]])/length(folds_tmp)))
    folds_mat <- cbind(fold_row_nums, folds_init)
    folds <- folds_mat[order(folds_mat[, 1]), 2][!is.na(folds_mat[, 1])]
    for (j in 1:length(binary_reduced_cv_sl_fits)) {
        
    
        binary_outcome_cv_vimp[[i]][[j]] <- vimp::cv_vim(Y = dat[, binary_outcomes[i]],
                                                             f1 = binary_cv_sl_fits[[i]],
                                                             f2 = binary_reduced_cv_sl_fits[[i]],
                                                             indx = which(pred_names %in% unlist(all_var_groups[grepl(binary_grps[j], names(all_var_groups))])),
                                                             run_regression = FALSE,
                                                             alpha = 0.05,
                                                             type = "auc",
                                                             folds = folds)
    }
}

## save them off
saveRDS(continuous_outcome_vimp, "/home/slfits/continuous_outcome_vimp.rds")
saveRDS(continuous_outcome_cv_vimp, "/home/slfits/continuous_outcome_cv_vimp.rds")

saveRDS(binary_outcome_vimp, "/home/slfits/binary_outcome_vimp.rds")
saveRDS(binary_outcome_vimp, "/home/slfits/binary_outcome_cv_vimp.rds")