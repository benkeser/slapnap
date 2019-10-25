#! /usr/bin/env Rscript

## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("vimp")
library("dplyr")
source("/home/lib/variable_groups.R")

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"
no_cv <- Sys.getenv("no_cv") == "TRUE"

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

# load all superlearner fits
sl_fit_names <- list.files("/home/slfits")
sl_fitted_names <- paste0("/home/slfits/", sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("cvfit_", sl_fit_names) & !grepl("slweights", sl_fit_names) & !grepl("report", sl_fit_names)])
all_sl_fits <- lapply(as.list(sl_fitted_names), readRDS)
names(all_sl_fits) <- sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("cvfit_", sl_fit_names) & !grepl("slweights", sl_fit_names) & !grepl("report", sl_fit_names)]
all_sl_nms <- names(all_sl_fits)
## subset to full, reduced
full_sl_fits <- all_sl_fits[grepl("fitted_", all_sl_nms) & !grepl("cvfitted_", all_sl_nms) & !grepl("minus", all_sl_nms) & !grepl("slweights", all_sl_nms)]
full_cv_sl_fits <- all_sl_fits[grepl("cvfitted_", all_sl_nms) & !grepl("minus", all_sl_nms) & !grepl("slweights", all_sl_nms)]
reduced_sl_fits <- all_sl_fits[grepl("fitted_", all_sl_nms) & !grepl("cvfitted_", all_sl_nms) & grepl("minus", all_sl_nms) & !grepl("slweights", all_sl_nms)]
reduced_cv_sl_fits <- all_sl_fits[grepl("cvfitted_", all_sl_nms) & grepl("minus", all_sl_nms) & !grepl("slweights", all_sl_nms)]
## get folds for cv
full_cv_sl_folds <- all_sl_fits[grepl("cvfolds_", all_sl_nms) & !grepl("minus", all_sl_nms)]


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
all_outcome_var_lst <- unique(gsub("__", "_", gsub(".rds", "", gsub("cv", "", gsub("fitted_", "", sl_fit_names[!grepl("fit_", sl_fit_names) & !grepl("folds", sl_fit_names) & !grepl("slweights", sl_fit_names) & !grepl(".RData", sl_fit_names)])))))
all_outcome_lst <- all_outcome_var_lst[!grepl("minus", all_outcome_var_lst)]
all_var_grp_lst <- all_outcome_var_lst[grepl("minus", all_outcome_var_lst)]

## ---------------------------------------------------------------------------
## get variable importance!
## ---------------------------------------------------------------------------

# for continuous outcomes, do r-squared
continuous_outcomes <- all_outcome_lst[grepl("ic50", all_outcome_lst) | grepl("ic80", all_outcome_lst) | grepl("iip", all_outcome_lst)]
continuous_grps <- unlist(lapply(strsplit(all_var_grp_lst[grepl("ic50", all_var_grp_lst) | grepl("ic80", all_var_grp_lst) | grepl("iip", all_var_grp_lst)], "_", fixed = TRUE), function(x) tail(x, n = 1)))
continuous_outcome_vimp <- vector("list", length = length(continuous_outcomes))
continuous_outcome_cv_vimp <- vector("list", length = length(continuous_outcomes))
continuous_nms_grid <- expand.grid(outcome = continuous_outcomes, grp = continuous_grps)
set.seed(474747)
for (i in 1:length(continuous_outcome_vimp)) {
    ## make sub-folds for non-cv
    sub_folds <- sample(1:2, length(dat[, continuous_outcomes[i]]), replace = TRUE, prob = c(0.5, 0.5))

    for (j in 1:length(continuous_reduced_sl_fits)) {
        eval(parse(text = paste0(paste0(unique(continuous_nms_grid$outcome)[i], "_", unique(continuous_nms_grid$grp)[j]),
                                 " <- vimp::vim(Y = dat[, continuous_outcomes[i]],
                                        f1 = continuous_sl_fits[[i]],
                                        f2 = continuous_reduced_sl_fits[[j]],
                                        indx = which(pred_names %in% unlist(all_var_groups[grepl(continuous_grps[j], names(all_var_groups))])),
                                        run_regression = FALSE,
                                        alpha = 0.05,
                                        type = 'r_squared',
                                        folds = sub_folds,
                                        na.rm = TRUE,
                                        scale = 'identity')")))
    }
    eval(parse(text = paste0("continuous_outcome_vimp[[i]] <- merge_vim(", paste(paste0(continuous_nms_grid$outcome[i], "_", unique(continuous_nms_grid$grp)), collapse = ", "), ")")))
    if (!no_cv) {
        ## make a vector of folds
        V <- length(continuous_cv_sl_folds[[i]])
        v_lst <- sapply(1:V, function(s) rep(s, length(continuous_cv_sl_folds[[i]][[s]])))
        joint_lst <- mapply(list, v_lst, continuous_cv_sl_folds[[i]], SIMPLIFY = FALSE)
        folds_mat <- do.call(rbind, lapply(joint_lst, function(x) cbind(x[[1]], x[[2]])))
        folds <- folds_mat[order(folds_mat[, 2]), 1]
        ## make a list of the outcome, fitted values
        full_lst <- lapply(as.list(1:length(unique(folds))), function(x) continuous_cv_sl_fits[[i]][folds == x])

        for (j in 1:length(continuous_reduced_cv_sl_fits)) {
            redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) continuous_reduced_cv_sl_fits[[j]][folds == x])

            eval(parse(text = paste0(paste0("cv_", unique(continuous_nms_grid$outcome)[i], "_", unique(continuous_nms_grid$grp)[j]),
                                     " <- vimp::cv_vim(Y = dat[, continuous_outcomes[i]],
                                            f1 = full_lst,
                                            f2 = redu_lst,
                                            indx = which(pred_names %in% unlist(all_var_groups[grepl(continuous_grps[j], names(all_var_groups))])),
                                            run_regression = FALSE,
                                            alpha = 0.05,
                                            type = 'r_squared',
                                            folds = folds,
                                            na.rm = TRUE,
                                            scale = 'identity')")))
        }
        eval(parse(text = paste0("continuous_outcome_cv_vimp[[i]] <- merge_vim(", paste(paste0("cv_", continuous_nms_grid$outcome[i], "_", unique(continuous_nms_grid$grp)), collapse = ", "), ")")))
    }
}
## save them off
saveRDS(continuous_outcome_vimp, "/home/slfits/continuous_outcome_vimp.rds")
if (!no_cv) {
    saveRDS(continuous_outcome_cv_vimp, "/home/slfits/continuous_outcome_cv_vimp.rds")
}

# for binary outcomes, do AUC (this only, for now)
binary_outcomes <- all_outcome_lst[grepl("dichotomous.1", all_outcome_lst) | grepl("dichotomous", all_outcome_lst)]
binary_grps <- unlist(lapply(strsplit(all_var_grp_lst[grepl("dichotomous", all_var_grp_lst)], "_", fixed = TRUE), function(x) tail(x, n = 1)))
binary_outcome_vimp <- vector("list", length = length(binary_outcomes))
binary_outcome_cv_vimp <- vector("list", length = length(binary_outcomes))
binary_nms_grid <- expand.grid(outcome = binary_outcomes, grp = binary_grps)
set.seed(474747)
if (!reduce_outcomes) {
    for (i in 1:length(binary_outcome_vimp)) {
        ## make sub-folds for non-cv
        sub_folds <- sample(1:2, length(dat[, binary_outcomes[i]]), replace = TRUE, prob = c(0.5, 0.5))

        for (j in 1:length(binary_reduced_sl_fits)) {
            eval(parse(text = paste0(paste0(unique(binary_nms_grid$outcome)[i], "_", unique(binary_nms_grid$grp)[j]),
                                     " <- vimp::vim(Y = dat[, binary_outcomes[i]],
                                            f1 = binary_sl_fits[[i]],
                                            f2 = binary_reduced_sl_fits[[j]],
                                            indx = which(pred_names %in% unlist(all_var_groups[grepl(binary_grps[j], names(all_var_groups))])),
                                            run_regression = FALSE,
                                            alpha = 0.05,
                                            type = 'auc',
                                            folds = sub_folds,
                                            na.rm = TRUE,
                                            scale = 'identity')")))
        }
        eval(parse(text = paste0("binary_outcome_vimp[[i]] <- merge_vim(", paste(paste0(binary_nms_grid$outcome[i], "_", unique(binary_nms_grid$grp)), collapse = ", "), ")")))
        if (!no_cv) {
            ## make a vector of folds
            V <- length(binary_cv_sl_folds[[i]])
            v_lst <- sapply(1:V, function(s) rep(s, length(binary_cv_sl_folds[[i]][[s]])))
            joint_lst <- mapply(list, v_lst, binary_cv_sl_folds[[i]], SIMPLIFY = FALSE)
            folds_mat <- do.call(rbind, lapply(joint_lst, function(x) cbind(x[[1]], x[[2]])))
            folds <- folds_mat[order(folds_mat[, 2]), 1]
            ## make a list of the outcome, fitted values
            full_lst <- lapply(as.list(1:length(unique(folds))), function(x) binary_cv_sl_fits[[i]][folds == x])

            for (j in 1:length(binary_reduced_cv_sl_fits)) {
                redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) binary_reduced_cv_sl_fits[[j]][folds == x])

                eval(parse(text = paste0(paste0("cv_", unique(binary_nms_grid$outcome)[i], "_", unique(binary_nms_grid$grp)[j]),
                                         " <- vimp::cv_vim(Y = dat[, binary_outcomes[i]],
                                                f1 = full_lst,
                                                f2 = redu_lst,
                                                indx = which(pred_names %in% unlist(all_var_groups[grepl(binary_grps[j], names(all_var_groups))])),
                                                run_regression = FALSE,
                                                alpha = 0.05,
                                                type = 'auc',
                                                folds = folds,
                                                na.rm = TRUE,
                                                scale = 'identity')")))
            }
            eval(parse(text = paste0("binary_outcome_cv_vimp[[i]] <- merge_vim(", paste(paste0("cv_", binary_nms_grid$outcome[i], "_", unique(binary_nms_grid$grp)), collapse = ", "), ")")))
        }
    }
}
## save them off
if (!reduce_outcomes) {
    saveRDS(binary_outcome_vimp, "/home/slfits/binary_outcome_vimp.rds")
    if (!no_cv) {
        saveRDS(binary_outcome_vimp, "/home/slfits/binary_outcome_cv_vimp.rds")
    }
}
