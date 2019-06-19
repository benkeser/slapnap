#! /usr/bin/env Rscript

## get variable importance!

# load libraries
library("SuperLearner")
library("vimp")
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

# get pre-defined groups
all_var_groups <- get_variable_groups(dat, pred_names)

# load all superlearner fits
sl_fit_names <- list.files("/home/slfits")
all_sl_fits <- lapply(as.list(sl_fit_names), readRDS)
full_sl_fits <- all_sl_fits[grepl("fit_", sl_fit_names)]
full_cv_sl_fits <- all_sl_fits[grepl("cvfit_", sl_fit_names)]
reduced_sl_fits <- all_sl_fits[grepl("fitted_", sl_fit_names)]
reduced_cv_sl_fits <- all_sl_fits[grepl("cvfitted_", sl_fit_names)]
full_cv_sl_folds <- all_sl_fits[grepl("cvfolds_", sl_fit_names)]

continuous_outcomes <- c("pc.ic50", "pc.ic80", "iip")
binary_outcomes <- c("dichotomous.1", "dichotomous.2")


# for continuous outcomes, do r-squared
continuous_outcome_vimp <- vector("list", length = length(continuous_outcomes))
continuous_outcome_cv_vimp <- vector("list", length = length(continuous_outcomes))
set.seed(474747)
for (i in 1:length(continuous_outcome_vimp)) {
    ## make sub-folds for non-cv
    sub_folds <- sample(1:2, length(dat[, continuous_outcomes[i]]), replace = TRUE, prob = c(0.5, 0.5))

    continuous_outcome_vimp[[i]] <- vector("list", length = length(reduced_sl_fits))
    for (j in 1:length(reduced_sl_fits)) {
        continuous_outcome_vimp[[i]][[j]] <- vimp::vimp_rsquared(Y = dat[, continuous_outcomes[i]], 
                                                            f1 = full_sl_fits[[i]],
                                                            f2 = reduced_sl_fits[[j]],
                                                            indx = ,
                                                            run_regression = FALSE,
                                                            alpha = 0.05,
                                                            folds = sub_folds)
    }
    continuous_outcome_cv_vimp[[i]] <- vector("list", length = length(reduced_cv_sl_fits))
    for (j in 1:length(reduced_cv_sl_fits)) {
        continuous_outcome_cv_vimp[[i]][[j]] <- vimp::cv_vim(Y = dat[, continuous_outcomes[i]],
                                                             f1 = full_cv_sl_fits[[i]],
                                                             f2 = reduced_cv_sl_fits[[i]],
                                                             indx = ,
                                                             run_regression = FALSE,
                                                             alpha = 0.05,
                                                             folds = full_cv_sl_folds[[i]])
    }
}

# for binary outcomes, do AUC, accuracy, deviance
binary_outcome_vimp <- vector("list", length = length(binary_outcomes))