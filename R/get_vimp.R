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

# for continuous outcomes, do r-squared

# for binary outcomes, do AUC, accuracy, deviance