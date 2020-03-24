#! /usr/bin/env Rscript
library("shiny")
source("/home/lib/utils.R")
source("/home/lib/check_opts_functions.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# --------------------------------------
# make sure that the options are in the form we expected
# --------------------------------------
# check the nab
check_opts_nab(opts$nab)
# check the outcome
all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
check_opts_outcomes(opts$outcomes, all_outcomes)
# check the learners
all_learners <- c("rf", "xgboost", "lasso")
check_opts_learners(opts$learners, all_learners)
# check cv tune
check_opts_cvtune(opts$cvtune)
# check cv perf
check_opts_cvperf(opts$cvperf)
# check importance
all_importance_grp <- c("marg", "cond")
all_importance_ind <- c("marg", "cond")
check_opts_vimp(opts$importance_grp, opts$importance_ind, all_importance_grp, all_importance_ind)
# check report name
check_opts_report_name(opts$report_name)
# check return object calls
check_opts_returns(opts$return_full_sl_obj, opts$return_analysis_dataset)
