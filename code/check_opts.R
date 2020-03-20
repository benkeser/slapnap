#! /usr/bin/env Rscript
library("shiny")
source("/home/lib/utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# --------------------------------------
# make sure that the options are in the form we expected
# --------------------------------------
# check the nab
shiny::validate(
    shiny::need(opts$nab, "Please enter at least one antibody with data in the CATNAP database. Multiple antibodies may be specified in a semicolon-separated list (e.g., 'VRC07-523-LS;PGT121')")
)
# check the outcome
all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
shiny::validate(
    shiny::need(opts$outcomes, "Please enter at least one outcome (one of 'ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
)
shiny::validate(
    shiny::need(opts$outcomes != '', "Please enter at least one outcome (one of 'ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80')."),
    shiny::need(length(setdiff(opts$outcomes, all_outcomes)) == 0, "You have entered one or more outcomes that are not supported at this time. Please enter at least one of the currently supported outcomes ('ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
)
# check the learners
all_learners <- c("rf", "xgboost", "lasso")
shiny::validate(
    shiny::need(opts$learners, "Please enter one or more learners ('rf', 'xgboost', 'lasso') in a semicolon-separated list (e.g., 'rf;lasso')."),
    shiny::need(length(setdiff(opts$learners, all_learners)) == 0, "You have entered one or more learners that are not supported at this time. Please enter at least one of the currently supported learners ('rf', 'lasso', 'xgboost') or a semicolon-separated list of learners (e.g., 'rf;lasso').")
)
# check cv tune
shiny::validate(
    shiny::need(opts$cvtune == TRUE | !opts$cvtune == TRUE, "cvtune must be logical (i.e., either 'TRUE' or 'FALSE').")
)
# check cv perf
shiny::validate(
    shiny::need(opts$cvperf == TRUE | !opts$cvperf == TRUE, "cvperf must be logical (i.e., either 'TRUE' or 'FALSE').")
)
# check importance
all_importance_grp <- c("marg", "cond")
all_importance_ind <- c("marg", "cond")
shiny::validate(
    shiny::need(length(setdiff(opts$importance_grp, all_importance_grp)) == 0 | opts$importance_grp == "", "Please enter either a semicolon-separated list of supported group importance types (i.e., 'marg', 'cond', or 'marg;cond') or an empty string ('')."),
    shiny::need(length(setdiff(opts$importance_ind, all_importance_ind)) == 0 | opts$importance_ind == "", "Please enter either a semicolon-separated list of supported individual importance types (i.e., 'marg', 'cond', or 'marg;cond') or an empty string ('').")
)
# check report name
shiny::validate(
    shiny::need(grepl(".html", opts$report_name) | opts$report_name == "", "Please either enter a desired name for your report (with the file extension .html) or a blank string ('').")
)
# check return object calls
shiny::validate(
    shiny::need(opts$return_full_sl_obj == TRUE | !opts$return_full_sl_obj == TRUE, "Please enter a logical value for whether or not the full regression object should be returned."),
    shiny::need(opts$return_analysis_dataset == TRUE | !opts$return_analysis_dataset == TRUE, "Please enter a logical value for whether or not the full regression object should be returned.")
)
