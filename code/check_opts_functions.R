#! /usr/bin/env Rscript
library("shiny")
source("/home/lib/utils.R")

# check nab
check_opts_nab <- function(nab_str) {
    shiny::validate(
        shiny::need(nab_str, "Please enter at least one antibody with data in the CATNAP database. Multiple antibodies may be specified in a semicolon-separated list (e.g., 'VRC07-523-LS;PGT121')")
    )
}
# check outcomes
check_opts_outcomes <- function(outcome_vec, all_outcomes) {
    shiny::validate(
        shiny::need(outcome_vec, "Please enter at least one outcome (one of 'ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
    )
    shiny::validate(
        shiny::need(outcome_vec != '', "Please enter at least one outcome (one of 'ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80')."),
        shiny::need(length(setdiff(outcome_vec, all_outcomes)) == 0, "You have entered one or more outcomes that are not supported at this time. Please enter at least one of the currently supported outcomes ('ic50', 'ic80', 'iip', 'sens1', 'sens2') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
    )
}
# check learners
check_opts_learners <- function(learner_vec, all_learners) {
    shiny::validate(
        shiny::need(learner_vec, "Please enter one or more learners ('rf', 'xgboost', 'lasso') in a semicolon-separated list (e.g., 'rf;lasso')."),
        shiny::need(length(setdiff(learner_vec, all_learners)) == 0, "You have entered one or more learners that are not supported at this time. Please enter at least one of the currently supported learners ('rf', 'lasso', 'xgboost') or a semicolon-separated list of learners (e.g., 'rf;lasso').")
    )
}
# check cvtune
check_opts_cvtune <- function(cvtune_str) {
    shiny::validate(
        shiny::need(cvtune_str == TRUE | !opts$cvtune == TRUE, "cvtune must be logical (i.e., either 'TRUE' or 'FALSE').")
    )

}
# check cvperf
check_opts_cvperf <- function(cvperf_str) {
    shiny::validate(
        shiny::need(cvperf_str == TRUE | !opts$cvperf == TRUE, "cvperf must be logical (i.e., either 'TRUE' or 'FALSE').")
    )
}
# check vimp
check_opts_vimp <- function(grp_str, ind_str, all_grp, all_ind) {
    shiny::validate(
        shiny::need(length(setdiff(grp_str, all_grp)) == 0 | grp_str == "", "Please enter either a semicolon-separated list of supported group importance types (i.e., 'marg', 'cond', or 'marg;cond') or an empty string ('')."),
        shiny::need(length(setdiff(ind_str, all_ind)) == 0 | ind_str == "", "Please enter either a semicolon-separated list of supported individual importance types (i.e., 'marg', 'cond', or 'marg;cond') or an empty string ('').")
    )
}
# check report name
check_opts_report_name <- function(report_name) {
    shiny::validate(
        shiny::need(grepl(".html", report_name, fixed = TRUE) | report_name == "", "Please either enter a desired name for your report (with the file extension .html) or a blank string ('').")
    )
}
# check object returns
check_opts_returns <- function(sl_obj, dataset) {
    shiny::validate(
        shiny::need(sl_obj == TRUE | !sl_obj == TRUE, "Please enter a logical value for whether or not the full regression object should be returned."),
        shiny::need(dataset == TRUE | !dataset == TRUE, "Please enter a logical value for whether or not the full regression object should be returned.")
    )
}
