#! /usr/bin/env Rscript

# check nab
check_opts_nab <- function(nab_str) {
    shiny::validate(
        shiny::need(nab_str, "Please enter at least one bNAb with data in the CATNAP database. Multiple bNAbs may be specified in a semicolon-separated list (e.g., 'VRC07-523-LS;PGT121')")
    )
}
# check outcomes
check_opts_outcomes <- function(outcome_vec, all_outcomes, n_abs) {
    shiny::validate(
        shiny::need(outcome_vec, "Please enter at least one outcome (one of 'ic50', 'ic80', 'iip', 'estsens', 'multsens', 'sens') or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
    )
    if(n_abs == 1){
        shiny::validate(
            shiny::need(!any(c("sens2") %in% outcome_vec), "If a single bNAb is used, please specify 'sens' as the outcome, instead of 'estsens' or 'multsens'.")
        )
    } else {
        shiny::validate(
            shiny::need(!("sens" %in% Sys.getenv("outcomes")), "If multiple bNAbs are used, please specify 'estsens' and/or 'multsens' as the outcome, instead of 'sens'.")
        )
    }
    shiny::validate(
        shiny::need(outcome_vec != '', "Please enter at least one outcome (one of 'ic50'; 'ic80'; 'iip'; 'estsens' and/or 'multsens' [for combination bnAb regimens]; or 'sens' [for single/multispecific bnAbs]) or a semicolon-separated list of outcomes (e.g., 'ic50;ic80')."),
        shiny::need(length(setdiff(outcome_vec, all_outcomes)) == 0, "You have entered one or more outcomes that are not supported at this time. Please enter at least one of the currently supported outcomes ('ic50'; 'ic80'; 'iip'; 'estsens' and/or 'multsens' [for combination bnAb regimens]; or 'sens' [for single/multispecific bnAbs]) or a semicolon-separated list of outcomes (e.g., 'ic50;ic80').")
    )
}
# check combination neutralization method
check_opts_method <- function(method_str, all_methods) {
    shiny::validate(
        shiny::need(!is.na(pmatch(method_str, all_methods)), "Please enter one of the supported methods for predicting combination neutralization. Currently, the supported methods are 'additive' and 'Bliss-Hill'.")
    )
}
# check binary outcome specification
check_opts_binary_outcome_specification <- function(bin_str, all_bins) {
    shiny::validate(
        shiny::need(bin_str %in% all_bins, "Please enter either 'ic50' or 'ic80' to use in defining binary outcomes.")
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
# check nfolds
check_opts_nfolds <- function(nfolds_str) {
    shiny::validate(
        shiny::need(nfolds_str == "" | is.numeric(as.numeric(nfolds_str)), "nfolds must be either a number (e.g., 5) or an empty string (in which case 2 folds will be used).")
    )
}
# check vimp
check_opts_vimp <- function(grp_str, ind_str, ind_type, all_grp, all_ind, all_types) {
    shiny::validate(
        shiny::need(length(setdiff(grp_str, all_grp)) == 0 | grp_str == "", "Please enter either a semicolon-separated list of supported group importance types (i.e., 'marg', 'cond', or 'marg;cond') or an empty string ('')."),
        shiny::need(length(setdiff(ind_str, all_ind)) == 0 | ind_str == "", "Please enter either a semicolon-separated list of supported individual importance types (i.e., 'marg', 'cond', 'pred', or, e.g., 'marg;cond;pred') or an empty string ('')."),
        shiny::need(length(setdiff(ind_type, all_types)) == 0 | ind_type == "",
        "Please enter either a supported individual-level intrinsic importance type (i.e., 'sitewise' or 'residuewise') or an empty string ('').")
    )
}
# check object returns
check_opts_returns <- function(return_str, all_returns) {
    shiny::validate(
        shiny::need(length(setdiff(return_str, all_returns)) == 0, "Please enter a semicolon-separated list of supported objects that you would like returned (i.e., 'report', 'learner', 'data', 'figures', 'vimp', or any semicolon-separated combination of these, e.g., 'report;learner;data;figures;vimp').")
    )
}

# check object returns
check_opts_sens_thresh <- function(sens_thresh_string) {
    sens_thresh_numeric <- suppressWarnings(as.numeric(sens_thresh_string))
    shiny::validate(
        shiny::need(!is.na(sens_thresh_numeric) & sens_thresh_numeric > 0, "Please enter a positive numeric value for sens_thresh.")
    )
}

# check object returns
check_opts_multsens_nab <- function(multsens_nab_string) {
    multsens_nab_numeric <- suppressWarnings(as.numeric(multsens_nab_string))
    shiny::validate(
        shiny::need(!is.na(multsens_nab_numeric) & multsens_nab_numeric > 0, "Please enter a positive numeric value for multsens_nab.")
    )
}

# check options
get_options_check <- function(opts) {
    # check the nab
    check_opts_nab(opts$nab)
    # check the outcome
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    check_opts_outcomes(opts$outcomes, all_outcomes, length(opts$nab))
    # check the method
    all_methods <- c("additive", "Bliss-Hill", "bliss-hill", "bh", "BH")
    check_opts_method(opts$combination_method, all_methods)
    # check the binary outcome specification (only if sens1 or sens2 is specified)
    all_bins <- c("ic50", "ic80")
    if ("sens1" %in% opts$outcomes || "sens2" %in% opts$outcomes) {
        check_opts_binary_outcome_specification(opts$binary_outcomes, all_bins)
    }
    # check the learners
    all_learners <- c("rf", "xgboost", "lasso")
    check_opts_learners(opts$learners, all_learners)
    # check cv tune
    check_opts_cvtune(opts$cvtune)
    # check cv perf
    check_opts_cvperf(opts$cvperf)
    # check importance
    all_importance_grp <- c("marg", "cond")
    all_importance_ind <- c("marg", "cond", "pred")
    all_ind_types <- c("sitewise", "residuewise")
    check_opts_vimp(opts$importance_grp, opts$importance_ind, opts$ind_importance_type, all_importance_grp, all_importance_ind, all_ind_types)
    # check objects requested for return
    all_returns <- c("report", "data", "learner", "figures", "vimp")
    check_opts_returns(opts$return, all_returns)
    check_opts_sens_thresh(opts$sens_thresh)
    check_opts_multsens_nab(opts$multsens_nab)
}
