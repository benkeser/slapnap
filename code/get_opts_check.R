#! /usr/bin/env Rscript

# check options
get_options_check(opts) {
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
    all_importance_ind <- c("marg", "cond", "pred")
    check_opts_vimp(opts$importance_grp, opts$importance_ind, all_importance_grp, all_importance_ind)
    # check objects requested for return
    all_returns <- c("report", "data", "learner", "figures", "vimp")
    check_opts_returns(opts$return, all_returns)
}
