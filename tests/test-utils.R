# test utility functions
library("testthat")
library("SuperLearner")
library("tidyr")
library("dplyr")
library("tibble")

source("code/utils.R")
source("code/super_learner_libraries.R")
source("code/plotting_functions.R")

# set up general test opts
opts <- list()
opts$nab <- "VRC07-523-LS"
opts$outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
opts$learners <- c("rf", "lasso", "xgboost")
opts$cvtune <- TRUE
opts$cvperf <- TRUE
opts$importance_grp <- ""
opts$importance_ind <- ""
opts$report_name <- ""
opts$return <- FALSE
opts_rf <- opts
opts_rf$learners <- "rf"
opts_rf_notune <- opts_rf
opts_rf_notune$cvtune <- FALSE
opts_lasso <- opts
opts_lasso$learners <- "lasso"
opts_lasso_notune <- opts_lasso
opts_lasso_notune$cvtune <- FALSE
opts_xgb <- opts
opts_xgb$learners <- "xgboost"
opts_xgb_notune <- opts_xgb
opts_xgb_notune$cvtune <- FALSE
# set up general importance df
imp_df <- tibble::tibble(algo = c("rf", "lasso", "xgboost"), rank = c(1, 2, 3), Importance = c(1, 2, 3))
imp_df_lasso <- tibble::tibble(algo = c("lasso", "rf", "xgboost"), rank = c(1, 2, 3), Importance = c(1, 2, 3))
imp_df_xgb <- tibble::tibble(algo = c("xgboost", "lasso", "rf"), rank = c(1, 2, 3), Importance = c(1, 2, 3))

Sys.setenv(nab = "VRC07-523-LS", outcomes = "ic50;ic80;iip;sens1;sens2", learners = "rf;xgboost;lasso", cvtune = TRUE, cvperf = TRUE, importance_grp = "", importance_ind = "", report_name = "", return = TRUE)
global_opts <- get_global_options()
test_that("Getting options works", {
    expect_match(global_opts$nab, "VRC07-523-LS")
    expect_equal(global_opts$outcomes, c("ic50", "ic80", "iip", "sens1", "sens2"))
    expect_true(global_opts$cvtune)
})
test_that("Importance text works", {
    # when SL is used
    expect_match(get_importance_text(opts, imp_df), "super learner", fixed = TRUE)
    # only rf used
    #   without tuning
    expect_match(get_importance_text(opts_rf_notune, imp_df), "random forest", fixed = TRUE)
    #   with tuning
    expect_match(get_importance_text(opts_rf, imp_df), "super learner ensemble", fixed = TRUE)
    # only xgb used
    #   without tuning
    expect_match(get_importance_text(opts_xgb_notune, imp_df_xgb), "tree", fixed = TRUE)
    #   with tuning
    expect_match(get_importance_text(opts_xgb, imp_df_xgb), "super learner ensemble", fixed = TRUE)
    # only lasso used
    #   without tuning
    expect_match(get_importance_text(opts_lasso_notune, imp_df_lasso), "lasso", fixed = TRUE)
    #   with tuning
    expect_match(get_importance_text(opts_lasso, imp_df_lasso), "cross-validation", fixed = TRUE)
    # rf wins
    expect_match(get_importance_text(opts, imp_df), "out-of-bag", fixed = TRUE)
    # lasso wins
    expect_match(get_importance_text(opts, imp_df_lasso), "coefficient in the final lasso fit", fixed = TRUE)
    # xgboost wins
    expect_match(get_importance_text(opts, imp_df_xgb), "xgboost gain importance measures", fixed = TRUE)
    # what happens with no importance df
    expect_equal(get_importance_text(opts), "")
})

opts_all_import <- opts
opts_all_import$importance_grp <- c("marg", "cond")
opts_all_import$importance_ind <- c("marg", "cond")
opts_marg <- opts
opts_marg$importance_grp <- "marg"
opts_marg$importance_ind <- "marg"
opts_cond <- opts
opts_cond$importance_grp <- "cond"
opts_cond$importance_ind <- "cond"
test_that("Biological importance text works", {
    # with marginal and conditional
    expect_match(get_biological_importance_plot_description(opts_all_import, grp = TRUE), "left-hand plot", fixed = TRUE)
    # marginal only
    expect_match(get_biological_importance_plot_description(opts_marg, grp = TRUE), "marginal biological importance", fixed = TRUE)
    # with conditional only
    expect_match(get_biological_importance_plot_description(opts_cond, grp = TRUE), "conditional biological importance", fixed = TRUE)
    # with no importance
    expect_equal(get_biological_importance_plot_description(opts, grp = TRUE), "")
})

ncomplete <- 100
num_obs_full <- 50
num_obs_red <- 50
test_that("Biological importance figure caption works", {
    # test outcomes
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = TRUE), "IC-50")
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic80", grp = TRUE), "IC-80")
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "iip", grp = TRUE), "IIP")
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "sens1", grp = TRUE), "estimated sensitivity")
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "sens2", grp = TRUE), "multiple sensitivity")
    # test group
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = TRUE), "Group")
    # test individual feature
    expect_match(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = FALSE), "Individual")
})

# generate some test data
set.seed(4747)
x <- replicate(10, rnorm(100, 0, 1))
y <- rbinom(100, 1, exp(x)/(1+exp(x)))
dat <- data.frame(dichotomous.1 = y, x)
dat_continuous <- data.frame(VRC07.523.LS.ic50.imputed = rnorm(100, 0, 1), x)
opts_dichot1 <- opts
opts_dichot1$outcomes <- "sens1"
opts_ic50 <- opts
opts_ic50$outcomes <- "ic50"
test_that("Individual nab summaries work", {
    nab_summ <- get_individual_nab_summaries(outcome = "ic50", opts = opts_ic50, dat = dat_continuous)
    expect_equal(length(nab_summ$hist[[1]]$labels), 4)
    expect_equal(length(nab_summ$summary[[1]]), 7)
})

test_that("Learner descriptions work", {
    expect_match(get_learner_descriptions(opts_rf), "with tuning parameters selected")
    expect_match(get_learner_descriptions(opts_rf_notune), "with tuning parameters set to their")
    expect_match(get_learner_descriptions(opts), "a super learner ensemble of several")
})

sl_fit <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = "SL.glmnet.0", cvControl = list(V = 5), method = "tmp_method.CC_nloglik")
sl_fit <- readRDS("tests/fit_dichotomous.1.rds")
sl_fitted <- readRDS("tests/fitted_dichotomous.1.rds")
test_that("Run SL works", {
    expect_equal(as.numeric(sl_fit$coef), 1)
})

sl_fit_lst <- load_cv_fits(opts_dichot1, "tests/")
test_that("CV outcome tables works", {
    expect_equal(length(sl_fit_lst$out), 1)
    cv_outcome_tbl_lst <- get_cv_outcomes_tables(sl_fit_lst, opts_dichot1)
    expect_equal(length(cv_outcome_tbl_lst$r2), 0)
    expect_output(print(cv_outcome_tbl_lst$auc), "CVAUC")
})

test_that("Outcome descriptions work", {
    expect_match(get_outcome_descriptions(opts), "IIP is calculated as")
    expect_match(get_outcome_descriptions(opts), "Estimated sensitivity is defined by")
    expect_match(get_outcome_descriptions(opts), "Since only one antibody was specified for this analysis")
    opts_2nab <- opts
    opts_2nab$nab <- c("VRC07-523-LS", "PGT121")
    expect_match(get_outcome_descriptions(opts_2nab), "Predicted IC-50")
    expect_match(get_outcome_descriptions(opts_2nab), "additive model of Wagh et al")
    expect_match(get_outcome_descriptions(opts_2nab), "Multiple sensitivity is defined as the binary indicator of")
})

test_that("Comma separated outcomes works", {
    expect_equal(get_comma_sep_outcomes(opts), "IC-50, IC-80, IIP, estimated sensitivity, multiple sensitivity.")
    expect_equal(get_comma_sep_outcomes(list(outcomes = "ic50")), "IC-50.")
})

set.seed(4747)
test_that("Making CV folds works", {
    y <- rbinom(100, 1, 0.5)
    V <- 5
    folds_1 <- make_folds(y, V, stratified = TRUE)
    folds_2 <- make_folds(y, V, stratified = FALSE)
    expect_equal(sum(folds_1 == 1), 21)
    expect_equal(sum(folds_2 == 1), 20)
})
test_that("Getting SL folds works", {
    folds_lst <- list(sample(1:5, 100, replace = TRUE), sample(1:5, 100, replace = TRUE))
    expect_length(get_cv_folds(folds_lst), 200)
})
test_that("SL opts work", {
    dichot_opts <- get_sl_options("dichot.1", 5)
    cont_opts <- get_sl_options("log10.pc.ic50", 5)
    expect_match(dichot_opts$method, "CC_nloglik")
    expect_match(cont_opts$method, "CC_LS")
    expect_equal(dichot_opts$fam, binomial())
    expect_equal(cont_opts$fam, gaussian())
})
test_that("VIM options work", {
    expect_match(get_vimp_options("ic50")$vimp_measure, "r_squared")
    expect_match(get_vimp_options("dichot.1")$vimp_measure, "auc")
})
test_that("VIM list creation works", {
    folds_lst <- list(sample(1:5, 100, replace = TRUE), sample(1:5, 100, replace = TRUE))
    expect_length(make_vimp_list(), 4)
    test_cv_lsts <- make_cv_lists(folds_lst, rnorm(200, 0, 1), rnorm(200, 0, 1))
    expect_length(test_cv_lsts, 3)
    expect_length(test_cv_lsts$folds, 200)
    expect_length(test_cv_lsts$full_lst, 2)
    expect_length(test_cv_lsts$redu_lst, 2)
})
test_that("Nice naming works for vimp", {
    expect_match(vimp_nice_group_names("cysteines"), "Cysteine counts")
    expect_length(vimp_nice_group_names("gp160"), 0)
    expect_equal(vimp_nice_ind_names("hxb2.K.160.1mer"), "K.160")
    expect_equal(vimp_plot_name("log10.pc.ic50"), "IC-50")
    expect_equal(vimp_plot_name("dichotomous.2"), "Multiple sensitivity")
    ## generate the data
    ## generate X
    p <- 10
    n <- 100
    x <- data.frame(replicate(p, stats::runif(n, -1, 1)))

    ## apply the function to the x's
    f <- function(x) 0.5 + 0.3*x[1] + 0.2*x[2]
    smooth <- apply(x, 1, function(z) f(z))

    ## generate Y ~ Normal (smooth, 1)
    y <- matrix(rbinom(n, size = 1, prob = smooth))

    ## set up a library for SuperLearner
    learners <- "SL.gam"

    ## using Y and X
    folds <- sample(rep(seq_len(2), length = length(y)))
    est <- vimp::vim(y, x, indx = c(1,2), type = "r_squared",
               alpha = 0.05, run_regression = TRUE,
               SL.library = learners, cvControl = list(V = 10),
               folds = folds)
    est_2 <- vimp::vim(y, x, indx = c(3, 4, 5, 6, 7, 8, 9, 10), type = "r_squared",
                     alpha = 0.05, run_regression = TRUE,
                     SL.library = learners, cvControl = list(V = 10),
                     folds = folds)
    all_est <- vimp::merge_vim(est, est_2)
    expect_equal(vimp_nice_rownames(all_est, cv = TRUE), c("NA_NA_2", "NA_NA_NA_est"))
})
