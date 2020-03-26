# test ml importance extraction
library("testthat")
library("SuperLearner")
library("tidyr")
library("dplyr")
library("tibble")

source("code/utils.R")
source("code/super_learner_libraries.R")
source("code/plotting_functions.R")
source("code/ml_var_importance_measures.R")

# set up general options
Sys.setenv(nab = "VRC07-523-LS", outcomes = "ic50", learners = "lasso", cvtune = FALSE, cvperf = FALSE, importance_grp = "", importance_ind = "", report_name = "", return_full_sl_obj = FALSE, return_analysis_dataset = FALSE)
opts <- get_global_options()
# generate some test data
set.seed(4747)
x <- replicate(10, rnorm(100, 0, 1))
y <- rbinom(100, 1, exp(x)/(1+exp(x)))
dat <- data.frame(dichotomous.1 = y, x)
dat_continuous <- data.frame(VRC07.523.LS.ic50.imputed = rnorm(100, 0, 1), x)
# !cvtune & !cvperf
sl_fit_nocv <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = "SL.glmnet.0", cvControl = list(V = 1), method = "tmp_method.CC_nloglik", fit_name = "fit_nocv.rds", cv_fit_name = "cvfit_nocv.rds")

test_that("!cvtune and !cvperf works with one learner", {
    fit_sl <- readRDS("tests/fit_nocv.rds")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})

Sys.setenv(learners="rf;lasso")
opts <- get_global_options()
# picks the first one in library
sl_fit_nocv_twolearners <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = c("SL.glmnet.0", "SL.ranger.small"), cvControl = list(V = 2), method = "tmp_method.CC_nloglik", fit_name = "fit_nocv_multilearner.rds", cv_fit_name = "cvfit_nocv_multilearner.rds")
test_that("!cvtune and !cvperf works with multiple learners", {
    fit_sl <- readRDS("tests/fit_nocv_multilearner.rds")
    expect_equal(class(fit_sl), "ranger")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})
Sys.setenv(learners="lasso")

Sys.setenv(cvperf = TRUE)
opts <- get_global_options()
sl_fit_nocvtune <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = "SL.glmnet.0", cvControl = list(V = 5), method = "tmp_method.CC_nloglik", fit_name = "fit_nocvtune.rds", cv_fit_name = "cvfit_nocvtune.rds")
test_that("!cvtune and cvperf works", {
    fit_sl <- readRDS("tests/fit_nocvtune.rds")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})

Sys.setenv(cvperf = FALSE, cvtune = TRUE)
opts <- get_global_options()
sl_fit_nocvperf <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = c("SL.glmnet.0", "SL.glmnet.25"), cvControl = list(V = 5), method = "tmp_method.CC_nloglik", fit_name = "fit_nocvperf.rds", cv_fit_name = "cvfit_nocvperf.rds")
test_that("cvtune and !cvperf works", {
    fit_sl <- readRDS("tests/fit_nocvperf.rds")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})

Sys.setenv(cvperf = TRUE)
opts <- get_global_options()
sl_fit_nocvtune <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = "SL.glmnet.0", cvControl = list(V = 5), method = "tmp_method.CC_nloglik", fit_name = "fit_cv.rds", cv_fit_name = "cvfit_cv.rds")
test_that("cvtune and cvperf works", {
    fit_sl <- readRDS("tests/fit_cv.rds")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})

Sys.setenv(learners="rf;lasso")
opts <- get_global_options()
sl_fit_twolearners <- sl_one_outcome(dat, outcome_name = "dichotomous.1", pred_names = paste0("X", 1:10), opts = opts, save_dir = "tests/", family = binomial(), SL.library = c("SL.glmnet.0", "SL.ranger.small"), cvControl = list(V = 2), method = "tmp_method.CC_nloglik", fit_name = "fit_cv_multilearner.rds", cv_fit_name = "cvfit_cv_multilearner.rds")
test_that("cvtune and cvperf works with multiple learners", {
    fit_sl <- readRDS("tests/fit_cv_multilearner.rds")
    import_df <- extract_importance(fit_sl, opts)
    expect_equal(ncol(import_df), 4)
    expect_equal(nrow(import_df), 10)
})
