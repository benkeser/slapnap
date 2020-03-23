# test utility functions
library("testthat")
library("SuperLearner")
library("tidyr")
library("dplyr")
library("tibble")

source("code/utils.R")

# set up general test opts
opts <- list()
opts$nab <- "VRC07-523-LS"
opts$outcomes <- "ic50;ic80;iip;sens1;sens2"
opts$learners <- "rf;xgboost;lasso"
opts$cvtune <- TRUE
opts$cvperf <- TRUE
opts$importance_grp <- ""
opts$importance_ind <- ""
opts$report_name <- ""
opts$return_full_sl_obj <- FALSE
opts$return_analysis_dataset <- FALSE
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
imp_df <- tibble::tibble(algo = c("rf", "lasso", "xgboost", rank = c(1, 2, 3), Importance = c(1, 2, 3)))
imp_df_lasso <- tibble::tibble(algo = c("lasso", "rf", "xgboost", rank = c(1, 2, 3), Importance = c(1, 2, 3)))
imp_df_xgb <- tibble::tibble(algo = c("xgboost", "lasso", "rf", rank = c(1, 2, 3), Importance = c(1, 2, 3)))

Sys.setenv(nab = "VRC07-523-LS", outcomes = "ic50;ic80;iip;sens1;sens2", learners = "rf;xgboost;lasso", cvtune = TRUE, cvperf = TRUE, importance_grp = "", importance_ind = "", report_name = "", return_full_sl_obj = TRUE, return_analysis_dataset = TRUE)
test_that("Getting options works", {
    global_opts <- get_global_options()
    expect_match(global_opts$nab, "VRC07-523-LS")
    expect_match(global_opts$outcomes, "ic50")
    expect_true(global_opts$cvtune)
})
test_that("Importance text works", {
    # when SL is used
    expect_output(get_importance_text(opts, imp_df), "super learner", fixed = TRUE)
    # only rf used
    #   without tuning
    expect_output(get_importance_text(opts_rf_notune, imp_df), "random forest", fixed = TRUE)
    #   with tuning
    expect_output(get_importance_text(opts_rf, imp_df), "cross-validation", fixed = TRUE)
    # only xgb used
    #   without tuning
    expect_output(get_importance_text(opts_xgb_notune, imp_df), "tree", fixed = TRUE)
    #   with tuning
    expect_output(get_importance_text(opts_xgb, imp_df), "cross-validation", fixed = TRUE)
    # only lasso used
    #   without tuning
    expect_output(get_importance_text(opts_lasso_notune, imp_df), "lasso", fixed = TRUE)
    #   with tuning
    expect_output(get_importance_text(opts_lasso, imp_df), "cross-validation", fixed = TRUE)
    # rf wins
    expect_output(get_importance_text(opts_rf, imp_df), "out-of-bag", fixed = TRUE)
    # lasso wins
    expect_output(get_importance_text(opts_lasso, imp_df_lasso), "coefficient in the final lasso fit", fixed = TRUE)
    # xgboost wins
    expect_output(get_importance_text(opts_xgb, imp_df_xgb), "xgboost gain importnce measures", fixed = TRUE)
    # what happens with no importance df
    expect_length(get_importance_text(opts), 0)
})

opts_all_import <- opts
opts_all_import$importance_grp <- "marg;cond"
opts_all_import$importance_ind <- "marg;cond"
opts_marg <- opts
opts_marg$importance_grp <- "marg"
opts_marg$importance_ind <- "marg"
opts_cond <- opts
opts_cond$importance_grp <- "cond"
opts_cond$importance_ind <- "cond"
test_that("Biological importance text works", {
    # with marginal and conditional
    expect_output(get_biological_importance_plot_description(opts_all_import, grp = TRUE), "left-hand plot", fixed = TRUE)
    # marginal only
    expect_output(get_biological_importance_plot_description(opts_marg, grp = TRUE), "marginal importance", fixed = TRUE)
    # with conditional only
    expect_output(get_biological_importance_plot_description(opts_cond, grp = TRUE), "conditional importance", fixed = TRUE)
    # with no importance
    expect_output(get_biological_importance_plot_description(opts, grp = TRUE), "")
})

ncomplete <- 100
num_obs_full <- 50
num_obs_red <- 50
test_that("Biological importance figure caption works", {
    # test outcomes
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = TRUE), "IC-50")
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic80", grp = TRUE), "IC-80")
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "iip", grp = TRUE), "IIP")
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "sens1", grp = TRUE), "estimated sensitivity")
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "sens2", grp = TRUE), "multiple sensitivity")
    # test group
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = TRUE), "Group")
    # test individual feature
    expect_output(biological_importance_figure_caption(ncomplete, num_obs_full, num_obs_red, "ic50", grp = FALSE), "Individual")
})
