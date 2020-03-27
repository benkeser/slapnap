# check options checking
library("testthat")
library("SuperLearner")
library("tidyr")
library("dplyr")
library("tibble")

source("code/utils.R")
source("code/check_opts_functions.R")

# set up general test opts
Sys.setenv(nab = "VRC07-523-LS", outcomes = "ic50;ic80;iip;sens1;sens2", learners = "rf;xgboost;lasso", cvtune = TRUE, cvperf = TRUE, importance_grp = "", importance_ind = "", report_name = "", return = "")
opts <- get_global_options()
test_that("Run with multiple learners and outcomes works", {
    expect_warning(get_options_check(opts), NA)
})

Sys.setenv(outcomes = "ic50", learners = "rf")
opts <- get_global_options()
test_that("Runs with single learners and outcomes work", {
    expect_warning(get_options_check(opts), NA)
})
Sys.unsetenv("nab")
opts <- get_global_options()
test_that("No NAb error message works", {
    expect_error(get_options_check(opts), "Please enter at least one antibody.*", class = "shiny.silent.error")
})
Sys.setenv(nab = "VRC07-523-LS")

Sys.unsetenv("outcomes")
opts <- get_global_options()
test_that("No outcome error message works", {
    expect_error(get_options_check(opts), "Please enter at least one outcome.*", class = "shiny.silent.error")
})
Sys.setenv(outcomes = "")
opts <- get_global_options()
test_that("Blank outcome error message works", {
    expect_error(get_options_check(opts), "Please enter at least one outcome.*", class = "shiny.silent.error")
})
Sys.setenv(outcomes = "not a real outcome")
opts <- get_global_options()
test_that("Not real outcome error works", {
    expect_error(get_options_check(opts), "You have entered one or more outcomes that are not supported.*", class = "shiny.silent.error")
})
Sys.setenv(outcomes = "ic50")
Sys.unsetenv("learners")
opts <- get_global_options()
test_that("No learners error message works", {
    expect_error(get_options_check(opts), "Please enter one or more learners.*", class = "shiny.silent.error")
})
Sys.setenv(learners = "glm")
opts <- get_global_options()
test_that("Unsupported learner error message works", {
    expect_error(get_options_check(opts), "You have entered one or more learners that are not supported.*", class = "shiny.silent.error")
})
Sys.setenv(learners = "rf")
Sys.unsetenv("importance_grp")
opts <- get_global_options()
test_that("No grp importance error works", {
    expect_warning(get_options_check(opts), NA)
})
Sys.setenv(importance_grp = "all")
opts <- get_global_options()
test_that("Unsupported grp importance error message works", {
    expect_error(get_options_check(opts), "Please enter either a semicolon-separated list of supported group importance types.*", class = "shiny.silent.error")
})
Sys.setenv(importance_grp = "marg")
Sys.unsetenv("importance_ind")
opts <- get_global_options()
test_that("Unset ind import works", {
    expect_warning(get_options_check(opts), NA)
})
Sys.setenv(importance_ind = "all")
opts <- get_global_options()
test_that("Unsupported in import error works", {
    expect_error(get_options_check(opts), "Please enter either a semicolon-separated list of supported individual importance types.*", class = "shiny.silent.error")
})
Sys.setenv(importance_ind = "marg")
Sys.unsetenv("report_name")
opts <- get_global_options()
test_that("Empty report works", {
    expect_warning(get_options_check(opts), NA)
})
