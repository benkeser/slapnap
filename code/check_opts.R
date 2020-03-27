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
get_options_check()
