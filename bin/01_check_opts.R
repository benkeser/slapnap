#!/usr/bin/env -S Rscript --vanilla

library("shiny")
library("slapnap")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

# --------------------------------------
# make sure that the options are in the form we expected
# --------------------------------------
get_options_check(opts)
