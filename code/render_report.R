#!/usr/bin/env Rscript

# get antibodies from environment
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]

# name report
report_name <- paste0 ("report_", paste (antibodies, collapse = "_"), "_", format (Sys.time (), "%d%b%Y"), ".html")

# rmarkdown::render('/home/lib/report.Rmd', output_file = paste0('/home/slfits/', report_name))
rmarkdown::render('/home/lib/new_report.Rmd', output_file = paste0('/home/slfits/', report_name))
