#!/usr/bin/env Rscript

# get antibodies from environment
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]

report_name <- Sys.getenv("report_name")
# get report name from environment
if(report_name == ""){
	report_name <- paste0("report_", paste (antibodies, collapse = "_"), "_", format (Sys.time (), "%d%b%Y"))
}

# name report
rmarkdown::render('/home/lib/new_report.Rmd', output_file = paste0('/home/output/', report_name, ".html"))