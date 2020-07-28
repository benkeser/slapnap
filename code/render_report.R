#!/usr/bin/env Rscript

# get antibodies from environment
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(gsub("/", "-", antibody_string), split = ";")[[1]]

report_name <- Sys.getenv("report_name")
current_date <- format(as.Date(Sys.getenv('current_date'), '%d%b%Y'), "%d%b%Y")# get report name from environment
if(report_name == ""){
	report_name <- paste0("report_", paste (antibodies, collapse = "_"), "_", current_date)
}

# get whether to return report from environment
objects_to_return <- Sys.getenv("return")
if (objects_to_return == "" | grepl("report", objects_to_return)) {
    # name report
    rmarkdown::render('/home/lib/new_report.Rmd', output_file = paste0('/home/output/', report_name, ".html"))
} else {
    rmarkdown::render('/home/lib/new_report.Rmd')
}
