#!/usr/bin/env Rscript

# get antibodies from environment
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(gsub("/", "-", antibody_string), split = ";")[[1]]
filename <- Sys.getenv("nab_str")

current_date <- format(as.Date(Sys.getenv('current_date'), '%d%b%Y'), "%d%b%Y")# get report name from environment
metadata_name <- paste0("metadata_", filename, "_", current_date)

# get whether to return report from environment
objects_to_return <- Sys.getenv("return")
if (objects_to_return == "" | grepl("data", objects_to_return)) {
    # name report
    rmarkdown::render('/home/lib/07_metadata.Rmd', output_file = paste0('/home/output/', metadata_name, ".html"))
} else {
    # do nothing
}
