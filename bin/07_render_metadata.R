#!/usr/bin/env -S Rscript --vanilla

library("slapnap")

path.out <- Sys.getenv("output")
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
    rmarkdown::render(system.file("rmd", "07_metadata.Rmd", package = "slapnap"), output_file = paste0(path.out, metadata_name, ".html"))
} else {
    # do nothing
}
