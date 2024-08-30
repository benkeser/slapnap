#!/usr/bin/env -S Rscript --vanilla

# get antibodies from environment
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(gsub("/", "-", antibody_string), split = ";")[[1]]
filename <- Sys.getenv("nab_str")
path.out <- Sys.getenv("output") #file.path(path.home, "output")

report_name <- Sys.getenv("report_name")
current_date <- format(as.Date(Sys.getenv('current_date'), '%d%b%Y'), "%d%b%Y")# get report name from environment
if(report_name == ""){
	report_name <- paste0(path.out, "report_", filename, "_", current_date)
}

# get whether to return report from environment
objects_to_return <- Sys.getenv("return")
if (objects_to_return == "" | grepl("report", objects_to_return)) {
    # name report
    rmarkdown::render(system.file("rmd", "05_report.Rmd", package = "slapnap"), output_file = paste0(report_name, ".html"))
} else {
    rmarkdown::render(system.file("rmd", "05_report.Rmd", package = "slapnap"))
}
