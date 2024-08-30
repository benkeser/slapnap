#!/usr/bin/env -S Rscript --vanilla

antibody_string <- Sys.getenv("nab")
report_name <- Sys.getenv("report_name")

current_date <- format(Sys.time(), "%d%b%Y")

antibodies <- strsplit(
    gsub("/", "-", antibody_string),
    split = ";"
)[[1]]

if (report_name == "") {
    report_name <- paste0(
        "report_",
        paste(antibodies, collapse = "_"),
        "_",
        current_date
    )
}
cat(report_name)