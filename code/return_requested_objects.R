#! /usr/bin/env Rscript

# if any objects are requested, return them
source("/home/lib/utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

#-------------------------------------------------
# save analysis dataset, if requested
#-------------------------------------------------
if (any(grepl("data", opts$return))) {
    data_dir <- "/home/dat/analysis/"
    analysis_data_names <- list.files(data_dir)
    analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
    init_data <- read.csv(paste0(data_dir, analysis_data_name), stringsAsFactors = FALSE)
    analysis_dataset <- clean_analysis_dataset(init_data, opts)
    write.csv(analysis_dataset, file = paste0("/home/output/", analysis_data_name), row.names = FALSE)
}
#--------------------------------------
# save full super learner, if requested
#--------------------------------------
if (any(grepl("learner", opts$return))) {
    all_fit_names <- list.files("/home/slfits")
    fit_names <- get_learner_fit_names(all_fit_names, opts)
    fit_renames <- gsub("dichotomous.1", ifelse(length(opts$nab) == 1, "sens", "estsens"), fit_names)
    fit_renames <- gsub("dichotomous.2", "multsens", fit_renames)
    fit_renames <- gsub("log10.pc.ic50", "ic50", fit_renames)
    fit_renames <- gsub("log10.pc.ic80", "ic80", fit_renames)
    fit_renames <- gsub("fit_", "learner_", fit_renames)
    file.rename(paste0("/home/slfits/", fit_names), paste0("/home/slfits/", fit_renames))
    file.copy(paste0("/home/slfits/", fit_renames), paste0("/home/output/", fit_renames))
}
#------------------------------------
# figures are saved directly from
# new_report.Rmd
#------------------------------------

#------------------------------------
# save variable importance objects, if requested
#------------------------------------
if (any(grepl("vimp", opts$return))) {
    all_fit_names <- list.files("/home/slfits")
    vimp_names <- get_vimp_object_names(all_fit_names, opts)
    file.copy(paste0("/home/slfits/", vimp_names), paste0("/home/output/", vimp_names))
}
