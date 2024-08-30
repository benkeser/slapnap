#! /usr/bin/env -S Rscript --vanilla

library("slapnap") # source our function library
path.data.catnap <- Sys.getenv("catnap") #file.path(path.data, "catnap")

# if any objects are requested, return them

path.ml.slfits <- Sys.getenv("slfits") #file.path(path.data, "analysis")
path.data.catnap <- Sys.getenv("catnap") #file.path(path.data, "catnap")
path.data.analysis <- Sys.getenv("analysis") #file.path(path.data, "analysis")
path.out <- Sys.getenv("output") #file.path(path.home, "output")

# source(file.path(path.lib, "00_utils.R"))

# source("/home/lib/00_utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()
filename <- Sys.getenv("nab_str")
postfix <- paste0(filename, "_", format(as.Date(Sys.getenv('current_date'), "%d%b%Y"), "%d%b%Y"))

h2o_here <- !(all(grepl("h2oboost", opts$learners) == FALSE))

#-------------------------------------------------
# save analysis dataset, if requested
#-------------------------------------------------
if (any(grepl("data", opts$return))) {
    data_dir <- file.path(path.data.analysis)
    analysis_data_names <- list.files(data_dir)
    analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
    init_data <- readRDS(paste0(data_dir, analysis_data_name))
    analysis_dataset <- clean_analysis_dataset(init_data, opts)
    readr::write_csv(analysis_dataset, file = paste0(path.out, gsub(".rds", ".csv", analysis_data_name)))
}
#--------------------------------------
# save full super learner, if requested
#--------------------------------------
if (any(grepl("learner", opts$return))) {
    all_fit_names <- list.files(path.ml.slfits)
    fit_names <- get_learner_fit_names(all_fit_names, opts)
    fit_renames <- gsub("dichotomous.1", ifelse(length(opts$nab) == 1, "sens", "estsens"), fit_names)
    fit_renames <- gsub("dichotomous.2", "multsens", fit_renames)
    fit_renames <- gsub("log10.pc.ic50", "ic50", fit_renames)
    fit_renames <- gsub("log10.pc.ic80", "ic80", fit_renames)
    fit_renames <- gsub("fit_", "learner_", fit_renames)
    fit_renames <- gsub(".rds", paste0("_", postfix, ".rds"), fit_renames)
    file.copy(paste0(path.ml.slfits, fit_names), paste0(path.out, fit_renames), overwrite = TRUE)
    if (h2o_here) {
        # copy h2o model files as well
        h2o_fit_names <- all_fit_names[grepl("GBM_model", all_fit_names)]
        file.copy(paste0(path.ml.slfits, h2o_fit_names), paste0(path.out, h2o_fit_names), overwrite = TRUE)
    }
}
#------------------------------------
# figures are saved directly from
# new_report.Rmd
#------------------------------------

#------------------------------------
# save variable importance objects, if requested
#------------------------------------
if (any(grepl("vimp", opts$return))) {
    all_fit_names <- list.files(path.ml.slfits)
    vimp_names <- get_vimp_object_names(all_fit_names, opts)
    vimp_renames <- gsub("dichotomous.1", ifelse(length(opts$nab) == 1, "sens", "estsens"), vimp_names)
    vimp_renames <- gsub("dichotomous.2", "multsens", vimp_renames)
    vimp_renames <- gsub("log10.pc.ic50", "ic50", vimp_renames)
    vimp_renames <- gsub("log10.pc.ic80", "ic80", vimp_renames)
    vimp_renames <- gsub(".rds", paste0("_", postfix, ".rds"), vimp_renames)
    file.copy(paste0(path.ml.slfits, vimp_names), paste0(path.out, vimp_renames), overwrite = TRUE)
}
