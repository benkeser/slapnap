#! /usr/bin/env Rscript

# if any objects are requested, return them
source("/home/lib/utils.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

#------------------------------------
# save analysis dataset, if requested
#------------------------------------
if (grepl("data", opts$return)) {
    analysis_data_names <- list.files("/home/dat/analysis")
    analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
    file.copy(paste0("/home/dat/analysis/", analysis_data_name), paste0("/home/output/", analysis_data_name))
}
#--------------------------------------
# save full super learner, if requested
#--------------------------------------
if (grepl("learner", opts$return)) {
    all_fit_names <- list.files("/home/slfits")
    fit_names <- get_learner_fit_names(all_fit_names, opts)
    file.copy(paste0("/home/slfits/", fit_names), paste0("/home/output/", fit_names))
}
#------------------------------------
# save figures, if requested
#------------------------------------
if (grepl("figures", opts$return)) {
    
}
#------------------------------------
# save variable importance objects, if requested
#------------------------------------
if (grepl("vimp", opts$return)) {
    all_fit_names <- list.files("/home/slfits")
    vimp_names <- get_vimp_object_names(all_fit_names, opts)
    file.copy(paste0("/home/slfits/", vimp_names), paste0("/home/output/", vimp_names))
}
