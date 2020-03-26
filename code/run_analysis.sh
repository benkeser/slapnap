#!/bin/bash

#-----------------------------------------------------
# This script is executed as entry point to container
#-----------------------------------------------------

# make sure that user-specified options match what we expect to see
printf "Checking options \n"
/home/lib/check_opts.R

# run script to build analytic data set
printf "Building analytic data set from CATNAP database \n"
/home/lib/merge_proc_v4.R

# run script to fit super learners
printf "Fitting super learners \n"
/home/lib/run_super_learners.R

# run script to get variable importance
printf "Estimating variable importance \n"
/home/lib/get_vimp.R

# run script to compile report
printf "Compiling results using R Markdown \n"
/home/lib/render_report.R

# run script
printf "Cleaning up \n"
rm -rf /home/out/*_files
