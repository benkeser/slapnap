#!/bin/bash 

#-----------------------------------------------------
# This script is executed as entry point to container
#-----------------------------------------------------

# run script to build analytic data set
printf "Building analytic data set from CATNAP database \n"
/home/lib/merge_proc_v2.R

# run script to fit super learners
printf "Fitting super learners \n"
/home/lib/run_super_learners.R

# run script to compile report
printf "Compiling results using R Markdown \n"
/home/lib/render_report.R

# run script
printf "Cleaning up"
rm -rf /home/out/*_files