#!/bin/bash

#-----------------------------------------------------
# This script is executed as entry point to container
#-----------------------------------------------------

# set up a log file to print out to
current_date=$(date "+%d%b%Y")
log_file_init=($(echo ${nab//'/'/'-'}"_"$current_date".log"))
log_file=($(echo "/home/output/"${log_file_init//';'/'_'}))
printf "Starting SLAPNAP \n"
printf "Messages, warnings, and errors (if any) will appear in your output directory under ${log_file//'/home/output/'/''} \n"

# make sure that user-specified options match what we expect to see
printf "Checking options \n"
echo "--- Checking options --- " > $log_file
Rscript /home/lib/check_opts.R >> $log_file 2>&1

# run script to build analytic data set
printf "Building analytic data set from CATNAP database \n"
echo "--- Building analytic data set from CATNAP database --- " >> $log_file
Rscript /home/lib/merge_proc_v4.R >> $log_file 2>&1

# run script to fit super learners
printf "Fitting super learners \n"
echo "--- Fitting super learners --- " >> $log_file
Rscript /home/lib/run_super_learners.R >> $log_file 2>&1

# run script to get variable importance
printf "Estimating variable importance \n"
echo "--- Estimating variable importance --- " >> $log_file
Rscript /home/lib/get_vimp.R >> $log_file 2>&1

# run script to compile report
printf "Compiling results using R Markdown \n"
echo "--- Compiling results using R Markdown --- " >> $log_file
Rscript /home/lib/render_report.R >> $log_file 2>&1

# return requested objects
printf "Returning requested objects \n"
echo "--- Returning requested objects --- " >> $log_file
Rscript /home/lib/return_requested_objects.R >> $log_file 2>&1

# run script
printf "Cleaning up \n"
echo "--- Cleaning up --- " >> $log_file
rm -rf /home/out/*_files

echo "--- END --- " >> $log_file
