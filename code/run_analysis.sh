#!/bin/bash

#-----------------------------------------------------
# This script is executed as entry point to container
#-----------------------------------------------------

# set up a log file to print out to
current_date=$(date "+%d%b%Y")
log_file=($(echo "/home/output/"$nab"_"$current_date".log"))
printf "Starting SLAPNAP \n"
printf "Messages, warnings, and errors (if any) will appear in $log_file \n"

# make sure that user-specified options match what we expect to see
printf "Checking options \n"
echo "--- Checking options --- \n" >> $log_file
/home/lib/check_opts.R > $log_file

# run script to build analytic data set
printf "Building analytic data set from CATNAP database \n"
echo "--- Building analytic data set from CATNAP database --- \n" >> $log_file
/home/lib/merge_proc_v4.R > $log_file

# run script to fit super learners
printf "Fitting super learners \n"
echo "--- Fitting super learners --- \n" >> $log_file
/home/lib/run_super_learners.R > $log_file

# run script to get variable importance
printf "Estimating variable importance \n"
echo "--- Estimating variable importance --- \n" >> $log_file
/home/lib/get_vimp.R > $log_file

# run script to compile report
printf "Compiling results using R Markdown \n"
echo "--- Compiling results using R Markdown --- \n" >> $log_file
/home/lib/render_report.R > $log_file

# return requested objects
printf "Returning requested objects \n"
echo "--- Returning requested objects --- \n" >> $log_file
/home/lib/return_requested_objects.R > $log_file

# run script
printf "Cleaning up \n"
echo "--- Cleaning up --- \n" >> $log_file
rm -rf /home/out/*_files

echo "--- END --- \n" >> $log_file
