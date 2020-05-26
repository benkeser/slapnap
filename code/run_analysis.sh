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
# but only fit if something other than just data is requested as output
if [[ "$return" == *"report"* ]] || [[ "$return" == *"learner"* ]] || [[ "$return" == *"vimp"* ]] || [[ "$return" == *"figures"* ]]
then
  printf "Fitting learners \n"
  echo "--- Fitting learners --- " >> $log_file
  Rscript /home/lib/run_super_learners.R >> $log_file 2>&1


  # run script to get variable importance
  printf "Estimating variable importance \n"
  echo "--- Estimating variable importance --- " >> $log_file
  Rscript /home/lib/get_vimp.R >> $log_file 2>&1
fi

if [[ "$return" == *"report"* ]] || [[ "$return" == *"figures"* ]]
then
  # run script to compile report
  printf "Compiling results using R Markdown \n"
  echo "--- Compiling results using R Markdown --- " >> $log_file
  Rscript /home/lib/render_report.R >> $log_file 2>&1
fi

# return requested objects
printf "Returning requested objects \n"
echo "--- Returning requested objects --- " >> $log_file
Rscript /home/lib/return_requested_objects.R >> $log_file 2>&1

# if requested, port
if [[ "$view_port" == "TRUE" ]] && [[ "$return" == *"report"* ]]
then
printf "Report can be viewed on localhost. To stop container, retrieve CONTAINER ID using 'docker container ps' and 'docker stop CONTAINER_ID'.  \n"
echo "--- Report can be viewed on localhost. ---" >> $log_file
cp /home/output/*.html /var/www/html/index.html
nginx -g "daemon off;"
fi

echo "--- END --- " >> $log_file
