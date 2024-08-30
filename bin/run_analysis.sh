#!/usr/bin/bash

#-----------------------------------------------------
# This script is executed as entry point to container
#-----------------------------------------------------
# check if output directory is mounted or view_port is turned on
# check_mount=$(mount | grep '/home/output')
# mount_val=$?
# if [[ mount_val -ne 0 && $view_port = "FALSE" ]]
# then
#     printf "No way to receive slapnap output. Please either mount a directory to the container directory /home/output (using -v) or set -e view_port='TRUE'. See documentation for further details."
#     exit 125
# fi
# allow errors to propagate up to container
set -e
# grab the current date, promote to system environment variable
current_date=$(date "+%d%b%Y")
export current_date


default_outdir=$(pwd)

export output="${default_outdir}/output/"
export analysis="${default_outdir}/dat/analysis/"
export catnap="${default_outdir}/dat/catnap/"
export slfits="${default_outdir}/slfits/"

# make the nab string suitable for naming files, promote to system env var
remove_slashes=${nab//'/'/'-'}
nab_str=${remove_slashes//';'/'_'}
export nab_str
# set up a log file to print out to
log_file_init=($(echo $nab_str"_"$current_date".log"))
log_file=($(echo "${output}${log_file_init}"))
# start SLAPNAP
printf "Starting SLAPNAP \n"
printf "Messages, warnings, and errors (if any) will appear in your output directory under $log_file \n"

# make sure that user-specified options match what we expect to see
printf "Checking options \n"
echo "--- Checking options --- " > $log_file


01_check_opts.R >> $log_file 2>&1

# run script to build analytic data set
printf "Building analytic data set from CATNAP database \n"
echo "--- Building analytic data set from CATNAP database --- " >> $log_file
02_compile_analysis_dataset.R >> $log_file 2>&1

# run script to fit super learners
# but only fit if something other than just data is requested as output
if [[ "$return" == *"report"* ]] || [[ "$return" == *"learner"* ]] || [[ "$return" == *"vimp"* ]] || [[ "$return" == *"figures"* ]]
then
    printf "Fitting learners \n"
    echo "--- Fitting learners --- " >> $log_file
    03_run_super_learners.R >> $log_file 2>&1

    VIMP_REGEX="^(marg|cond|marg;cond)"
    if [[ "$importance_grp" =~ $VIMP_REGEX ]] || [[ "$importance_ind" =~ $VIMP_REGEX ]]
    then
        # run script to get variable importance
        printf "Estimating variable importance \n"
        echo "--- Estimating variable importance --- " >> $log_file
        04_get_vimp.R >> $log_file 2>&1
    fi
fi

if [[ "$return" == *"report"* ]] || [[ "$return" == *"figures"* ]]
then
    # run script to compile report
    printf "Compiling results using R Markdown \n"
    echo "--- Compiling results using R Markdown --- " >> $log_file
    05_render_report.R >> $log_file 2>&1
fi

# return requested objects
printf "Returning requested objects \n"
echo "--- Returning requested objects --- " >> $log_file
06_return_requested_objects.R >> $log_file 2>&1

if [[ "$return" == *"data"* ]]
then
    printf "Generating metadata using R Markdown \n"
    echo "--- Generating metadata using R Markdown --- " >> $log_file
    07_render_metadata.R >> $log_file 2>&1
fi

# if requested, port
if [[ "$view_port" == "TRUE" ]] && [[ "$return" == *"report"* ]]
then
    printf "Report can be viewed on localhost. To stop container, retrieve CONTAINER ID using 'docker container ps' and 'docker stop CONTAINER_ID'.  \n"
    echo "--- Report can be viewed on localhost. ---" >> $log_file
    # copy report to www folder for viewing
    name_of_report=$(08_generate_report_name.R)
    cp ${output}${name_of_report}.html /var/www/html/index.html
    nginx -g "daemon off;"
fi

printf "Closing down SLAPNAP. Check your output directory for requested objects.\n"
echo "--- END --- " >> $log_file
