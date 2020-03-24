#!/bin/sh

# post-build test script

# check options checking
Rscript tests/test-check_opts.R

# check utility functions
Rscript tests/test-utils.R

# set up local directory arg
brian_computer=true
if $brian_computer; then
    local_dir=~/Projects/VIDD/hvtn/slapnap/
else
    # for David to set up
    local_dir=~/Dropbox/Emory/AMP/slapnap/
fi

# only run if necessary
# sudo docker build -t slapnap .
# ------------------------------------------
# Single nAb, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50; prints only descriptive statistics (no predictive performance) -- do we want that?
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_ic50_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_ic50_nocv_novimp" \
    slapnap

# no cv, no vimp, ic80
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_ic80_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic80" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_ic80_nocv_novimp" \
    slapnap

# no cv, no vimp, iip
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_iip_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="iip" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_iip_nocv_novimp" \
    slapnap

# no cv, no vimp, sens1
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_sens1_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="sens1" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_sens1_nocv_novimp" \
    slapnap

# no cv, no vimp, sens2
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_sens2_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="sens2" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_sens2_nocv_novimp" \
    slapnap

# no cv, vimp, ic50
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_ic50_nocv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_ic50_nocv_vimp" \
    slapnap

# cv, vimp, ic50 (~ 20 minutes)
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_ic50_cv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_ic50_cv_vimp" \
    slapnap

# cv, vimp, sens1 (~ 20 minutes)
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_sens1_cv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="sens1" \
    -e learners="lasso" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS_sens1_cv_vimp" \
    slapnap

# ------------------------------------------
# Multiple nAbs, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_ic50_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_nocv_novimp" \
    slapnap

# no cv, no vimp, iip
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_iip_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="iip" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_iip_nocv_novimp" \
    slapnap

# no cv, vimp
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_ic50_nocv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_nocv_vimp" \
    slapnap

# cv, vimp
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_ic50_cv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_cv_vimp" \
    slapnap

# ------------------------------------------
# Multiple nAbs, single outcome, multiple learners
# ------------------------------------------
# no cv, no vimp, ic50
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_multilearner_ic50_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_nocv_novimp_multilearn" \
    slapnap

# no cv, no vimp, iip
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_multilearner_iip_nocv_novimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="iip" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_iip_nocv_novimp_multilearn" \
    slapnap

# no cv, vimp
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_multilearner_ic50_nocv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_nocv_vimp_multilearn" \
    slapnap

# cv, vimp
sudo docker run \
    -v $local_dir/sandbox/dat/analysis:/home/dat/analysis/ \
    -v $local_dir/code:/home/lib/ \
    -v $local_dir/sandbox/slfits_multinab_multilearner_ic50_cv_vimp:/home/slfits/ \
    -v $local_dir/sandbox/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    -e report_name="VRC07-523-LS;PGT121_ic50_cv_vimp_multilearn" \
    slapnap
