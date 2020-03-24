#!/bin/sh

# post-build test script

# check options checking
Rscript tests/test-check_opts.R

# check utility functions
Rscript tests/test-utils.R

# ------------------------------------------
# Single nAb, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50
sudo docker run \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_ic50_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, no vimp, sens1
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_sens1_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="sens1" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, vimp, ic50
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_ic50_nocv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# cv, vimp, ic50
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_ic50_cv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# ------------------------------------------
# Multiple nAbs, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_ic50_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, no vimp, iip
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_iip_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="iip" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, vimp
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_ic50_nocv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# cv, vimp
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_ic50_cv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# ------------------------------------------
# Multiple nAbs, single outcome, multiple learners
# ------------------------------------------
# no cv, no vimp, ic50
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_multilearner_ic50_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, no vimp, iip
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_multilearner_iip_nocv_novimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="iip" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="" \
    -e importance_ind="" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# no cv, vimp
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_multilearner_ic50_nocv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="FALSE" \
    -e cvperf="FALSE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap

# cv, vimp
sudo docker run -v ~/Projects/VIDD/hvtn/slapnap/sandbox:/home/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/dat/analysis:/home/dat/analysis/ \
    -v ~/Projects/VIDD/hvtn/slapnap/code:/home/lib/ \
    -v ~/Projects/VIDD/hvtn/slapnap/sandbox/slfits_multinab_multilearner_ic50_cv_vimp:/home/slfits/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap
