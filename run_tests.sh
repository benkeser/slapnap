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
docker run -e nab="VRC07-523-LS" \
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
docker run -e nab="VRC07-523-LS" \
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
docker run -e nab="VRC07-523-LS" \
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
docker run -e nab="VRC07-523-LS" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
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
docker run -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50" \
    -e learners="lasso;rf" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="pred" \
    -e return_full_sl_obj="FALSE" \
    -e return_analysis_dataset="FALSE" \
    slapnap
