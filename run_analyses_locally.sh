#!/bin/sh

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

# -------------------------------------------
# Single-bnAb regimens
# -------------------------------------------
# replicate the Magaret et al. (2019) analysis; for HVTN 703, HVTN704
sudo docker run -it \
    -v $local_dir/local_production_runs/slfits_vrc01:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC01" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap bash

# update to VRC07-523LS (HVTN ##)
sudo docker run -it \
    -v $local_dir/local_production_runs/slfits_vrc07_523ls:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap bash

# PGT121 (HVTN 136)
sudo docker run -it \
    -v $local_dir/local_production_runs/slfits_pgt121:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="PGT121" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap bash

# -------------------------------------------
# 2-bnAb regimens
# -------------------------------------------
# HVTN 130 comparison 1
sudo docker run -it \
    -v $local_dir/local_production_runs/slfits_vrc07-523-ls_pgt121:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap bash

# HVTN 130 comparison 2
sudo docker run \
    -v $local_dir/local_production_runs/slfits_vrc07-523-ls_pgdm1400:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS;PGDM1400" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap

# HVTN 130 comparison 3
sudo docker run \
    -v $local_dir/local_production_runs/slfits_vrc07-523-ls_10-1074:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS;10-1074" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap


# -------------------------------------------
# 3-bnAb regimens
# -------------------------------------------
# HVTN 129
sudo docker run \
    -v $local_dir/local_production_runs/slfits_vrc01-10e8v4-pgdm1400:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC01/PGDM1400-10E8v4" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap


# HVTN 130
sudo docker run \
    -v $local_dir/local_production_runs/slfits_vrc07-523-ls_pgt121_pgdm1400:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121;PGDM1400" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return="report;data" \
    slapnap
