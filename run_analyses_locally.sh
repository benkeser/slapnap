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
sudo docker run \
    -v $local_dir/local_production_runs/slfits_vrc01:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC01" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return_full_sl_obj="TRUE" \
    -e return_analysis_dataset="TRUE" \
    slapnap

# update to VRC07-523LS
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
    -e return_full_sl_obj="TRUE" \
    -e return_analysis_dataset="TRUE" \
    slapnap bash

# -------------------------------------------
# 2-bnAb regimens
# -------------------------------------------
# HVTN 130
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
    -e return_full_sl_obj="TRUE" \
    -e return_analysis_dataset="TRUE" \
    slapnap bash

# -------------------------------------------
# 3-bnAb regimens
# -------------------------------------------
# HVTN 130
sudo docker run -it \
    -v $local_dir/local_production_runs/slfits_vrc07-523-ls_pgt121_pgdm1400:/home/slfits/ \
    -v $local_dir/local_production_runs/output:/home/output/ \
    -e nab="VRC07-523-LS;PGT121;PGDM1400" \
    -e outcomes="ic50;ic80;iip;sens1;sens2" \
    -e learners="rf;lasso;xgboost" \
    -e cvtune="TRUE" \
    -e cvperf="TRUE" \
    -e importance_grp="marg" \
    -e importance_ind="marg;pred" \
    -e return_full_sl_obj="TRUE" \
    -e return_analysis_dataset="TRUE" \
    slapnap bash
