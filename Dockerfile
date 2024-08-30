# Start by pulling ubuntu image
FROM ubuntu:latest

# update libraries
RUN apt-get clean && apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y apt-transport-https

# non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

#-----------------------
# Installing software
#-----------------------
# install:
#   wget to scrape from web
#   software-properties-common to help manage repos
RUN apt-get update && apt-get install -y \
  wget \
  software-properties-common

# install R from command line; get >= R-3.5
RUN add-apt-repository -y ppa:marutter/rrutter4.0
# RUN add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

# install:
#   curl
#   libcurl, Java (for h20)
#   r and r-dev
#   pandoc (for Rmarkdown conversions)
#   vim (for editing while in container)
#   nginx (for static website hosting)
#   ffmpeg (for animating figures)
# No need to install pandoc-citeproc -> https://github.com/jgm/pandoc-citeproc?tab=readme-ov-file#pandoc-citeproc
RUN apt-get update && apt-get install -y -qq --no-install-recommends --purge \
  curl \
  libcurl4-openssl-dev \
  openjdk-8-jdk \
  r-base \
  r-base-dev \
  pandoc \
  cmake \
  nginx \
  ffmpeg
#   r-cran-devtools \
  #   vim \

RUN rm /var/www/html/index.nginx-debian.html

# install R libraries needed for analysis
RUN Rscript -e 'install.packages("nloptr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("rmarkdown", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("bookdown", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("seqinr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("SuperLearner", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("quadprog", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("remotes", repos="https://cran.rstudio.com")'
RUN Rscript -e 'remotes::install_github("benkeser/cvma")'
# get ggplot2, dplyr, tidyr, readr, tibble, stringr, forcats (and purrr for free)
RUN Rscript -e 'install.packages("tidyverse", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("cowplot", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("glmnet", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("ranger", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("xgboost", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("sandwich", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("ggseqlogo", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("shiny", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("testthat", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("RCurl", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("bit64", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("jsonlite", repos="https://cran.rstudio.com")'
RUN Rscript -e 'options(timeout=600); install.packages("h2o", type="source", repos="http://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/4/R"); library(h2o)'
RUN Rscript -e 'library(h2o); h2o.init()'

# RUN wget "https://h2o-release.s3.amazonaws.com/h2o/latest_stable_R"
# RUN Rscript -e 'install.packages("h20", type = "source"); library(h2o)'
RUN Rscript -e 'install.packages("gam", repos="https://cran.rstudio.com")'

RUN Rscript -e 'remotes::install_github(repo = "bdwilliamson/vimp"); library(vimp)'
# RUN Rscript -e 'remotes::install_version("vimp", version="2.1.9", repos="https://cran.rstudio.com"); library(vimp)'

ADD slapnap slapnap/
RUN R CMD INSTALL -dc slapnap
# RUN Rscript -e 'setwd("slapnap/"); devtools::check()'
# RUN Rscript -e 'install.packages("slapnap_0.1.0.tar.gz", type="source")'
RUN Rscript -e 'library("slapnap")'

RUN rm -rf slapnap_0.1.0.tar.gz slapnap/

# Create default user and workdir
RUN useradd -ms /bin/bash slapnap
WORKDIR /home/slapnap/

RUN chown -R slapnap:slapnap /home/slapnap/

#---------------------------------------------------------------------
# Permanent options
#---------------------------------------------------------------------
RUN mkdir -p slfits output dat/analysis dat/catnap

# Slapnap data analysis directory
ENV analysis="/home/slapnap/dat/analysis/"

# Slapnap slfits directory
ENV slfits="/home/slapnap/slfits/"

# Slapnap output directory
ENV output="/home/slapnap/output/"

# which antibody to analyze
#   "VRC01" is arbitrarily selected as default
ENV nab="VRC01"

# which outcomes to include in the analysis
#   possible outcomes include "ic50", "ic80",
#   "iip", "sens", "estsens", "multsens" and semicolon-separated
#   combinations of these
#   For a single/multispecific bnAb, enter "sens".
#   For a bnAb combination, enter "estsens" or "multsens".
ENV outcomes="ic50"

# which method to use for predicting combination IC-50 and IC-80
#   possible methods are "additive" and "Bliss-Hill". For "Bliss-Hill",
#   "bliss-hill", "bh", or "BH" may also be entered.
ENV combination_method="additive"

# whether or not to use IC-50 or IC-80 to define binary outcomes
#   possible values are "ic50" or "ic80"
ENV binary_outcomes="ic50"

# which learners are included by default
#  if more than a single algorithm is listed, then super learner is used
#  if a single algorithm is listed, then the boolean `cvtune` variable can be used
#  to determine if default tuning parameters are selected or if a small grid
#  search is performed to select tuning parameters.
#
#  rf = random forest
#  xgboost = eXtreme gradient boosting
#  lasso = elastic net regression
#  h2oboost = gradient boosting using H2O.ai
ENV learners="rf"

# should cv be used to select tuning parameters?
#   if TRUE, then a small grid search is performed to select tuning parameters
#   if FALSE, then the "default" tuning parameters of the respective R packages are used
#   note: if more than one learner, then this option controls whether a single version of each
#    algorithm is included in the super learner, or multiple.
ENV cvtune="FALSE"

# should cv be used to measure performance?
#   if TRUE, then cross-validation is used to validate the performance of the prediction
#     algorithm in predicting the selected outcomes
#   if FALSE, then the learner is trained on each outcome, but nothing more is performed
ENV cvperf="TRUE"

# how many folds should be used for cross-validation?
#   only has an effect if cvtune=TRUE or cvperf=TRUE
#   note that intrinsic importance is based on (nfolds / 2)-fold cross-validation
ENV nfolds="2"

# what group-level importance measures should be computed?
#   possible values are 'marg' (for marginal), 'cond' (for conditional), 'marg;cond' (for both marginal and conditional), or none (input "")
ENV importance_grp=""
# what individual-level importance measures should be computed?
#   possible values are 'marg' (for marginal), 'cond' (for conditional), 'pred' (for ML-specific predictive importance), a semicolon-separated combination of these three, or none (input "")
ENV importance_ind=""
# should individual-level intrinsic importance be measured on a site-wise or residue-wise basis?
#   possible values are 'sitewise' or 'residuewise' (only used if importance_ind contains 'marg' or 'cond')
Env ind_importance_type="sitewise"

# set the name of the saved report
#  if set to "", then will default to report_[_-separated list of nabs]_[date].html
ENV report_name=""

# output to save in addition to the report
#   a semicolon-separated list of items,
#   including all possible combinations of
#   "report" (default, return the report)
#   "learner" (return the fitted R object)
#   "data" (return the analysis dataset)
#   "figures" (return the figures from the report)
#   "vimp" (return the R variable importance objects)
#   if set to "", then will default to returning only the report
ENV return="report"

# option to control sensitivity threshold for defining dichotomous
# endpoints as (estimated/multiple) IC-50 < sens_thresh
ENV sens_thresh="1"

# number of sensitive abs needed to declare a pseudovirus sensitive
ENV multsens_nab="1"

# option to view output on exposed port
ENV view_port="FALSE"

# option to subset to pseudoviruses that have all measured outcomes
ENV same_subset="FALSE"

# option to set minimum variability threshold for binary features
ENV var_thresh="0"

# add an argument to bust the cache, so that data are downloaded
# fresh every build. taken from this SO answer:
# https://stackoverflow.com/questions/35134713/disable-cache-for-specific-run-commands
ARG CACHEBUST=1
RUN echo "$CACHEBUST"

COPY bin/* /bin/

RUN chown -R slapnap:slapnap /home/slapnap/

# USER slapnap

# pull CATNAP data from LANL
RUN wget -O dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
# RUN wget -O dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
COPY virseqs_aa.fasta dat/catnap/virseqs_aa.fasta
RUN wget -O dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# entry point to container runs run_analysis.sh
CMD run_analysis.sh
