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
RUN apt-get update

# make sure we have wget
RUN apt-get install -y wget

# install R from command line; get >= R-3.5
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:marutter/rrutter3.5
RUN apt-get install -y r-base && apt-get install -y r-base-dev

# put vim on for ease of editing docs inside container
RUN apt-get install -y vim

# install pandoc (for Rmarkdown conversions)
RUN apt-get install -y pandoc

# install R libraries needed for analysis
RUN Rscript -e 'install.packages("nloptr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("rmarkdown", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("bookdown", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("seqinr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("SuperLearner", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("quadprog", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("tidyr", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("cowplot", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("glmnet", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("ranger", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("xgboost", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("sandwich", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("ggseqlogo", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("forcats", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("tibble", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("shiny", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("testthat", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("remotes", repos="https://cran.rstudio.com")'
RUN Rscript -e 'remotes::install_github("bdwilliamson/vimp@v2.0.1")'

# install nginx for static website hosting
RUN apt-get install -y nginx
RUN rm /var/www/html/index.nginx-debian.html

# make directories
# lib contains R source files
# dat contains data
# dat/catnap contains original catnap data
# dat/analysis contains analysis data
RUN mkdir /home/dat /home/dat/catnap /home/dat/analysis /home/out
RUN mkdir /home/slfits /home/output

# install ffmpeg for animating figures
RUN apt-get update
RUN apt-get install -y ffmpeg

# copy R scripts to do do data pull, check options, run analysis, and return requested objects (and make executable)
COPY code/multi_ab_v5.Rlib /home/lib/multi_ab_v5.Rlib
COPY code/merge_proc_v4.R /home/lib/merge_proc_v4.R
COPY code/variable_groups.R /home/lib/variable_groups.R
COPY code/utils.R /home/lib/utils.R
COPY code/run_super_learners.R /home/lib/run_super_learners.R
COPY code/get_vimp.R /home/lib/get_vimp.R
COPY code/super_learner_libraries.R /home/lib/super_learner_libraries.R
COPY code/plotting_functions.R /home/lib/plotting_functions.R
COPY code/check_opts.R /home/lib/check_opts.R
COPY code/check_opts_functions.R /home/lib/check_opts_functions.R
COPY code/return_requested_objects.R /home/lib/return_requested_objects.R
COPY code/ml_var_importance_measures.R /home/lib/ml_var_importance_measures.R
COPY code/plot_one_vimp.R /home/lib/plot_one_vimp.R
COPY code/outcome_dist_plot.R /home/lib/outcome_dist_plot.R
COPY code/pred_importance.R /home/lib/pred_importance.R
COPY code/var_import_plot.R /home/lib/var_import_plot.R
COPY code/vimp_executive_summary_table.R /home/lib/vimp_executive_summary_table.R
COPY code/biological_importance.R /home/lib/biological_importance.R

RUN chmod +x /home/lib/merge_proc_v4.R /home/lib/run_super_learners.R /home/lib/get_vimp.R
RUN chmod +x /home/lib/check_opts.R /home/lib/return_requested_objects.R

# copy report Rmd
COPY code/new_report.Rmd /home/lib/new_report.Rmd
COPY code/run_analysis.sh /home/lib/run_analysis.sh
COPY code/render_report.R /home/lib/render_report.R
COPY code/report_preamble.R /home/lib/report_preamble.R
RUN chmod +x /home/lib/run_analysis.sh /home/lib/render_report.R /home/lib/report_preamble.R


#---------------------------------------------------------------------
# Permanent options
#---------------------------------------------------------------------
# which antibody to analyze
#   "VRC01" is arbitrarily selected as default
ENV nab="VRC01"

# which outcomes to include in the analysis
#   possible outcomes include "ic50", "ic80",
#   "iip", "sens", "estsens", "multsens" and semicolon-separated
#   combinations of these
#   For a single/multispecific bnAb, enter "sens".
#   For a bnAb combination, enter "estsens" or "multsens".
ENV outcomes="ic50;sens"

# which learners are included by default
#  if more than a single algorithm is listed, then super learner is used
#  if a single algorithm is listed, then the boolean `cvtune` variable can be used
#  to determine if default tuning parameters are selected or if a small grid
#  search is performed to select tuning parameters.
#
#  rf = random forest
#  xgboost = eXtreme gradient boosting
#  lasso = elastic net regression
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
ENV nfolds="2"

# what group-level importance measures should be computed?
#   possible values are 'marg' (for marginal), 'cond' (for conditional), 'marg;cond' (for both marginal and conditional), or none (input "")
ENV importance_grp=""
# what individual-level importance measures should be computed?
#   possible values are marg (for marginal), cond (for conditional), pred (for ML-specific predictive importance), a semicolon-separated combination of these three, or none (input "")
ENV importance_ind=""

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
# endpoints as (estimated/multiple) IC-50 > sens_thresh
ENV sens_thresh="1"

# option to control multiple
ENV multsens_nab="2"

# option to view output on exposed port
ENV view_port="FALSE"

# add an argument to bust the cache, so that data are downloaded
# fresh every build. taken from this SO answer:
# https://stackoverflow.com/questions/35134713/disable-cache-for-specific-run-commands
ARG CACHEBUST=1
RUN echo "$CACHEBUST"

# pull CATNAP data from LANL
RUN wget -O /home/dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O /home/dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
RUN wget -O /home/dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
RUN wget -O /home/dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# entry point to container runs run_analysis.sh
CMD /home/lib/run_analysis.sh
