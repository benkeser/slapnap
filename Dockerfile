# Start by pulling ubuntu image
FROM ubuntu:latest

# update libraries
RUN apt-get clean && apt-get update && apt-get upgrade -y && apt-get install -y apt-transport-https

# non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

#----------------- 
# Temp options
#----------------- 
# eventually we will get rid of these options
ENV reduce_covs=FALSE
ENV reduce_outcomes=FALSE
ENV reduce_library=FALSE
ENV reduce_groups=FALSE
ENV no_cv=FALSE

#---------------------
# Permanent options
#---------------------
# which antibody to analyze
#   "VRC07-523-LS" is arbitrarily selected as default
ENV nab="VRC07-523-LS"

# which outcomes to include in the analysis
ENV outcomes="ic50;ic80;iip;sens1;sens2"
# which learners are included by default
#  if more than a single algorithm is listed, then super learner is used
#  if a single algorithm is listed, then the boolean `cvtune` variable can be used
#  to determine if default tuning parameters are selected or if a small grid
#  search is performed to select tuning parameters. 
# 
#  rf = random forest
#  xgboost = eXtreme gradient boosting
#  lasso = elastic net regression
ENV learners="rf;xgboost;lasso"

# should cv be used to select tuning parameters?
#   if TRUE, then a small grid search is performed to select tuning parameters
#   if FALSE, then the "default" tuning parameters of the respective R packages are used
#   note: if more than one learner, then this option controls whether a single version of each
#    algorithm is included in the super learner, or multiple.
ENV cvtune="TRUE"

# should cv be used to measure performance?
#   if TRUE, then cross-validation is used to validate the performance of the prediction 
#     algorithm in predicting the selected outcomes
#   if FALSE, then the learner is trained on each outcome, but nothing more is performed
ENV cvperf="TRUE"

# what group-level importance measures should be computed?
#   possible values are marg (for marginal) cond (for conditional) or none (input "")
ENV importance_grp="marg;cond"
# what individual-level importance measures should be computed?
#   possible values are marg (for marginal), cond (for conditional), pred (for ML-specific
#   predictive importance), or none (input "")
ENV importance_ind="marg;cond;pred"


#-----------------------
# Installing software
#-----------------------
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

# make directories
# lib contains R source files
# dat contains data
# dat/catnap contains original catnap data
# dat/analysis contains analysis data
RUN mkdir /home/dat /home/dat/catnap /home/dat/analysis /home/out
RUN mkdir /home/slfits

# install ffmpeg for animating figures
RUN apt-get update
RUN apt-get install -y ffmpeg

# pull CATNAP data from LANL
RUN wget -O /home/dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O /home/dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
RUN wget -O /home/dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
RUN wget -O /home/dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# copy R package (only until new version gets to GitHub)
COPY vimp_2.0.0.tar.gz /home/lib/vimp_2.0.0.tar.gz
RUN Rscript -e 'suppressMessages(install.packages("/home/lib/vimp_2.0.0.tar.gz", type = "source", repos = NULL))'
# copy R scripts to do data pull and make executable
COPY code/multi_ab_v4.Rlib /home/lib/multi_ab_v4.Rlib
COPY code/merge_proc_v4.R /home/lib/merge_proc_v4.R
COPY code/variable_groups.R /home/lib/variable_groups.R
COPY code/utils.R /home/lib/utils.R
COPY code/run_super_learners.R /home/lib/run_super_learners.R
COPY code/get_vimp.R /home/lib/get_vimp.R
COPY code/super_learner_libraries.R /home/lib/super_learner_libraries.R
COPY code/plotting_functions.R /home/lib/plotting_functions.R

RUN chmod +x /home/lib/merge_proc_v4.R /home/lib/run_super_learners.R /home/lib/get_vimp.R

# copy report Rmd
COPY code/report.Rmd /home/lib/report.Rmd
COPY code/run_analysis.sh /home/lib/run_analysis.sh
COPY code/render_report.R /home/lib/render_report.R
RUN chmod +x /home/lib/run_analysis.sh /home/lib/render_report.R

# entry point to container runs run_analysis.sh
CMD /home/lib/run_analysis.sh
