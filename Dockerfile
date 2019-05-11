# Start by pulling r-base image
FROM ubuntu:latest

# update libraries
RUN apt-get update -y

# non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# by default non short run mode
# pass environment variable at run time to turn on short-
# run mode
ENV reduce_covs=FALSE
ENV reduce_outcomes=FALSE
ENV reduce_library=FALSE

# install R from command line
RUN apt-get install -y r-base

# install pandoc (for Rmarkdown conversions)
RUN apt-get install -y pandoc

# make directories
# lib contains R source files
# dat contains data 
# dat/catnap contains original catnap data
# dat/analysis contains analysis data
RUN mkdir /home/lib /home/dat /home/dat/catnap /home/dat/analysis /home/out
RUN mkdir /home/slfits

# install R libraries needed for analysis
COPY R/r_package_installs.R /home/lib/r_package_installs.R
RUN chmod +x /home/lib/r_package_installs.R && /home/lib/r_package_installs.R

# make sure we have wget
RUN apt-get install wget

# pull CATNAP data from LANL
RUN wget -O /home/dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O /home/dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
RUN wget -O /home/dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
RUN wget -O /home/dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# copy R scripts to do data pull and make executable
COPY R/multi_ab_v2.Rlib /home/lib/multi_ab_v2.Rlib
COPY R/merge_proc_v2.R /home/lib/merge_proc_v2.R
COPY R/run_super_learners.R /home/lib/run_super_learners.R
RUN chmod +x /home/lib/merge_proc_v2.R /home/lib/run_super_learners.R

# copy report Rmd 
COPY R/report.Rmd /home/lib/report.Rmd
COPY sh/run_analysis.sh /home/lib/run_analysis.sh
COPY R/render_report.R /home/lib/render_report.R
RUN chmod +x /home/lib/run_analysis.sh /home/lib/render_report.R

# e="VRC07-523-LS;PGT121"
# -e  "Nab=VRC07-523-LS;PGT121;PGDM1400"
# # requirements.R file is where we'll define 
# # all the packages that we need to install. 
# COPY requirements.R /tmp/requirements.R
# # test_rscript.R is the R script that executes the job
# COPY test_rscript.R /tmp/test_rscript.R
# # run_remote.sh is a shell script that runs all R code and 
# # copies results back up to s3
# COPY run_remote.sh /tmp/run_remote.sh

# # make run_remote executable
# RUN chmod +x /tmp/run_remote.sh
# # define a save directory that is passed to run_remote.sh
# ENV SAVE_DIR="/tmp/save_dir/"
# # actually make the save directory
# RUN mkdir "$SAVE_DIR"

# # install required libs on container
# RUN Rscript /tmp/requirements.R

# entry point to container runs run_analysis.sh and passes through 
CMD /home/lib/run_analysis.sh