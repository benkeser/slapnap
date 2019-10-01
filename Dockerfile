# Start by pulling ubuntu image
FROM ubuntu:latest

# update libraries
RUN apt-get clean && apt-get update && apt-get upgrade -y 

# non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# by default non short run mode
# pass environment variable at run time to turn on short-
# run mode
ENV reduce_covs=FALSE
ENV reduce_outcomes=FALSE
ENV reduce_library=FALSE
ENV reduce_groups=FALSE

# install R from command line
RUN apt-get install -y r-base

# put vim on for ease of editing docs inside container
RUN apt-get install -y vim 

# install pandoc (for Rmarkdown conversions)
RUN apt-get install -y pandoc

# install R libraries needed for analysis
# copy R package (only until new version gets to GitHub)
COPY vimp_1.3.0.tar.gz /home/lib/vimp_1.3.0.tar.gz
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
RUN Rscript -e 'suppressMessages(install.packages("/home/lib/vimp_1.3.0.tar.gz", type = "source", repos = NULL))'
RUN Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("sandwich", repos="https://cran.rstudio.com")'

# make directories
# lib contains R source files
# dat contains data 
# dat/catnap contains original catnap data
# dat/analysis contains analysis data
RUN mkdir /home/dat /home/dat/catnap /home/dat/analysis /home/out
RUN mkdir /home/slfits

# install ffmpeg for animating figures
RUN apt-get install -y ffmpeg

# make sure we have wget
RUN apt-get install -y wget

# pull CATNAP data from LANL
RUN wget -O /home/dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O /home/dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
RUN wget -O /home/dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
RUN wget -O /home/dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# copy R scripts to do data pull and make executable
COPY code/multi_ab_v3.Rlib /home/lib/multi_ab_v3.Rlib
COPY code/merge_proc_v3.R /home/lib/merge_proc_v3.R
COPY code/variable_groups.R /home/lib/variable_groups.R
COPY code/run_super_learners.R /home/lib/run_super_learners.R
COPY code/get_vimp.R /home/lib/get_vimp.R
COPY code/super_learner_libraries.R /home/lib/super_learner_libraries.R
COPY code/plotting_functions.R /home/lib/plotting_functions.R

RUN chmod +x /home/lib/merge_proc_v3.R /home/lib/run_super_learners.R /home/lib/get_vimp.R

# add option to avoid cv super learner fitting (for debugging)
ENV no_cv=FALSE

# copy report Rmd 
COPY code/report.Rmd /home/lib/report.Rmd
COPY code/run_analysis.sh /home/lib/run_analysis.sh
COPY code/render_report.R /home/lib/render_report.R
RUN chmod +x /home/lib/run_analysis.sh /home/lib/render_report.R

# entry point to container runs run_analysis.sh
CMD /home/lib/run_analysis.sh