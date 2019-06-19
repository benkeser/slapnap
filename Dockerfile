# Start by pulling ubuntu image
FROM ubuntu:latest

# update libraries
RUN apt-get update && apt-get upgrade -y

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
RUN apt-get install -y wget

# pull CATNAP data from LANL
RUN wget -O /home/dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
RUN wget -O /home/dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
RUN wget -O /home/dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
RUN wget -O /home/dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

# copy R scripts to do data pull and make executable
COPY R/multi_ab_v2.Rlib /home/lib/multi_ab_v2.Rlib
COPY R/merge_proc_v2.R /home/lib/merge_proc_v2.R
COPY R/variable_groups.R /home/lib/variable_groups.R
COPY R/run_super_learners.R /home/lib/run_super_learners.R
COPY R/get_vimp.R /home/lib/get_vimp.R

# copy R package (only until new version gets to GitHub)
COPY vimp_

RUN chmod +x /home/lib/merge_proc_v2.R /home/lib/run_super_learners.R /home/lib/get_vimp.R


# copy report Rmd 
COPY R/report.Rmd /home/lib/report.Rmd
COPY sh/run_analysis.sh /home/lib/run_analysis.sh
COPY R/render_report.R /home/lib/render_report.R
RUN chmod +x /home/lib/run_analysis.sh /home/lib/render_report.R

# entry point to container runs run_analysis.sh
CMD /home/lib/run_analysis.sh




