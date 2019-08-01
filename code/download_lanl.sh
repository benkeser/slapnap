#!/usr/bin/tcsh

# STEP -1 establish some filesystem constants
set HOME="~/compbio/projectus/cam117_multiAb"

# STEP 0:  download data from LANL
cd $HOME/dat/catnap
wget -O assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
wget -O viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
wget -O virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
wget -O abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"

