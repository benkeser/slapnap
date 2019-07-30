#!/usr/bin/env Rscript

# attempt to quiet down the package install outputs

# store a copy of system2
assign("system2.default", base::system2, baseenv())

# create a quiet version of system2
assign("system2.quiet", function(...){
  dots <- list(...)
  dots$stdout <- FALSE
  do.call("system2.default", dots)
}, baseenv())
# overwrite system2 with the quiet version
assignInNamespace("system2", system2.quiet, "base")

grbg <- function(...){
  a <- list(...)
}


# install packages
pkg <- c(
  "rmarkdown",      
  "seqinr",
  "bookdown",
  # "glmnet",
  # "ranger",
  "SuperLearner",
  "nloptr",
  "quadprog",
  "dplyr",
  "tidyr",
  "ggplot2",
  "cowplot"
)

for(p in pkg){
	suppressMessages(install.packages(p))
}
suppressMessages(install.packages("/home/lib/vimp_1.3.0.tar.gz", type = "source", repos = NULL))
