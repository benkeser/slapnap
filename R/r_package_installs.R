#!/usr/bin/env Rscript

pkg <- c(
  "rmarkdown",      
  "seqinr",
  "bookdown",
  "glmnet",
  "ranger",
  "SuperLearner",
  "nloptr",
  "quadprog",
  "ggplot2"
)

for(p in pkg){
	install.packages(p)	
}
