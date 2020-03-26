## ----------------------------------------------------------------------------
## Knitr options, packages, files
## ----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

options(stringsAsFactors = FALSE)

library(grid)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(tidyr)
library(ROCR)
library(knitr)
library(gridExtra)
library(xgboost); library(ranger); library(glmnet)
library(vimp)
# code_dir, slfits_dir, data_dir specified in new_report.Rmd
source(paste0(code_dir, "variable_groups.R"))
source(paste0(code_dir, "super_learner_libraries.R"))
source(paste0(code_dir, "utils.R"))
source(paste0(code_dir, "plotting_functions.R"))
source(paste0(code_dir, "ml_var_importance_measures.R"))
source(paste0(code_dir, "var_import_plot.R"))
source(paste0(code_dir, "vimp_executive_summary_table.R"))
source(paste0(code_dir, "plot_one_vimp.R"))
source(paste0(code_dir, "variable_groups.R"))

#---------------------
# Permanent options
#---------------------
# read in options
if (in_container) {
  opts <- get_global_options()
} else {

}

# --------------------
# Antibodies
# --------------------
# get NAbs names from ENV variable
antibody_string <- Sys.getenv("nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]
n_ab <- length(antibodies)

## ----------------------------------------------------------------------------
## Read in the data, do some set-up
## ----------------------------------------------------------------------------
get_dat <- function(){
	# load data
    analysis_data_names <- list.files("/home/dat/analysis")
    # if more than one analysis dataset, use the most recent one
    if (length(analysis_data_names) > 1) {
        analysis_data_name <- analysis_data_names[length(analysis_data_names)]
    } else {
        analysis_data_name <- analysis_data_names
    }
	dat <- read.csv(paste0(data_dir, "analysis/", analysis_data_name), header = TRUE)

	# check missing values
	n_row_prev <- nrow(dat)
	dat <- dat[complete.cases(dat),]
	n_row_now <- nrow(dat)
	return(list(dat = dat, n_row_prev = n_row_prev,
	            n_row_now = n_row_now))
}
dat_list <- get_dat()
dat <- dat_list$dat; n_row_prev <- dat_list$n_row_prev; n_row_now = dat_list$n_row_now

my_webm <- function (x, options) {
	library(knitr)
	format <- "webm"
    x = c(xfun::sans_ext(x), xfun::file_ext(x))
    fig.num = options$fig.num
    format = sub("^[.]", "", format)
    base = sub(paste0(fig.num, "$"), "", x[1])
    fig.fname = paste0(base, "%d", ".", x[2])
    mov.fname = paste0(sub("-$", "", base), ".", format)
    extra = switch(format, webm = paste("-b:v", knitr:::`%n%`(options$ffmpeg.bitrate,
        "1M"), "-crf 10"), mp4 = "-pix_fmt yuv420p")
    more_extra <- paste0('-vf "movie=', fig.fname, ' [fix]; [in] setsar=sar=1,format=rgba [inf]; [inf][fix] blend=all_mode=overlay:all_opacity=1,format=yuva422p10le [out]"')
    ffmpeg.cmd = paste("ffmpeg", "-y", "-r", 1/options$interval,
        "-i", fig.fname, extra, more_extra, mov.fname)
    if (Sys.which("ffmpeg") == "")
        stop2("Could not find ffmpeg command. You should either change the animation.fun ",
            "hook option or install ffmpeg with libvpx enabled.")
    message("executing: ", ffmpeg.cmd)
    system(ffmpeg.cmd, ignore.stdout = TRUE)

    opts = paste(knitr:::sc_split(options$aniopts), collapse = " ")
    opts = paste(sprintf("width=\"%s\"", options$out.width),
        sprintf("height=\"%s\"", options$out.height), opts)
    cap = knitr:::.img.cap(options, alt = TRUE)
    if (cap != "") cap = sprintf("<p>%s</p>", cap)
    message(cap)
    sprintf("<video %s><source src=\"%s\" />%s</video>", trimws(opts),
        paste0(opts_knit$get("base.url"), mov.fname), cap)
}

path.home <- "/home/slfits"
path.home <- "/home/dat/analysis/"

# get names of predictors
geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]

# get names of outcomes
outcome_names <- c(
    switch("ic50" %in% opts$outcomes, "log10.pc.ic50", NULL),
    switch("ic80" %in% opts$outcomes, "log10.pc.ic80", NULL),
    switch("iip" %in% opts$outcomes, "iip", NULL),
    switch("sens1" %in% opts$outcomes, "dichotomous.1", NULL),
    switch("sens2" %in% opts$outcomes, "dichotomous.2", NULL)
)
all_outcome_names <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]

## get biological importance
if (!(opts$importance_grp == "") | ("marg" %in% opts$importance_ind) | ("cond" %in% opts$importance_ind)) {
    source(paste0(code_dir, "biological_importance.R"), local = TRUE)
}
