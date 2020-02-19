## ----------------------------------------------------------------------------
## Knitr options, packages, files
## ----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

options(stringsAsFactors = FALSE)

library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ROCR)
library(knitr)
library(gridExtra)
library(xgboost); library(ranger); library(glmnet)
library(vimp)
source("/home/lib/variable_groups.R")
source("/home/lib/super_learner_libraries.R")
source("/home/lib/utils.R")
source("/home/lib/plotting_functions.R")
source("/home/lib/ml_var_importance_measures.R")
source("/home/lib/var_import_plot.R")
source("/home/lib/vimp_executive_summary_table.R")
source("/home/lib/plot_one_vimp.R")
source("/home/lib/variable_groups.R")

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()

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
	analysis_data_name <- list.files("/home/dat/analysis")
	dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)

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

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]

## ----------------------------------------------------------------------------
## Individual-level algorithm-specific importance
## ----------------------------------------------------------------------------
# get individual p-value-based importance measures; create tables
col_idx <- geog_idx:ncol(dat)
max_features <- 15
if ("log10.pc.ic50" %in% outcome_names) {
    imp_ic50 <- get_all_importance("log10.pc.ic50", binary_outcome = FALSE,
                                   dir_loc = "/home/slfits/",
                                   dat = dat, which_cols = col_idx, opts = opts)
   ic50_tab <- get_importance_table(imp_ic50, max_features = max_features)
   direction_resis <- get_importance_resis(ic50_tab, dat = dat, which_outcome = "log10.pc.ic50")
   ic50_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("log10.pc.ic80" %in% outcome_names) {
    imp_ic80 <- get_all_importance("log10.pc.ic80", binary_outcome = FALSE,
                                   dir_loc = "/home/slfits/",
                                   dat = dat, which_cols = col_idx, opts = opts)
   ic80_tab <- get_importance_table(imp_ic80, max_features = max_features)
   direction_resis <- get_importance_resis(ic80_tab, dat = dat, which_outcome = "log10.pc.ic80")
   ic80_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("iip" %in% outcome_names) {
    imp_iip <- get_all_importance("iip", binary_outcome = FALSE,
                                   dir_loc = "/home/slfits/",
                                   dat = dat, which_cols = col_idx, opts = opts)
   iip_tab <- get_importance_table(imp_iip, max_features = max_features)
   direction_resis <- get_importance_resis(iip_tab, dat = dat, which_outcome = "iip")
   iip_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("dichotomous.1" %in% outcome_names) {
    imp_dichot1 <- get_all_importance("dichotomous.1", binary_outcome = TRUE,
                                   dir_loc = "/home/slfits/",
                                   dat = dat, which_cols = col_idx, opts = opts)
   dichot1_tab <- get_importance_table(imp_dichot1, max_features = max_features)
   direction_resis <- get_importance_resis(dichot1_tab, dat = dat, which_outcome = "dichotomous.1")
   dichot1_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("dichotomous.2" %in% outcome_names) {
    imp_dichot2 <- get_all_importance("dichotomous.2", binary_outcome = TRUE,
                                   dir_loc = "/home/slfits/",
                                   dat = dat, which_cols = col_idx, opts = opts)
   dichot2_tab <- get_importance_table(imp_dichot2, max_features = max_features)
   direction_resis <- get_importance_resis(dichot2_tab, dat = dat, which_outcome = "dichotomous.2")
   dichot2_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
continuous_imp_lst <- list(
    switch("log10.pc.ic50" %in% outcome_names, ic50_tab, NULL),
    switch("log10.pc.ic80" %in% outcome_names, ic80_tab, NULL),
    switch("iip" %in% outcome_names, iip_tab, NULL)
)
dichot_imp_lst <- list(
    switch("dichotomous.1" %in% outcome_names, dichot1_tab, NULL),
    switch("dichotomous.2" %in% outcome_names, dichot2_tab, NULL)
)
imp_continuous <- combine_importance(continuous_imp_lst, out_names = c("IC50", "IC80", "IIP"))
imp_dichot <- combine_importance(dichot_imp_lst, out_names = c("Estimated Sens.", "Multiple Sens."))
imp_overall <- combine_importance(c(continuous_imp_lst, dichot_imp_lst))

## ----------------------------------------------------------------------------
## Population variable importance
## ----------------------------------------------------------------------------
num_pop_import <- 20  # the number of individual features to display in plots

## plotting things
x_lab_continuous <- expression(paste("Difference in ", R^2, sep = ""))
x_lim_continuous <- c(0, 1)
x_lab_binary <- expression(paste("Difference in ", AUC, sep = ""))
x_lim_binary <- c(0, 1)
## read in importance results for each outcome, create a plot for each
## only return non-cv plots if cv = FALSE
imp_nms <- list(all_var_groups, all_var_groups, var_inds)
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    if (grepl("dichot", this_outcome_name)) {
        this_x_lab <- x_lab_binary
        this_x_lim <- x_lim_binary
    } else {
        this_x_lab <- x_lab_continuous
        this_x_lim <- x_lim_continuous
    }
    ## importance results
    eval(parse(text = paste0(this_outcome_name, "_vimp_lst <- readRDS(file = '/home/slfits/", this_outcome_name, "_vimp.rds')")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst <- readRDS(file = '/home/slfits/", this_outcome_name, "_cv_vimp.rds')")))
    eval(parse(text = paste0(this_outcome_name, "_outer_folds <- readRDS(file = '/home/slfits/", this_outcome_name, "_outer_folds.rds')")))
    eval(parse(text = paste0("num_obs_full <- sum(", this_outcome_name, "_outer_folds == 1)")))
    eval(parse(text = paste0("num_obs_red <- sum(", this_outcome_name, "_outer_folds == 2)")))
    ## make plots
    eval(parse(text = paste0("current_vimp_lst <- ", this_outcome_name, "_vimp_lst")))
    eval(parse(text = paste0("current_cv_vimp_lst <- ", this_outcome_name, "_cv_vimp_lst")))
    vimp_plot_titles <- paste0(vimp_plot_name(this_outcome_name), ": ", names(current_vimp_lst))
    eval(parse(text = paste0(this_outcome_name, "_vimp_plots <- mapply(function(x, y) plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = this_x_lim, cv = FALSE, num_plot = num_pop_import), current_vimp_lst, vimp_plot_titles, SIMPLIFY = FALSE)")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_plots <- mapply(function(x, y) plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = this_x_lim, cv = TRUE, num_plot = num_pop_import), current_cv_vimp_lst, vimp_plot_titles, SIMPLIFY = FALSE)")))
}

## make table for executive summary
vimp_threshold <- 0.05
if (opts$cvperf) {
    vimp_summary_tbl <- make_vimp_executive_summary_table(
        switch("log10.pc.ic50" %in% outcome_names, log10.pc.ic50_cv_vimp_lst, NULL),
        switch("log10.pc.ic80" %in% outcome_names, log10.pc.ic80_cv_vimp_lst, NULL),
        switch("iip" %in% outcome_names, iip_cv_vimp_lst, NULL),
        switch("dichotomous.1" %in% outcome_names, dichotomous.1_cv_vimp_lst, NULL),
        switch("dichotomous.2" %in% outcome_names, dichotomous.2_cv_vimp_lst, NULL),
        threshold = vimp_threshold, outcome_names = outcome_names, cv = TRUE, opts = opts)
} else {
    vimp_summary_tbl <- make_vimp_executive_summary_table(switch("log10.pc.ic50" %in% outcome_names, log10.pc.ic50_vimp_lst, NULL),
    switch("log10.pc.ic80" %in% outcome_names, log10.pc.ic80_vimp_lst, NULL),
    switch("iip" %in% outcome_names, iip_vimp_lst, NULL),
    switch("dichotomous.1" %in% outcome_names, dichotomous.1_vimp_lst, NULL),
    switch("dichotomous.2" %in% outcome_names, dichotomous.2_vimp_lst, NULL),
    threshold = vimp_threshold, outcome_names = outcome_names, cv = FALSE, opts = opts)
}