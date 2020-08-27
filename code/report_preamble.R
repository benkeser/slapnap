## ----------------------------------------------------------------------------
## Knitr options, packages, files
## ----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

options(stringsAsFactors = FALSE)

#---------------------
# Some set up
#---------------------
slfits_dir <- "/home/slfits/"
data_dir <- "/home/dat/"

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

# get biological importance
vimp_threshold <- 0.05

# set number of predictive importance features to show
n_imp_ft <- 20

#---------------------
# Permanent options
#---------------------
# read in options
opts <- get_global_options()
filename <- Sys.getenv("nab_str")

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
    analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
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
outcome_names <- get_outcome_names(opts)
all_outcome_names <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")

# get variable groups
all_var_groups <- get_variable_groups(dat, pred_names)
all_geog_vars <- pred_names[grepl("geog", pred_names)]
num_covs <- length(pred_names) - length(all_geog_vars)
var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]

# set number of CV folds
V <- as.numeric(opts$nfolds)

# check the outcomes to see if we can run them or not
run_sl_vimp_bools <- check_outcomes(dat, outcome_names, V)
run_sl_vimp_bools2 <- lapply(check_outcomes(dat, outcome_names, V), function(x){
    x[c("ic50", "ic80", "iip", "sens1", "sens2") %in% opts$outcomes]
})

## get biological importance
if (!(all(opts$importance_grp == "")) | ("marg" %in% opts$importance_ind) | ("cond" %in% opts$importance_ind)) {
    source(paste0(code_dir, "biological_importance.R"), local = TRUE)
} else {
    ran_vimp_dichot1 <- FALSE
    ran_vimp_dichot2 <- FALSE
    any_signif_all <- FALSE
    log10.pc.ic50_any_signif <- log10.pc.ic80_any_signif <- iip_any_signif <- dichotomous.1_any_signif <- dichotomous.2_any_signif <- FALSE
    log10.pc.ic50_any_signif_cv <- log10.pc.ic80_any_signif_cv <- iip_any_signif_cv <- dichotomous.1_any_signif_cv <- dichotomous.2_any_signif_cv <- FALSE
}

#-------------------------------------------------------
# other variables and objects needed throughout report
#-------------------------------------------------------

# boolean of single nab
one_nab <- length(opts$nab) == 1
est_fillin <- ifelse(one_nab, "", "estimated ")
any_cont <- any(c("ic50", "ic80") %in% opts$outcomes)
any_dich <- any(c("sens1", "sens2") %in% opts$outcomes)
# read in number of observations
nprevious <- readRDS(paste0(slfits_dir, "nprevious.rds"))
ncomplete_features <- readRDS(paste0(slfits_dir,"ncomplete_features.rds"))
ncomplete_ic50 <- readRDS(paste0(slfits_dir,"ncomplete_ic50.rds"))
ncomplete_ic80 <- readRDS(paste0(slfits_dir,"ncomplete_ic80.rds"))
ncomplete_ic5080 <- readRDS(paste0(slfits_dir,"ncomplete_ic5080.rds"))

# read in data
analysis_data_names <- list.files("/home/dat/analysis")
analysis_data_name <- get_analysis_dataset_name(analysis_data_names, opts)
dat <- read.csv(paste0(data_dir, "analysis/", analysis_data_name), header = TRUE)

# first covariate column
min_cov_col_idx <- min(grep("geographic", colnames(dat)))
ncol_data_final <- ncol(dat)
# number with complete sequence/geog information
complete_features_idx <- complete.cases(dat[ , min_cov_col_idx:ncol_data_final])

# data with complete features
dat_comp_ft <- dat[complete_features_idx, ]

# iip column
iip_col_idx <- which(colnames(dat_comp_ft) == "iip")
# data used for ic50-related analyses
if(opts$same_subset){
  dat_ic50 <- dat_comp_ft[complete.cases(dat_comp_ft[,-iip_col_idx]), ]
}else{
  dat_ic50 <- dat_comp_ft[!is.na(dat_comp_ft$pc.ic50), ]
}
# data used for ic80-related analyses
if(!opts$same_subset | !("ic80" %in% opts$outcomes & length(opts$outcomes) > 1)){
  dat_ic80 <- dat_comp_ft[!is.na(dat_comp_ft$pc.ic80), ]
}else{
  dat_ic80 <- dat_comp_ft[complete.cases(dat_comp_ft[,-iip_col_idx]), ]
}

dat_ic5080 <- dat_comp_ft[complete.cases(dat_comp_ft[,-iip_col_idx]), ]

# compute number of sequences with same ic50 and ic80,
# these guys do not have IIP defined
iip_undef <- sum(dat_ic5080$log10.pc.ic80 == dat_ic5080$log10.pc.ic50)

# make nice outcome names
all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
all_ncomplete <- c(ncomplete_ic50, ncomplete_ic80, ncomplete_ic5080, ncomplete_ic50, ncomplete_ic50)
names(all_ncomplete) <- all_outcomes
all_labels <- c("IC$_{50}$", "IC$_{80}$", "IIP", ifelse(one_nab, "Sensitivity", "Estimated sensitivity"), "Multiple sensitivity")
nice_outcomes <- opts$outcomes
for(i in seq_along(all_outcomes)){
    nice_outcomes <- gsub(all_outcomes[i], all_labels[i], nice_outcomes)
}
outcome_names <- get_outcome_names(opts)
ncompletes <- all_ncomplete[names(all_ncomplete) %in% opts$outcomes]
run_sl_vimp_bools <- check_outcomes(dat, outcome_names, V)
# now format continuous outcomes table
cont_idx <- which(opts$outcomes %in% c("ic50", "ic80", "iip"))
bin_idx <- which(opts$outcomes %in% c("sens1", "sens2"))
cont_nms <- nice_outcomes[cont_idx]
bin_nms <- nice_outcomes[bin_idx]
# postfix for naming plots
postfix <- paste0(filename, "_", format(as.Date(Sys.getenv('current_date'), "%d%b%Y"), "%d%b%Y"))

# for plotting IC50, IC80
if (length(opts$nab) == 1) {
    ic50_lab <- bquote(IC[50])
    ic50_loglab <- bquote(log[10]~"(IC"[50]*")")
    ic80_lab <- bquote(IC[80])
    ic80_loglab <- bquote(log[10]~"(IC"[80]*")")
} else {
    ic50_lab <- bquote("Estimated"~IC[50])
    ic50_loglab <- bquote(log[10]~"(Estimated IC"[50]*")")
    ic80_lab <- bquote("Estimated"~IC[80])
    ic80_loglab <- bquote(log[10]~"(Estimated IC"[80]*")")
}

# for if we ran sl for dichotomous endpoints
ran_sl_dichot1 <- run_sl_vimp_bools2$run_sl[grepl("dichotomous.1", names(run_sl_vimp_bools2$run_sl))]
ran_sl_dichot2 <- run_sl_vimp_bools2$run_sl[grepl("dichotomous.2", names(run_sl_vimp_bools2$run_sl))]
ran_sl_dichot1 <- ifelse(length(ran_sl_dichot1) == 0, FALSE, ran_sl_dichot1)
ran_sl_dichot2 <- ifelse(length(ran_sl_dichot2) == 0, FALSE, ran_sl_dichot2)
