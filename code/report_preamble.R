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

# attempt to read in environment variables
reduce_covs <- Sys.getenv("reduce_covs") == "TRUE"
reduce_outcomes <- Sys.getenv("reduce_outcomes") == "TRUE"
reduce_library <- Sys.getenv("reduce_library") == "TRUE"
reduce_groups <- Sys.getenv("reduce_groups") == "TRUE"


source("/home/lib/plotting_functions.R")
# get NAbs names from ENV variable
antibody_string <- Sys.getenv("Nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]
n_ab <- length(antibodies)

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
    if (cap != "") 
        cap = sprintf("<p>%s</p>", cap)
    sprintf("<video %s><source src=\"%s\" />%s</video>", trimws(opts), 
        paste0(opts_knit$get("base.url"), mov.fname), cap)
}

path.home <- "/home/slfits"
path.home <- "/home/dat/analysis/"

geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
pred_names <- colnames(dat)[geog_idx:ncol(dat)]

source("/home/lib/ml_var_importance_measures.R")
source("/home/lib/var_import_plot.R")



# get importance measures
imp_ic50 <- get_all_importance("log10.pc.ic50", binary_outcome = FALSE,
                               dir_loc = "/home/slfits/",
                               dat = dat, which_cols = geog_idx:ncol(dat))
imp_ic80 <- get_all_importance("log10.pc.ic80", binary_outcome = FALSE,
                               dir_loc = "/home/slfits/",
                               dat = dat, which_cols = geog_idx:ncol(dat))
imp_iip <- get_all_importance("iip", binary_outcome = FALSE,
                               dir_loc = "/home/slfits/",
                               dat = dat, which_cols = geog_idx:ncol(dat))
imp_dichot1 <- get_all_importance("dichotomous.1", binary_outcome = TRUE,
                               dir_loc = "/home/slfits/",
                               dat = dat, which_cols = geog_idx:ncol(dat))
imp_dichot2 <- get_all_importance("dichotomous.2", binary_outcome = TRUE,
                               dir_loc = "/home/slfits/",
                               dat = dat, which_cols = geog_idx:ncol(dat))


