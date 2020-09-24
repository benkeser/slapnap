# # identify our filesystem locations; you will need to set "path.home" to the
# # path of the repository on your local filesystem
# path.home <- "/home"
# path.home <- "~/Dropbox/Emory/AMP/slapnap/slfits"
# path.lib <- paste0(path.home, "/code")
# path.data <- paste0(path.home, "/data")
# path.rdata <- paste0(path.home, "/data")

# # load data
# analysis_data_name <- list.files("/home/dat/analysis")
# analysis_data_name <- "/multiab_catnap_VRC07-523-LS_PGT121_PGDM1400_23Sep2019.csv"
# dat <- read.csv(paste0(path.home, analysis_data_name), header = TRUE)
# dat <- dat[complete.cases(dat),]

make_ml_import_plot <- function(dat, imp_df, max_import = 25, 
                                plot_y_lim = 650 * 0.3, 
                                dist_btw_axis_labs = 55 * 0.5,
                                add_to_pos_for_region_label = 25 * 0.5,
                                dist_btw_main_region_axis = 50 * 0.5,
                                tick_length = 50 * 0.5,
                                dist_between_feature_grps = 15 * 0.5,
                                ...){
	# pdf("~/Dropbox/var_imp_test.pdf", height = 5, width = 10)
	
	# set this importance variable
	imp_df$any_imp <- imp_df[[paste0("any_imp_", max_import)]]
	par(bg="white")

	par(oma=c(0, 0, 0, 0), mar=c(3.1, 2.1, 5.2, 8), 
	    mgp = c(5.2, 0.5, 0),
	    xpd=TRUE)

	# set up plot
	plot(x = 0, y = 0, pch = "", xlab="Env Position(HXB2)", ylab="", ylim=c(0, plot_y_lim), 
	     xlim = c(0, 856),
	     type="n", bty="n", xaxt="n", yaxt="n",
	     main = paste0("Top ", max_import, " most important features"))

	# add top and bottom axes
	axis(side=1, at=c(1, seq(100, 800, 100), 856))
	axis(side=3, at=c(1, 30.5, 511.5, 856), pos=plot_y_lim, labels=F)
	text(15, plot_y_lim + dist_btw_axis_labs, labels="Signal\nPeptide")
	text(271, plot_y_lim + dist_btw_axis_labs, labels="gp120")
	text(684, plot_y_lim + dist_btw_axis_labs, labels="gp41")

	# add axes for regions
	axis(side=3, at=c(131, 157.5, 196), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	axis(side=3, at=c(275, 283), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	axis(side=3, at=c(296, 331), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	axis(side=3, at=c(353, 357), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	axis(side=3, at=c(385, 418), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	axis(side=3, at=c(460, 469), pos=plot_y_lim - dist_btw_main_region_axis, labels=F)
	text(144, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "V1", cex=0.7)
	text(177, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "V2", cex=0.7)
	#text(279, plot_y_lim - dist_btw_axis_labs * 0.7 + dist_btw_axis_labs * 0.7, "Loop D", cex=0.7)
	text(270, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "Loop D", cex=0.7)
	text(314, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "V3", cex=0.7)
	text(355, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "Loop E", cex=0.7)
	text(402, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "V4", cex=0.7)
	text(465, plot_y_lim - dist_btw_main_region_axis + dist_btw_axis_labs * 0.7, "V5", cex=0.7)

	# add ticks for important individual amino acid features
	# get all amino acid numbers
	my_filter <- rep (FALSE, nrow (imp_df))
	position <- NULL

	feature.names <- strsplit (imp_df$variable, fixed=TRUE, split=".")
	for (item in 1:length (feature.names)) {
	  if (feature.names[[item]][1] == "hxb2" & feature.names[[item]][3] != "sequon_actual") {
	    my_filter[item] <- TRUE
	    position <- append (position, feature.names[[item]][2])
	  }
	}
	# remove any lingering letters
	lets <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
	          "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
	          "aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj", "kk",
	          "ll", "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt", "uu", "vv",
	          "ww", "xx", "yy", "zz")
	for(l in lets){
		position <- gsub(l, "", position)
	}

	position <- as.numeric(position)

	# go through and add tick if important
	aa_posn_y <- plot_y_lim - dist_btw_main_region_axis - tick_length - dist_between_feature_grps
	# tick_length <- 50
	# dist_between_feature_grps <- 15
	ct <- 0
	for(i in seq_along(my_filter)){	
		if(my_filter[i]){
			ct <- ct + 1
			if(imp_df$any_imp[i]){
				segments(x0 = position[ct], y0 = aa_posn_y, y1 = aa_posn_y + tick_length, col = "chartreuse4")
			}
		}
	}

	# add ticks for sequons 
	my_filter <- rep (FALSE, nrow (imp_df))
	position <- NULL

	feature.names <- strsplit (imp_df$variable, fixed=TRUE, split=".")
	for (item in 1:length (feature.names)) {
	  if (feature.names[[item]][1] == "hxb2" & feature.names[[item]][3] == "sequon_actual") {
	    my_filter[item] <- TRUE
	    position <- append (position, feature.names[[item]][2])
	  }
	}
	# remove any lingering numerics
	lets <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
	          "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
	for(l in lets){
		position <- gsub(l, "", position)
	}

	position <- as.numeric(position)

	# go through and add tick if important
	aa_posn_sequon <- aa_posn_y - (dist_between_feature_grps + tick_length)
	ct <- 0
	for(i in seq_along(my_filter)){	
		if(my_filter[i]){
			ct <- ct + 1
			if(imp_df$any_imp[i]){
				segments(x0 = position[ct], y0 = aa_posn_sequon, y1 = aa_posn_sequon + tick_length, col = "goldenrod4")
			}
		}
	}

	# LENGTH VARIABLES
	length_var_idx <- grep("length\\.", imp_df$variable)

	# find total number of length variables
	n_sig_length <- sum(imp_df$any_imp[grep("length\\.", imp_df$variable)])
	vert_start <- aa_posn_sequon - dist_between_feature_grps
	vert_stop <- aa_posn_sequon - dist_between_feature_grps - tick_length
	ylocs <- seq(vert_start, vert_stop, length = n_sig_length + 2)

	ct <- 1
	if(imp_df$any_imp[which(imp_df == "length.env")]){
		ct <- ct + 1
		arrows(x0 = 1, x1 = 856, y0 = ylocs[ct], col = "blue4",
		       code = 3, length = 0.05)
	}
	if(imp_df$any_imp[which(imp_df == "length.gp120")]){
		ct <- ct + 1
		arrows(x0 = 30, x1 = 511, y0 = ylocs[ct], col = "blue4",
		       code = 3, length = 0.05)
	}
	if(imp_df$any_imp[which(imp_df == "length.v2")]){
		ct <- ct + 1
		arrows(x0 = 158, x1 = 196, y0 = ylocs[ct], col = "blue4",
		       code = 3, length = 0.05)
	}
	if(imp_df$any_imp[which(imp_df == "length.v3")]){
		ct <- ct + 1
		arrows(x0 = 296, x1 = 331, y0 = ylocs[ct], col = "blue4",
		       code = 3, length = 0.05)
	}
	if(imp_df$any_imp[which(imp_df == "length.v5")]){
		ct <- ct + 1
		arrows(x0 = 460, x1 = 469, y0 = ylocs[ct], col = "blue4",
		       code = 3, length = 0.05)
	}

	# Numer of cysteines
	cysteine_var_idx <- grep("cysteine\\.", imp_df$variable)

	# find total number of cysteine variables
	n_sig_cysteine <- sum(imp_df$any_imp[grep("cysteine\\.", imp_df$variable)])
	vert_start <- aa_posn_sequon - 2*dist_between_feature_grps - tick_length
	vert_stop <- aa_posn_sequon - 2*dist_between_feature_grps - 2 * tick_length
	ylocs <- seq(vert_start, vert_stop, length = n_sig_cysteine + 2)

	barbell <- function(x0, x1, y0, col, pch = 19){
		segments(x0 = x0, x1 = x1, y0 = y0, col = col)
		points(x = c(x0,x1), y = c(y0, y0), pch = pch, col = col, cex = 0.75)
	}

	ct <- 1
	if(imp_df$any_imp[which(imp_df == "num.cysteine.env")]){
		ct <- ct + 1
		barbell(x0 = 1, x1 = 856, y0 = ylocs[ct], col = "darkred")
	}
	if(imp_df$any_imp[which(imp_df == "num.cysteine.gp120")]){
		ct <- ct + 1
		barbell(x0 = 30, x1 = 511, y0 = ylocs[ct], col = "darkred")
	}
	if(imp_df$any_imp[which(imp_df == "num.cysteine.v2")]){
		ct <- ct + 1
		barbell(x0 = 158, x1 = 196, y0 = ylocs[ct], col = "darkred")
	}
	if(imp_df$any_imp[which(imp_df == "num.cysteine.v5")]){
		ct <- ct + 1
		barbell(x0 = 460, x1 = 469, y0 = ylocs[ct], col = "darkred")
	}

	# Numer of sequons
	sequons_var_idx <- grep("num\\.sequons\\.", imp_df$variable)

	# find total number of sequons variables
	n_sig_sequons <- sum(imp_df$any_imp[grep("sequons\\.", imp_df$variable)])
	vert_start <- aa_posn_sequon - 3*dist_between_feature_grps - 2*tick_length
	vert_stop <- aa_posn_sequon - 3*dist_between_feature_grps - 3 * tick_length
	ylocs <- seq(vert_start, vert_stop, length = n_sig_sequons + 2)
	# browser()
	ct <- 1
	if(imp_df$any_imp[which(imp_df == "num.sequons.env")]){
		ct <- ct + 1
		barbell(x0 = 1, x1 = 856, y0 = ylocs[ct], col = "cadetblue4", pch = 17)
	}
	if(imp_df$any_imp[which(imp_df == "num.sequons.gp120")]){
		ct <- ct + 1
		barbell(x0 = 30, x1 = 511, y0 = ylocs[ct], col = "cadetblue4", pch = 17)
	}
	if(imp_df$any_imp[which(imp_df == "num.sequons.v2")]){
		ct <- ct + 1
		barbell(x0 = 158, x1 = 196, y0 = ylocs[ct], col = "cadetblue4", pch = 17)
	}
	if(imp_df$any_imp[which(imp_df == "num.sequons.v5")]){
		ct <- ct + 1
		barbell(x0 = 460, x1 = 469, y0 = ylocs[ct], col = "cadetblue4", pch = 17)
	}

	# labels
	text(x = 857, y = aa_posn_y + 0.5 * tick_length, "AA positions", pos=4, col="chartreuse4", cex=0.8)
	text(x = 857, y = aa_posn_sequon + 0.5 * tick_length, "PNGS", pos=4, col="goldenrod4", cex=0.8)
	text(x = 857, y = aa_posn_sequon - dist_between_feature_grps - tick_length * 0.5, 
	     "Length", pos=4, col="blue4", cex=0.8)
	text(x = 857, y = aa_posn_sequon - 2*dist_between_feature_grps - 1.5*tick_length, 
	     "Num. cysteine", pos=4, col="darkred", cex=0.8)
	text(x = 857, y = aa_posn_sequon - 3*dist_between_feature_grps - 2.5*tick_length, 
	     "Num. PNGS", pos=4, col="cadetblue4", cex=0.8)

	# dev.off()

}