# A function that extracts importance measures for each feature
# from a fitted super learner object. We get the best fitting ranger, xgboost,
# and LASSO model from the library and extract importance features from those.
extract_importance <- function(fit_sl, data, outcome_name, ...){
	# find name of best fitting ranger model
	ranger_fit_idx <- grep("SL.ranger", fit_sl$libraryNames)
	xgboost_fit_idx <- grep("SL.xgboost", fit_sl$libraryNames)
	glmnet_fit_idx <- grep("SL.glmnet", fit_sl$libraryNames)

	# best of each algorithm
	min_ranger_idx <- which.min(fit_sl$cvRisk[ranger_fit_idx])
	best_ranger <- names(fit_sl$cvRisk)[ranger_fit_idx][min_ranger_idx]
	min_xgboost_idx <- which.min(fit_sl$cvRisk[xgboost_fit_idx])
	best_xgboost <- names(fit_sl$cvRisk)[xgboost_fit_idx][min_xgboost_idx]
	min_glmnet_idx <- which.min(fit_sl$cvRisk[glmnet_fit_idx])
	best_glmnet <- names(fit_sl$cvRisk)[glmnet_fit_idx][min_glmnet_idx]

	# idx of columns of Z matrix corresponding to best fit of each algo
	col_Z_idx <- c(ranger_fit_idx[min_ranger_idx], xgboost_fit_idx[min_xgboost_idx], 
	               glmnet_fit_idx[min_glmnet_idx])

	# check that best one predicts better than chance
	if(grepl("dichotomous", outcome_name)){
		# compute CV-AUC for each algorithm
		V <- length(fit_sl$validRows)
		cv_auc <- NULL
		for(j in col_Z_idx){
			split_Z <- vector(mode = "list", length = V)		
			split_Y <- vector(mode = "list", length = V)
			for(v in seq_len(V)){
				split_Z[[v]] <- fit_sl$Z[fit_sl$validRows[[v]], j]
				split_Y[[v]] <- data[fit_sl$validRows[[v]], outcome_name]
			}
			this_cv_auc <- cvAUC::cvAUC(predictions = split_Z, labels = split_Y)$cvAUC
			cv_auc <- c(cv_auc, this_cv_auc)
		}
		include_ranger <- cv_auc[1] > 0.5
		include_xgboost <- cv_auc[2] > 0.5
		include_glmnet <- cv_auc[3] > 0.5
	}else{
		# compute CV-R^2 for each algorithm
		include_ranger <- 1 - fit_sl$cvRisk[ranger_fit_idx][min_ranger_idx]/var(data[,outcome_name]) > 0
		include_xgboost <- 1 - fit_sl$cvRisk[xgboost_fit_idx][min_xgboost_idx]/var(data[,outcome_name]) > 0
		include_glmnet <- 1 - fit_sl$cvRisk[glmnet_fit_idx][min_glmnet_idx]/var(data[,outcome_name]) > 0
	}

	# get importance of best ranger
	
		best_ranger_fit <- fit_sl$fitLibrary[[best_ranger]]$object
		ranger_imp <- importance(best_ranger_fit)
		ranger_imp_ranks <- rank(-ranger_imp)
	if(!include_ranger){
		# replace ranks with Inf
		ranger_imp_ranks <- rep(Inf, length(ranger_imp_ranks))
	}

	# get importance of best xgboost
	best_xgboost_fit <- fit_sl$fitLibrary[[best_xgboost]]$object
	xgboost_imp_dt <- xgb.importance(model = best_xgboost_fit)
	colnames(xgboost_imp_dt)[1] <- "variable"
	xgboost_imp_dt$rank <- seq_len(nrow(xgboost_imp_dt))
	if(!include_xgboost){
		xgboost_imp_dt$rank <- rep(Inf, length(xgboost_imp_dt$rank))
	}

	# get importance of best glmnet
	best_glmnet_fit <- fit_sl$fitLibrary[[best_glmnet]]$object
	glmnet_coef <- best_glmnet_fit$glmnet.fit$beta[, which(best_glmnet_fit$lambda == best_glmnet_fit$lambda.min)]
	glmnet_imp_rank <- rank(-abs(glmnet_coef))
	if(!include_glmnet){
		glmnet_imp_rank <- rep(Inf, length(glmnet_imp_rank))
	}

	# make a data frame to hold results
	imp_df_init <- data.frame(variable = names(ranger_imp),
	                     rf_rank = ranger_imp_ranks,
	                     glmnet_rank = glmnet_imp_rank)
	# imp_df_init$order <- seq_len(nrow(imp_df_init))
	imp_df <- merge(x = xgboost_imp_dt[,c("variable", "rank")], y = imp_df_init, all.y = TRUE, all.x = TRUE)
	# imp_df <- imp_df[order(imp_df$order), ]
	colnames(imp_df) <- c("variable", "xgboost_rank", "rf_rank", "glmnet_rank")
	return(imp_df)
}

# Function to compute univariate measures of variable importance
get_uni_pvals <- function(outcome_name, binary_outcome = FALSE, dat,
                          which_cols){
	pred_names <- colnames(dat)[which_cols]
	geographic_vars <- paste0(colnames(dat)[grepl("geographic", colnames(dat))], collapse = " + ")
	pvals <- sapply(pred_names, USE.NAMES = FALSE, function(p){
		fit <- glm(as.formula(paste0(outcome_name, "~", p, "+", geographic_vars)),
			           family = ifelse(binary_outcome, "binomial", "gaussian"),
			           data = dat)
		if(!is.na(fit$coefficients[2])){
			coef_cov <- sandwich::vcovHC(fit, type = "HC0")[2,2] # robust variance
			coef_pred <- fit$coef[2]
			wald_stat <- abs(coef_pred / sqrt(coef_cov))
			pval <- 2*pnorm(-wald_stat)
		}else{
			pval <- 1
		}
		return(pval)
	})
	adj_pvals <- p.adjust(pvals, method = "holm")
	return(adj_pvals)
}

get_all_importance <- function(outcome_name, binary_outcome,
                               dir_loc = "/home/slfits/",
                               max_rank_thresh = 1:50,
                               n_importance_measures = 1,
                               dat, which_cols, reduce_covs){
 	fit_sl <- readRDS(paste0(dir_loc, "fit_", outcome_name, ".rds"))
 	imp_df <- extract_importance(fit_sl, data = dat, outcome_name = outcome_name)
    imp_df$variable <- as.character(imp_df$variable)


 	# get p-value importance
 	uni_pvals <- get_uni_pvals(outcome_name = outcome_name, binary_outcome = binary_outcome,
 	                           dat = dat, which_cols = which_cols)
    ## if reduce_covs, need to add on blank rows
    if (reduce_covs) {
        tmp_extra_covs <- data.frame(variable = names(dat)[which_cols], xgboost_rank = NA, rf_rank = NA, glmnet_rank = NA, stringsAsFactors = FALSE)
        tmp_joined_df <- suppressWarnings(dplyr::anti_join(tmp_extra_covs, imp_df, by = "variable"))
        tmp_df <- rbind.data.frame(imp_df, tmp_joined_df)
        imp_df <- tmp_df
    }
 	# init_ranks <- rank(uni_pvals)
 	imp_df$pval_imp <- uni_pvals < 0.05
 	imp_df$pval <- uni_pvals

 	any_below_max_rank_thresh <- sapply(max_rank_thresh, function(m){
	 	apply(imp_df[,!"variable", with = FALSE], 1, function(x){
	 		sum(x[1:3] <= m, na.rm = TRUE) >= n_importance_measures & x[4]
	 	})
 	})
 	imp_df <- cbind(imp_df, any_below_max_rank_thresh)
 	colnames(imp_df)[7:ncol(imp_df)] <- paste0("any_imp_", max_rank_thresh)
 	# make indicator that it's in top
 	return(imp_df)
}

# make a table
get_importance_table <- function(imp_df, max_features = 10){
	imp_df_df <- data.frame(imp_df)
	n_ft <- apply(imp_df_df[,grep("any_imp", colnames(imp_df))], 2, function(x){
		sum(x)
	})
	col_idx <- max(which(n_ft < max_features))
	col_name <- paste0("any_imp_", col_idx)
	feature_tab <- imp_df_df[,c(1:4,6)][imp_df_df[,col_name],]
	feature_tab$pval <- round(feature_tab$pval, 3)
	feature_tab$pval[feature_tab$pval < 0.001] <- "<0.001"
	return(feature_tab)
}


get_importance_resis <- function(imp_df, dat, which_outcome){
	Y <- dat[,which_outcome]
	out <- apply(imp_df, 1, function(x){
		this_x <- dat[,x[1]]
		cor(Y, this_x) >= 0
	})
	return(out)
}

# imp_list is a list of importance tables across different outcomes
combine_importance <- function(imp_list,
                               out_names = c("IC50", "IC80", "IIP", "Estimated sens.", "Multiple sens.")){
	full_list <- Reduce("rbind", imp_list)
	imp_ft_more_than_one_outcome <- unique(full_list$variable[duplicated(full_list$variable)])
	# find which outcomes each one is important for
	if(length(imp_ft_more_than_one_outcome) > 0){
		imp <- NULL
		for(v in imp_ft_more_than_one_outcome){
			this_imp <- NULL
			for(o in seq_along(out_names)){
				this_imp <- c(this_imp, v %in% imp_list[[o]]$variable)
			}
			imp <- c(imp, paste0(out_names[this_imp], collapse = ", "))
		}
		out <- data.frame(Variable = imp_ft_more_than_one_outcome,
	                  Outcomes = imp)
	}else{
		out <- data.frame(Variable = "None", Outcomes = "N/A")
	}
	
	return(out)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get variable importance measures and save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# outcome_names <- c("pc.ic50", "pc.ic80", "iip",
#                    "dichotomous.1", "dichotomous.2")
# library(ranger); library(xgboost); library(glmnet)
# for(o in outcome_names){
# 	imp_df <- get_all_importance(o, binary_outcome = grepl("dichot", o))
# 	save(imp_df, file = paste0("/home/slfits/imp_", o, ".rds"))
# }

# geog_idx <- min(grep("geographic.region.of", colnames(dat))) # geography seems to be first column of relevant data
# pred_names <- colnames(dat)[geog_idx:ncol(dat)]

# library(xgboost); library(ranger); library(glmnet)
# imp_ic50 <- get_all_importance("log10.pc.ic50", binary_outcome = FALSE,
#                                dir_loc = "~/Dropbox/Emory/AMP/slapnap/slfits/",
#                                dat = dat, which_cols = geog_idx:ncol(dat))
