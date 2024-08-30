# path.home <- "/Users/rudolph/Documents/projects/slapnap_nonsudo"

# A function that extracts importance measures for each feature
# from a fitted super learner object. We get the best fitting ranger, xgboost,
# h2oboost, and/or LASSO model from the library and extract importance features from those.
#' @param fit_sl the fitted regression function. If opts$learners is a single learner (e.g., xgboost),
#' then this is an object with class equal to the class of that learner. If opts$learners has multiple learners,
#' then this is a SuperLearner object.
#' @param data the dataset
#' @param outcome_name the outcome of interest
#' @param opts the options
#' @param ... other args
#' @importFrom ranger importance
#' @importFrom sandwich vcovHC
#' @export
extract_importance <- function(fit_sl, opts, ...){
    if (length(opts$learners) > 1) {
        # find index of one with highest weight
        biggest_weight_idx <- which.max(fit_sl$coef)
        fit_object <- fit_sl$fitLibrary[[biggest_weight_idx]]$object
    } else if (length(opts$learners) == 1 & opts$cvtune) {
        # find best fit
        best_fit_idx <- which.min(fit_sl$cvRisk)
        fit_object <- fit_sl$fitLibrary[[best_fit_idx]]$object
    } else if (length(opts$learners) == 1 & opts$cvperf) {
        fit_object <- fit_sl$fitLibrary[[1]]$object
    } else {
        fit_object <- fit_sl
    }

    if ("ranger" %in% class(fit_object)) {
        # get importance of best ranger
    	ranger_imp <- importance(fit_object)
    	ranger_imp_ranks <- rank(-ranger_imp)
        imp_dt <- data.frame(algo = "rf", Feature = names(ranger_imp),
                             rank = ranger_imp_ranks,
                             Importance = ranger_imp)
    }
    if (("xgboost" %in% class(fit_object)) | ("xgb.Booster" %in% class(fit_object))) {
        # get importance of best xgboost
    	xgboost_imp_dt_init <- xgb.importance(model = fit_object)
    	colnames(xgboost_imp_dt_init)[1] <- "Feature"
    	xgboost_imp_dt_init$rank <- seq_len(nrow(xgboost_imp_dt_init))
        imp_dt <- data.frame(algo = "xgboost", Feature = xgboost_imp_dt_init$Feature,
                             rank = xgboost_imp_dt_init$rank,
                             Importance = xgboost_imp_dt_init$Feature)
    }
    if (any(grepl("H2O", class(fit_object)))) {
        # get importance of best h2oboost
        h2oboost_imp_dt_init <- fit_object@model$variable_importances
        h2oboost_imp_dt_init$rank <- seq_len(nrow(h2oboost_imp_dt_init))
        imp_dt <- data.frame(algo = "h2oboost", Feature = h2oboost_imp_dt_init$variable, rank = h2oboost_imp_dt_init$rank, Importance = h2oboost_imp_dt_init$relative_importance)
    }
    if ("cv.glmnet" %in% class(fit_object)) {
        # get importance of best glmnet
    	glmnet_coef <- fit_object$glmnet.fit$beta[, which(fit_object$lambda == fit_object$lambda.min)]
    	glmnet_imp_rank <- rank(-abs(glmnet_coef))
        imp_dt <- data.frame(algo = "lasso", Feature = names(glmnet_coef),
                             rank = glmnet_imp_rank,
                             Importance = glmnet_coef)
    }
    if ("numeric" %in% class(fit_object)) {
        imp_dt <- data.frame(algo = "mean", Feature = fit_sl$varNames, rank = mean(1:length(fit_sl$varNames)), Importance = 0)
    }
    row.names(imp_dt) <- NULL
	return(imp_dt[order(imp_dt$rank), ])
}

# Function to compute univariate measures of variable importance
#' @export
get_uni_pvals <- function(outcome_name, binary_outcome = FALSE, dat,
                          which_cols){
	pred_names <- colnames(dat)[which_cols]
	geographic_vars <- paste0(colnames(dat)[grepl("geographic", colnames(dat))], collapse = " + ")
	pvals <- sapply(pred_names, USE.NAMES = FALSE, function(p){
		fit <- glm(as.formula(paste0(outcome_name, "~", p, "+", geographic_vars)),
			           family = ifelse(binary_outcome, "binomial", "gaussian"),
			           data = dat)
		if(!is.na(fit$coefficients[2])){
			coef_cov <- vcovHC(fit, type = "HC0")[2,2] # robust variance
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

#' @export
get_all_importance <- function(outcome_name, binary_outcome,
                               dir_loc = Sys.getenv('slfits'),
                               max_rank_thresh = 1:50,
                               n_importance_measures = 1,
                               dat, which_cols, opts){
 	fit_sl <- readRDS(paste0(dir_loc, "fit_", outcome_name, ".rds"))
 	imp_df <- extract_importance(fit_sl, data = dat, outcome_name = outcome_name, opts = opts)
    imp_df$variable <- as.character(imp_df$variable)


 	# get p-value importance
 	uni_pvals <- get_uni_pvals(outcome_name = outcome_name, binary_outcome = binary_outcome, dat = dat, which_cols = which_cols)
 	imp_df$pval_imp <- uni_pvals < 0.05
 	imp_df$pval <- uni_pvals
    pval_idx <- max(which(grepl("pval", names(imp_df))))

 	any_below_max_rank_thresh <- sapply(max_rank_thresh, function(m){
	 	apply(imp_df %>% select(-variable), 1, function(x){
	 		sum(x[grepl("rank", names(x))] <= m, na.rm = TRUE) >= n_importance_measures & x[grepl("pval_imp", names(x))]
	 	})
 	})
 	imp_df <- cbind(imp_df, any_below_max_rank_thresh)
 	colnames(imp_df)[(pval_idx + 1):ncol(imp_df)] <- paste0("any_imp_", max_rank_thresh)
 	# make indicator that it's in top
 	return(imp_df)
}

# make a table
#' @export
get_importance_table <- function(imp_df, max_features = 10){
	imp_df_df <- data.frame(imp_df)
	n_ft <- apply(imp_df_df[,grep("any_imp", colnames(imp_df))], 2, function(x){
		sum(x)
	})
	col_idx <- max(which(n_ft < max_features))
	col_name <- paste0("any_imp_", col_idx)
	feature_tab <- imp_df_df[,c("variable", names(imp_df_df)[grepl("rank", names(imp_df_df))], "pval")][imp_df_df[,col_name],]
	feature_tab$pval <- round(feature_tab$pval, 3)
	feature_tab$pval[feature_tab$pval < 0.001] <- "<0.001"
	return(feature_tab)
}

#' @export
get_importance_resis <- function(imp_df, dat, which_outcome){
	Y <- dat[,which_outcome]
	out <- apply(imp_df, 1, function(x){
		this_x <- dat[,x[1]]
		cor(Y, this_x) >= 0
	})
	return(out)
}

# imp_list is a list of importance tables across different outcomes
#' @export
combine_importance <- function(imp_list,
                               out_names = c("IC50", "IC80", "IIP", "Estimated sens.", "Multiple sens.")){
	full_list <- Reduce("rbind", imp_list)
	imp_ft_more_than_one_outcome <- unique(full_list$variable[duplicated(full_list$variable)])
	# find which outcomes each one is important for
	if (length(imp_ft_more_than_one_outcome) > 0) {
		imp <- NULL
		for (v in imp_ft_more_than_one_outcome) {
			this_imp <- NULL
			for (o in seq_along(out_names)) {
				this_imp <- c(this_imp, v %in% imp_list[[o]]$variable)
			}
			imp <- c(imp, paste0(out_names[this_imp], collapse = ", "))
		}
		out <- data.frame(Variable = imp_ft_more_than_one_outcome,
	                  Outcomes = imp)
	} else {
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
