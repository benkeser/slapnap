# A function that extracts importance measures for each feature
# from a fitted super learner object. We get the best fitting ranger, xgboost,
# and LASSO model from the library and extract importance features from those.
extract_importance <- function(fit_sl, ...){
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
	
	# get importance of best ranger
	best_ranger_fit <- fit_sl$fitLibrary[[best_ranger]]$object
	ranger_imp <- importance(best_ranger_fit)	
	ranger_imp_ranks <- rank(-ranger_imp) 

	# get importance of best xgboost
	best_xgboost_fit <- fit_sl$fitLibrary[[best_xgboost]]$object
	xgboost_imp_dt <- xgb.importance(model = best_xgboost_fit)
	colnames(xgboost_imp_dt)[1] <- "variable"
	xgboost_imp_dt$rank <- seq_len(nrow(xgboost_imp_dt))
	xgboost_imp_ranks <- rank(-xgboost_imp) 

	# get importance of best glmnet
	best_glmnet_fit <- fit_sl$fitLibrary[[best_glmnet]]$object
	glmnet_coef <- best_glmnet_fit$glmnet.fit$beta[, which(best_glmnet_fit$lambda == best_glmnet_fit$lambda.min)]
	glmnet_imp_rank <- rank(-abs(glmnet_coef))

	# make a data frame to hold results
	imp_df_init <- data.frame(variable = names(ranger_imp),
	                     rf_rank = ranger_imp_ranks, 
	                     glmnet_rank = glmnet_imp_rank)
	imp_df_init$order <- seq_len(nrow(imp_df_init))
	imp_df <- merge(xgboost_imp_dt[,c("variable", "rank")], imp_df_init)
	imp_df <- imp_df[order(imp_df$order), ]
	return(imp_df[,-5])
}

# Function to compute univariate measures of variable importance
get_uni_pvals <- function(outcome_name, binary_outcome = FALSE){
	analysis_data_name <- list.files("/home/dat/analysis")
	dat <- read.csv(paste0("/home/dat/analysis/", analysis_data_name), header = TRUE)
	dat <- dat[complete.cases(dat),]
	non_pred_names <- c("pc.ic50", "pc.ic80", "iip",
	                    "dichotomous.1", "dichotomous.2",
	                    "seq.id.lanl","seq.id.catnap")
	pred_names <- colnames(dat)[!(colnames(dat) %in% non_pred_names)]

	pvals <- sapply(pred_names, USE.NAMES = FALSE, function(p){
		fit <- glm(as.formula(paste0(outcome_name, "~", p)),
		           family = ifelse(binary_outcome, "binomial", "gaussian"),
		           data = dat)
		coef_cov <- sandwich::vcovHC(fit, type = "HC0")[2,2] # robust variance
		coef_pred <- fit$coef[2]
		wald_stat <- abs(coef_pred / sqrt(coef_cov))
		pval <- 2*pnorm(-wald_stat)
		return(pval)
	})
	return(pvals)
}

get_all_importance <- function(outcome_name, binary_outcome){
 	fit_sl <- get(load(paste0("/home/slfits/fitted_", outcome_name, ".rds")))
 	imp_df <- extract_importance(fit_sl)
 	uni_pvals <- get_uni_pvals(outcome_name = outcome_name, binary_outcome = binary_outcome)
 	imp_df$pval_imp <- rank(uni_pvals)
 	return(imp_df)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get variable importance measures and save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outcome_names <- c("pc.ic50", "pc.ic80", "iip",
                   "dichotomous.1", "dichotomous.2")
for(o in outcome_names){
	imp_df <- get_all_importance(o, binary_outcome = grepl("dichot", o))
	save(imp_df, file = paste0("/home/slfits/imp_", o, ".rds"))
}