## various utility functions
#' @param opts options
#' @param imp_df importance data.frame
#' @param n_ft number of features shown
get_importance_text <- function(opts, imp_df, n_ft = 20){
    algo_with_highest_wt <- imp_df$algo[1]

    # check if super learner
    is_sl <- length(opts$learners > 1)
    # check if tuning parameters varied
    is_tuned <- opts$cvtune

    # if random forest
    if(!is_sl & "rf" %in% opts$learners){
        text_out <- paste0("Specifically, random forest variable permutation-based variable importance measures were computed and the top ", n_ft, " features are shown. The permutation-based importance measures the decrease in predictive accuracy when making out-of-bag predictions and randomly permuting a given feature from its original values.")
        if(is_tuned){
            text_out <- paste0(text_out, " These measures are shown for the choice of tuning parameters with the best model fit, as chosen by cross-validation.")
        }
    }else if(!is_sl & "xgboost" %in% opts$learners){
        text_out <- paste0("Specifically, xgboost gain importance measures were computed and the top ", n_ft, " features are shown. Gain measures the improvement in accuracy brought by a given feature to the tree branches on which it appears. The essential idea is that before adding a split on a given feature to the branch, there may be some observations that were poorly predicted, while after adding an additional split on this feature, and each resultant branch is more accurate. Gain measures this change in accuracy.")
        if(is_tuned){
            text_out <- paste0(text_out, " These measures are shown for the choice of tuning parameters with the best model fit, as chosen by cross-validation.")
        }
    }else if(!is_sl & "lasso" %in% opts$learners){
        text_out <- paste0("Specifically, lasso variable importance is taken to be the magnitude of the coefficient for the model with $lambda$ chosen via cross-validation, and the top ", n_ft, " are shown.")
        if(is_tuned){
            text_out <- paste0(text_out, " These ranks are shown for the choice of alpha that resulted in the best model fit, as chosen by cross-validation.")
        }
        text_out <- paste0(text_out, " Overall, there were ", sum(abs(imp_df$value) > 0), " features that had non-zero coefficient in the final fit.")
    }else{
        text_out <- "Specifically, the algorithm with the largest weight in the super learner ensemble was selected and associated variable importance metrics for this algorithm are shown."
        text_out <- paste0(text_out, " In this case, the highest weight was assigned to a ", algo_with_highest_wt, " algorithm, and thus the variable importance measures presented correspond to ")
        if(algo_with_highest_wt == "rf"){
            text_out <- paste0(text_out, "random forest variable permutation-based variable importance measures were computed and are shown by their rank. The permutation-based importance measures the decrease in predictive accuracy when making out-of-bag predictions and randomly permuting a given feature from its original values.")
        }else if(algo_with_highest_wt == "lasso"){
            text_out <- paste0(text_out, "the magnitude of the coefficient for the model with $lambda$ chosen via cross-validation.")
            text_out <- paste0(text_out, " Overall, there were ", sum(abs(imp_df$value) > 0), " features that had non-zero coefficient in the final lasso fit.")
        }else if(algo_with_highest_wt == "xgboost"){
            text_out <- paste0(text_out, "xgboost gain importance measures were computed and are shown by their rank. Gain measures the improvement in accuracy brought by a given feature to the tree branches on which it appears. The essential idea is that before adding a split on a given feature to the branch, there may be some observations that are poorly predicted, while after adding an additional split on this feature, and each resultant branch is more accurate. Gain measures this change in accuracy.")
        }
    }
    return(text_out)
}

# for a given outcome make a panel histogram of the individual
# nabs and a summary table
get_individual_nab_summaries <- function(outcome = "ic50", opts, dat){
    out_hist <- list()
    out_summary <- list()
    # re-label
    outcome_label <- if(outcome == "ic50"){
        "IC-50"
    }else if(outcome == "ic80"){
        "IC-80"
    }else if(outcome == "iip"){
        "IIP"
    }

    ct <- 0
    for(i in seq_along(opts$nab)){
        ct <- ct + 1
        this_name <- gsub("-", ".", paste0(opts$nab[i], ".ic50.imputed"))
        out_hist[[ct]] <- make_hist_plot(dat, var_name = this_name,
                                          x_lab = paste0(outcome_label, opts$nab[i]),
                                          y_lab = "Density")
        tmp_sum <- summary(dat[, this_name])[1:6] # to ignore NA columns
        tmp_sum <- c(tmp_sum[1:3], 10^mean(log10(dat[, this_name])), tmp_sum[4:6])
        names(tmp_sum)[4] <- "Geom. Mean"
        out_summary[[i]] <- tmp_sum
        ct <- ct+1
        dat[,paste0("log10_",this_name)] <- log10(dat[, this_name])
        out_hist[[ct]] <- make_hist_plot(dat, var_name = paste0("log10_",this_name),
                                          x_lab = bquote(log[10]~"(IC-50 "~.(opts$nab[i])~")"),
                                          y_lab = "")
    }
    return(list(hist = out_hist, summary = out_summary))
}

get_learner_descriptions <- function(opts){

    if(length(opts$learners) == 1){
        learner_label <- if(opts$learners == "rf"){
            "random forest"
        }else if(opts$learners == "xgboost"){
            "extreme gradient boosting"
        }else if(opts$learners == "lasso"){
            "elastic net regression"
        }
        tmp <- paste(learner_label,
                     ifelse(opts$cvtune,
                            "with tuning parameters selected using a limited grid search and cross-validation.",
                            "with tuning parameters set to their 'default' values."))
    }else{
        lib_label <- NULL
        if("rf" %in% opts$learners){
            lib_label <- c(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "random forest", ifelse(opts$cvtune, "s with varied tuning parameters", "")))
        }
        if("xgboost" %in% opts$learners){
            if("rf" %in% opts$learners){
                if(!("lasso" %in% opts$learners)){
                    lib_label <- paste0(lib_label, " and ")
                }else{
                    lib_label <- paste0(lib_label, ", ")
                }
            }
            lib_label <- paste0(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "gradient boosted tree", ifelse(opts$cvtune, "s with varied tuning parameters", "")))
        }
        if("lasso" %in% opts$learners){
            if("rf" %in% opts$learners | "xgboost" %in% opts$learners){
                lib_label <- paste0(lib_label, " and ")
            }
            lib_label <- paste0(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "elastic net regression", ifelse(opts$cvtune, "s with varied tuning parameters", "")))
        }
        tmp <- paste0("a super learner ensemble of ", lib_label, ".")
    }
    return(tmp)
}
# given options and n_row_now (for captions) load appropriate SuperLearner/
# CV.SuperLearner fits and create a list of CV results for all continuous
# valued outcomes ([[1]] of output) and dichotomous outcomes ([[2]] of output)

# each entry in the output list is a kable that should be properly labeled.
get_cv_outcomes_tables <- function(fit_list, opts){
    fit_list <- fit_list_out$out
    V <- fit_list_out$V
    n_row_now <- fit_list_out$n_row_now
    table_list <- lapply(fit_list, summary.myCV.SuperLearner, opts = opts)

    # re-label
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    all_labels <- c("IC-50", "IC-80", "IIP", "Estimated", "Multiple")
    tmp <- opts$outcomes
    for(i in seq_along(all_outcomes)){
        tmp <- gsub(all_outcomes[i], all_labels[i], tmp)
    }

    # now format continuous outcomes table
    cont_idx <- which(opts$outcomes %in% c("ic50", "ic80", "iip"))
    rsq_kab <- NULL
    if(length(cont_idx) > 0){
        list_rows <- sapply(cont_idx, get_est_and_ci, fit_list = table_list, Rsquared = TRUE, simplify = FALSE)
        rsqtab <- Reduce(rbind, lapply(list_rows, unlist, use.names = FALSE))
        if(is.null(dim(rsqtab))) rsqtab <- matrix(rsqtab, nrow = 1)
        row.names(rsqtab) <- tmp[cont_idx]
        rsq_kab <- kable(rsqtab, col.names = c(expression(CV-R^2), "Lower 95% CI", "Upper 95% CI"),
              digits = 3, row.names = TRUE,
              caption = paste0("Estimates of ", V, "-fold cross-validated R-squared for super learner predictions ",
                               "of the three continuous-valued outcomes (n = ", n_row_now,
                               " observations with complete sequence data)."))
    }
    # now format dichotomous outcomes table
    dich_idx <- which(opts$outcomes %in% c("sens1", "sens2"))
    auc_kab <- NULL
    if(length(dich_idx) > 0){
        list_rows <- sapply(dich_idx, get_est_and_ci, fit_list = table_list, Rsquared = FALSE, simplify = FALSE)
        auctab <- Reduce(rbind, lapply(list_rows, unlist, use.names = FALSE))
        if(is.null(dim(auctab))) auctab <- matrix(auctab, nrow = 1)
        row.names(auctab) <- tmp[dich_idx]
        auc_kab <- kable(auctab, col.names = c("CVAUC", "Lower 95% CI", "Upper 95% CI"),
              digits = 3, row.names = TRUE,
              caption = paste0("Estimates of ", V, "-fold cross-validated R-squared for super learner predictions ",
                               "of the three continuous-valued outcomes (n = ", n_row_now,
                               " observations with complete sequence data)."))
    }
    return(list(r2 = rsq_kab, auc = auc_kab))
}

# load cv_fits for given set of opts, needed since the naming convention
# is different if length(opts$learners) == 1 and opts$cvtune == FALSE
load_cv_fits <- function(opts, code_dir){
    if(!opts$cvtune & !opts$cvperf){
        stop("no cross-validated fit for these options")
    }
    out_list <- vector(mode = "list", length = length(opts$outcomes))
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    all_file_labels <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")
    if(length(opts$learners) == 1 & !opts$cvtune){
        ct <- 0
        for(i in seq_along(all_outcomes)){
            if(all_outcomes[i] %in% opts$outcomes){
                ct <- ct + 1
                out_list[[ct]] <- readRDS(paste0(code_dir, "fit_", all_file_labels[i], ".rds"))
                class(out_list[[ct]]) <- c("myCV.SuperLearner", class(out_list[[ct]]))
                names(out_list)[ct] <- all_outcomes[i]
            }
        }
    }else{
        ct <- 0
        for(i in seq_along(all_outcomes)){
            if(all_outcomes[i] %in% opts$outcomes){
                ct <- ct + 1
                out_list[[ct]] <- readRDS(paste0(code_dir, "cvfit_", all_file_labels[i], ".rds"))
                class(out_list[[ct]]) <- c("myCV.SuperLearner", class(out_list[[ct]]))
                names(out_list)[ct] <- all_outcomes[i]
            }
        }
    }
    if("SuperLearner" %in% class(out_list[[1]])){
        V <- length(out_list[[1]]$validRows)
    }else{
        V <- out_list[[1]]$V
    }
    # figure out number
    return(list(out = out_list, V = V, n_row_now = length(out_list[[1]]$Y)))
}


# get descriptions of outcomes
get_outcome_descriptions <- function(opts){
    # first describe IC-50, IC-80 if present
    tmp_text <- NULL
    ic50_pres <- "ic50" %in% opts$outcomes
    ic80_pres <- "ic80" %in% opts$outcomes
    iip_pres <- "iip" %in% opts$outcomes
    sens1_pres <- "sens1" %in% opts$outcomes
    sens2_pres <- "sens2" %in% opts$outcomes
    if(length(opts$nab) > 1){
        if(ic50_pres | ic80_pres | iip_pres | sens1_pres){
            tmp <- paste0("Predicted ",
                          ifelse(ic50_pres | iip_pres | sens1_pres, "IC-50 ", ""),
                          ifelse((ic50_pres & ic80_pres) | iip_pres, "and ", ""),
                          ifelse(ic80_pres | iip_pres, "IC-80 ", ""), collapse = "")
            tmp1_5 <- ifelse((ic50_pres & ic80_pres) | iip_pres, "were ", "was ")
            tmp2 <- paste0("computed based on the additive model of Wagh et al. (2016); ",
                           "for $J$ antibodies, it is computed as \\[ \\mbox{predicted IC} = \\left( \\sum_{j=1}^J \\mbox{IC}_j^{-1} \\right)^{-1} \\ , \\]",
                           " where $\\mbox{IC}_j$ denotes the measured ",
                           paste0(ifelse(ic50_pres, "IC-50 ", ""),
                                  ifelse(ic50_pres & ic80_pres, "or ", ""),
                                  ifelse(ic80_pres, "IC-80 ", ""), collapse = ""),
                           "for antibody $j$.")
            tmp_text <- c(tmp_text, paste0(tmp, tmp1_5, tmp2, collapse = ""))
        }
        if(iip_pres){
            tmp <- paste0("IIP is calculated as ",
                          "\\[ \\frac{10^m}{\\mbox{predicted IC-50}^m + 10^m} \ , \\]",
                          "where $m = \\mbox{log}_{10}(4) / (\\mbox{predicted IC-80} - \\mbox{predicted IC-50})$ ",
                          "and predicted IC-50 and IC-80 are computed as described above.",
                          collapse = "")
            tmp_text <- c(tmp_text, tmp)
        }
        if(sens1_pres){
            tmp_text <- c(tmp_text, "Estimated sensitivity is defined by the binary indicator that predicted IC-50 $> 1$.")
        }
        if(sens2_pres){
            tmp_text <- c(tmp_text, "Multiple sensitivity is defined as the binary indicator of having measured IC50 $> 1$ for at least two antibodies.")
        }
    } else {
        if(iip_pres){
            tmp <- paste0("IIP is calculated as ",
                          "\\[ \\frac{10^m}{\\mbox{IC-50}^m + 10^m} \ , \\]",
                          "where $m = \\mbox{log}_{10}(4) / (\\mbox{IC-80} - \\mbox{IC-50})$.",
                          collapse = "")
            tmp_text <- c(tmp_text, tmp)
        }
        if(sens1_pres | sens2_pres){
            tmp_text <- c(tmp_text, "Estimated sensitivity is defined by the binary indicator that IC-50 $> 1$.")
        }
        if(sens2_pres){
            tmp_text <- c(tmp_text, "Since only one antibody was specified for this analysis, multiple sensitivity is the same as estimated sensitivity.")
        }
    }
    return(paste0(tmp_text, collapse = " "))
}


# get a comma separated list of outcomes for report
get_comma_sep_outcomes <- function(opts){
    tmp <- paste0(opts$outcomes, collapse = ", ")
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    all_labels <- c("IC-50", "IC-80", "IIP", "estimated sensitivity", "multiple sensitivity")
    for(i in seq_along(all_outcomes)){
        tmp <- gsub(all_outcomes[i], all_labels[i], tmp)
    }
    tmp <- paste0(tmp, ".")
    return(tmp)
}
## read in a system variable
## use boolean for boolean options
## use !boolean for semi-colon-separated list options
get_sys_var <- function(option = "nab", boolean = FALSE){
    read_string <- Sys.getenv(option)
    if(boolean){
        out <- read_string == "TRUE"
    }else{
        out <- strsplit(read_string, split = ";")[[1]]
    }
    return(out)
}

## read in permanent options
get_global_options <- function(options = c("nab","outcomes", "learners", "cvtune", "cvperf",
                                           "importance_grp", "importance_ind"),
                               options_boolean = c(FALSE, FALSE, FALSE, TRUE,
                                                   TRUE, FALSE, FALSE)){
    out <- mapply(option = options, boolean = options_boolean,
                  FUN = get_sys_var, SIMPLIFY = FALSE)
    return(out)
}

## make outer folds for VIM hypothesis test (based on sample splitting)
make_folds <- function(y, V, stratified = TRUE) {
  if (stratified) {
    y_1 <- y == 1
    y_0 <- y == 0
    folds_1 <- rep(seq_len(V), length = sum(y_1))
    folds_1 <- sample(folds_1)
    folds_0 <- rep(seq_len(V), length = sum(y_0))
    folds_0 <- sample(folds_0)
    folds <- vector("numeric", length(y))
    folds[y_1] <- folds_1
    folds[y_0] <- folds_0
  } else {
    folds <- rep(seq_len(V), length = length(y))
    folds <- sample(folds)
  }
  return(folds)
}
## get cv folds from a list created by CV.SuperLearner
get_cv_folds <- function(folds_lst) {
    V <- length(folds_lst)
    v_lst <- sapply(1:V, function(s) rep(s, length(folds_lst[[s]])), simplify = FALSE)
    joint_lst <- mapply(list, v_lst, folds_lst, SIMPLIFY = FALSE)
    folds_mat <- do.call(rbind, lapply(joint_lst, function(x) cbind(x[[1]], x[[2]])))
    folds <- folds_mat[order(folds_mat[, 2]), 1]
    return(folds)
}
## determine SL options based on outcome name
get_sl_options <- function(outcome_name, V) {
    if (grepl("dichot", outcome_name)) {
        sl_fam <- binomial()
        cv_ctrl_lst <- list(V = V, stratifyCV = TRUE)
        sl_method <- "tmp_method.CC_nloglik"
    } else {
        sl_fam <- gaussian()
        cv_ctrl_lst <- list(V = V)
        sl_method <- "tmp_method.CC_LS"
    }
    return(list(fam = sl_fam, ctrl = cv_ctrl_lst, method = sl_method))
}
## Determine vimp options based on outcome name
get_vimp_options <- function(outcome_name) {
    if (grepl("dichot", outcome_name)) {
        vimp_measure <- "auc"
    } else {
        vimp_measure <- "r_squared"
    }
    return(list(vimp_measure = vimp_measure))
}
## Make list of vimp objects
make_vimp_list <- function(var_groups, var_inds) {
    list_names <- c("grp_conditional", "grp_marginal", "ind_conditional", "ind_marginal")
    lst <- sapply(list_names, function(x) NULL, simplify = FALSE)
    return(lst)
}
## Make lists of cv objects
make_cv_lists <- function(folds_lst, full_vec, redu_vec) {
    folds <- get_cv_folds(folds_lst)
    ## make lists of the fitted values
    full_lst <- lapply(as.list(1:length(unique(folds))), function(x) full_vec[folds == x])
    redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) redu_vec[folds == x])
    return(list(folds = folds, full_lst = full_lst, redu_lst = redu_lst))
}
## Nice group names for vimp
vimp_nice_group_names <- function(nm_vec) {
    nice_names <- c("Cysteine counts", "Viral geometry", "Region-specific counts of PNG sites", "gp120 CD4 binding sites", "gp120 V2", "gp120 V3", "gp41 MPER")
    reference_nm_vec <- c("cysteines", "geometry", "glyco", "gp120_cd4bs", "gp120_v2", "gp120_v3", "gp41_mper")
    reference_positions <- apply(as.matrix(nm_vec), 1, function(x) grep(x, reference_nm_vec))
    return(nice_names[reference_positions])
}
vimp_nice_ind_names <- function(nm_vec) {
    no_hxb2 <- gsub("hxb2.", "", nm_vec)
    no_1mer <- gsub(".1mer", "", nm_vec)
    return(no_1mer)
}
## nice plotting names
vimp_plot_name <- function(vimp_str) {
    plot_nms <- rep(NA, length(vimp_str))
    plot_nms[grepl("iip", vimp_str)] <- "IIP"
    plot_nms[grepl("pc.ic50", vimp_str)] <- "IC-50"
    plot_nms[grepl("pc.ic80", vimp_str)] <- "IC-80"
    plot_nms[grepl("dichotomous.1", vimp_str)] <- "Estimated sensitivity"
    plot_nms[grepl("dichotomous.2", vimp_str)] <- "Multiple sensitivity"
    return(plot_nms)
}
vimp_plot_type <- function(str) {
    grp_nm <- rev(unlist(strsplit(str, "_", fixed = TRUE)))
    nice_grp <- gsub("grp", "Group Variable Importance", grp_nm)
    nice_grp_ind <- gsub("ind", "Individual Variable Importance", nice_grp)
    nice_grp_ind_marg <- gsub("marginal", "Marginal ", nice_grp_ind)
    nice_grp_ind_marg_cond <- gsub("conditional", "Conditional ", nice_grp_ind_marg)
    return(paste(nice_grp_ind_marg_cond, collapse = ""))
}
vimp_nice_rownames <- function(vimp_obj, cv = FALSE) {
    mat_s <- vimp_obj$mat$s
    lst_s <- vimp_obj$s
    indx_mat <- sapply(1:length(mat_s), function(x) which(mat_s[x] == lst_s))
    paste_ind <- 3
    if (cv) {
        paste_ind <- 4
    }
    tmp_nms <- unlist(lapply(strsplit(names(lst_s), "_", fixed = TRUE), function(x) paste(x[paste_ind:length(x)], collapse = "_")))
    return(tmp_nms[indx_mat])
}
