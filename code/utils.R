## various utility functions

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
    list_names <- c("conditional", "marginal", "individual")
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
