## various utility functions

## make outer folds for VIM hypothesis test (based on sample splitting)
make_folds <- function(y, V, stratified = TRUE) {
  if (stratified) {
    y_1 <- y == 1
    y_0 <- y == 0
    folds_1 <- rep(seq_len(V), length = sum(y_1))
    folds_1 <- sample(folds_1)
    folds_0 <- rep(seq_len(V), length = sum(y_0))
    folds_0 <- sample(folds_0)
    folds <- vector("numeric", length(draw$y_cat))
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
get_sl_options <- function(outcome_name) {
    if (grepl("dichot", outcome_name)) {
        sl_fam <- "binomial"
        cv_ctrl_lst <- list(V = V, stratifyCV = TRUE)
        sl_method <- "tmp_method.CC_nloglik"
    } else {
        sl_fam <- "gaussian"
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
