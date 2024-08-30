#!/usr/bin/env -S Rscript --vanilla

## ----------------------------------------------------------------------------
## Individual-level algorithm-specific importance
## ----------------------------------------------------------------------------
# get individual p-value-based importance measures; create tables
col_idx <- geog_idx:ncol(dat)
max_features <- 15
if ("log10.pc.ic50" %in% outcome_names) {
    imp_ic50 <- get_all_importance("log10.pc.ic50", binary_outcome = FALSE,
                                   dir_loc = slfits_dir,
                                   dat = dat, which_cols = col_idx, opts = opts)
   ic50_tab <- get_importance_table(imp_ic50, max_features = max_features)
   direction_resis <- get_importance_resis(ic50_tab, dat = dat, which_outcome = "log10.pc.ic50")
   ic50_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("log10.pc.ic80" %in% outcome_names) {
    imp_ic80 <- get_all_importance("log10.pc.ic80", binary_outcome = FALSE,
                                   dir_loc = slfits_dir,
                                   dat = dat, which_cols = col_idx, opts = opts)
   ic80_tab <- get_importance_table(imp_ic80, max_features = max_features)
   direction_resis <- get_importance_resis(ic80_tab, dat = dat, which_outcome = "log10.pc.ic80")
   ic80_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("iip" %in% outcome_names) {
    imp_iip <- get_all_importance("iip", binary_outcome = FALSE,
                                   dir_loc = slfits_dir,
                                   dat = dat, which_cols = col_idx, opts = opts)
   iip_tab <- get_importance_table(imp_iip, max_features = max_features)
   direction_resis <- get_importance_resis(iip_tab, dat = dat, which_outcome = "iip")
   iip_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("dichotomous.1" %in% outcome_names) {
    imp_dichot1 <- get_all_importance("dichotomous.1", binary_outcome = TRUE,
                                   dir_loc = slfits_dir,
                                   dat = dat, which_cols = col_idx, opts = opts)
   dichot1_tab <- get_importance_table(imp_dichot1, max_features = max_features)
   direction_resis <- get_importance_resis(dichot1_tab, dat = dat, which_outcome = "dichotomous.1")
   dichot1_tab$direction <- ifelse(direction_resis, "Resistant", "Sensitive")
}
if ("dichotomous.2" %in% outcome_names) {
    imp_dichot2 <- get_all_importance("dichotomous.2", binary_outcome = TRUE,
                                   dir_loc = slfits_dir,
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
