
## ----------------------------------------------------------------------------
## Population variable importance
## ----------------------------------------------------------------------------
num_pop_import <- 20  # the number of individual features to display in plots

## plotting things
x_lab_continuous <- expression(paste("Difference in ", R^2, sep = ""))
x_lim_continuous <- c(0, 1.2)
x_lab_binary <- expression(paste("Difference in ", AUC, sep = ""))
x_lim_binary <- c(0, 1.2)
## read in importance results for each outcome, create a plot for each
## only return non-cv plots if cv = FALSE
imp_nms <- list(all_var_groups, all_var_groups, var_inds)
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    if (grepl("dichot", this_outcome_name)) {
        this_x_lab <- x_lab_binary
        this_x_lim <- x_lim_binary
    } else {
        this_x_lab <- x_lab_continuous
        this_x_lim <- x_lim_continuous
    }
    ## importance results
    eval(parse(text = paste0(this_outcome_name, "_vimp_lst <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_vimp.rds'))")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_cv_vimp.rds'))")))
    eval(parse(text = paste0(this_outcome_name, "_outer_folds <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_outer_folds.rds'))")))
    eval(parse(text = paste0("num_obs_full <- sum(", this_outcome_name, "_outer_folds == 1)")))
    eval(parse(text = paste0("num_obs_red <- sum(", this_outcome_name, "_outer_folds == 2)")))
    ## make plots
    eval(parse(text = paste0("current_vimp_lst <- ", this_outcome_name, "_vimp_lst")))
    eval(parse(text = paste0("current_cv_vimp_lst <- ", this_outcome_name, "_cv_vimp_lst")))
    vimp_plot_titles <- paste0(vimp_plot_name(this_outcome_name), ": ", unlist(lapply(as.list(names(current_vimp_lst)), vimp_plot_type)))
    eval(parse(text = paste0(this_outcome_name, "_vimp_plots <- mapply(function(x, y) plot_one_vimp(x, title = y, x_lab = this_x_lab, cv = FALSE, num_plot = num_pop_import), current_vimp_lst, vimp_plot_titles, SIMPLIFY = FALSE)")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_plots <- mapply(function(x, y) plot_one_vimp(x, title = y, x_lab = this_x_lab, cv = TRUE, num_plot = num_pop_import), current_cv_vimp_lst, vimp_plot_titles, SIMPLIFY = FALSE)")))
}

## make table for executive summary
vimp_threshold <- 0.05
if (opts$cvperf) {
    vimp_summary_tbl <- make_vimp_executive_summary_table(
        switch("log10.pc.ic50" %in% outcome_names, log10.pc.ic50_cv_vimp_lst, NULL),
        switch("log10.pc.ic80" %in% outcome_names, log10.pc.ic80_cv_vimp_lst, NULL),
        switch("iip" %in% outcome_names, iip_cv_vimp_lst, NULL),
        switch("dichotomous.1" %in% outcome_names, dichotomous.1_cv_vimp_lst, NULL),
        switch("dichotomous.2" %in% outcome_names, dichotomous.2_cv_vimp_lst, NULL),
        threshold = vimp_threshold, outcome_names = outcome_names, cv = TRUE, opts = opts)
} else {
    vimp_summary_tbl <- make_vimp_executive_summary_table(switch("log10.pc.ic50" %in% outcome_names, log10.pc.ic50_vimp_lst, NULL),
    switch("log10.pc.ic80" %in% outcome_names, log10.pc.ic80_vimp_lst, NULL),
    switch("iip" %in% outcome_names, iip_vimp_lst, NULL),
    switch("dichotomous.1" %in% outcome_names, dichotomous.1_vimp_lst, NULL),
    switch("dichotomous.2" %in% outcome_names, dichotomous.2_vimp_lst, NULL),
    threshold = vimp_threshold, outcome_names = outcome_names, cv = FALSE, opts = opts)
}
