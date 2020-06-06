
## ----------------------------------------------------------------------------
## Population variable importance
## ----------------------------------------------------------------------------
num_pop_import <- 20  # the number of individual features to display in plots

## plotting things
x_lab_continuous <- expression(paste("Difference in ", R^2, sep = ""))
x_lab_binary <- expression(paste("Difference in ", AUC, sep = ""))
## read in importance results for each outcome, create a plot for each
## only return non-cv plots if cv = FALSE
imp_nms <- list(all_var_groups, all_var_groups, var_inds)
num_obs_fulls <- vector("numeric", length(outcome_names))
num_obs_reds <- vector("numeric", length(outcome_names))
names(num_obs_fulls) <- opts$outcomes
names(num_obs_reds) <- opts$outcomes
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    if (run_sl_vimp_bools$run_vimp[i]) {
        if (grepl("dichot", this_outcome_name)) {
            this_x_lab <- x_lab_binary
        } else {
            this_x_lab <- x_lab_continuous
        }
        ## importance results
        eval(parse(text = paste0(this_outcome_name, "_vimp_lst <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_vimp.rds'))")))
        eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_cv_vimp.rds'))")))
        eval(parse(text = paste0(this_outcome_name, "_outer_folds <- readRDS(file = paste0(slfits_dir, '", this_outcome_name, "_outer_folds.rds'))")))
        eval(parse(text = paste0("num_obs_fulls[i] <- sum(", this_outcome_name, "_outer_folds == 1)")))
        eval(parse(text = paste0("num_obs_reds[i] <- sum(", this_outcome_name, "_outer_folds == 2)")))
        ## make plots
        eval(parse(text = paste0("current_vimp_lst <- ", this_outcome_name, "_vimp_lst")))
        eval(parse(text = paste0("current_cv_vimp_lst <- ", this_outcome_name, "_cv_vimp_lst")))
        plot_title_expr <- vimp_plot_name_expr(this_outcome_name, one_nab = n_ab < 1)
        vimp_type_expr <- unlist(lapply(as.list(names(current_vimp_lst)), vimp_plot_type_expr))
        vimp_plot_titles <- unlist(lapply(vimp_type_expr, function(x) eval(bquote(expression(.(plot_title_expr)*.(x))))))
        grp_bool_lst <- as.list(grepl("grp", names(current_vimp_lst)))
        eval(parse(text = paste0(this_outcome_name, "_vimp_plots <- mapply(function(x, y, z) plot_one_vimp(x, title = y, x_lab = this_x_lab, cv = FALSE, grp = z, threshold = vimp_threshold, num_plot = num_pop_import, opts = opts), current_vimp_lst, vimp_plot_titles, grp_bool_lst, SIMPLIFY = FALSE)")))
        eval(parse(text = paste0(this_outcome_name, "_cv_vimp_plots <- mapply(function(x, y, z) plot_one_vimp(x, title = y, x_lab = this_x_lab, cv = TRUE, grp = z, threshold = vimp_threshold, num_plot = num_pop_import, opts = opts), current_cv_vimp_lst, vimp_plot_titles, grp_bool_lst, SIMPLIFY = FALSE)")))
    }
}

## make table for executive summary
ran_vimp_dichot1 <- run_sl_vimp_bools2$run_vimp[grepl("dichotomous.1", names(run_sl_vimp_bools2$run_vimp))]
ran_vimp_dichot2 <- run_sl_vimp_bools2$run_vimp[grepl("dichotomous.2", names(run_sl_vimp_bools2$run_vimp))]
ran_vimp_dichot1 <- ifelse(length(ran_vimp_dichot1) == 0, FALSE, ran_vimp_dichot1)
ran_vimp_dichot2 <- ifelse(length(ran_vimp_dichot2) == 0, FALSE, ran_vimp_dichot2)
if (opts$cvperf) {
    vimp_summary_tbl <- make_vimp_executive_summary_table(
        switch("log10.pc.ic50" %in% outcome_names + 1, NULL, log10.pc.ic50_cv_vimp_lst),
        switch("log10.pc.ic80" %in% outcome_names + 1, NULL, log10.pc.ic80_cv_vimp_lst),
        switch("iip" %in% outcome_names + 1, NULL, iip_cv_vimp_lst),
        switch((("dichotomous.1" %in% outcome_names) & (ran_vimp_dichot1)) + 1, NULL, dichotomous.1_cv_vimp_lst),
        switch((("dichotomous.2" %in% outcome_names) & (ran_vimp_dichot2)) + 1, NULL, dichotomous.2_cv_vimp_lst),
        threshold = vimp_threshold, outcome_names = all_outcome_names, cv = TRUE, opts = opts)
} else {
    vimp_summary_tbl <- make_vimp_executive_summary_table(switch("log10.pc.ic50" %in% outcome_names + 1, NULL, log10.pc.ic50_vimp_lst),
    switch("log10.pc.ic80" %in% outcome_names + 1, NULL, log10.pc.ic80_vimp_lst),
    switch("iip" %in% outcome_names + 1, NULL, iip_vimp_lst),
    switch((("dichotomous.1" %in% outcome_names) & (ran_vimp_dichot1)) + 1, NULL, dichotomous.1_vimp_lst),
    switch((("dichotomous.2" %in% outcome_names) & (ran_vimp_dichot2)) + 1, NULL, dichotomous.2_vimp_lst),
    threshold = vimp_threshold, outcome_names = all_outcome_names, cv = FALSE, opts = opts)
}
