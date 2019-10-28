#!/usr/bin/env Rscript

## plot vimp for a single outcome/group combination
## @param vimp_obj the variable importance object
## @param title the title of the plot
## @param x_lim the x-axis limits
## @param x_lab the x-axis label
## @param lgnd_pos the legend position
## @param point_size the point size
## @param main_font_size the size of text
## @param cv whether or not this is cv importance
## @param num_plot the number of features to plot (in descending order)
plot_one_vimp <- function(vimp_obj, title = "Variable importance", x_lim = c(0, 1), x_lab = expression(paste(R^2)), lgnd_pos = c(0.1, 0.3), point_size = 5, main_font_size = 20, cv = FALSE, num_plot = 50) {
    if (!is.null(vimp_obj)) {
        ## get the variable importances
        if (!is.null(vimp_obj$mat)) {
            vimp_est <- vimp_obj$mat
            vimp_est$group <- vimp_nice_rownames(vimp_obj, cv = cv)
        } else {
            vimp_est <- cbind(est = vimp_obj$est, se = vimp_obj$se, cil = vimp_obj$cil, ciu = vimp_obj$ciu)
            print_s <- ifelse(length(vimp_obj$s) <= 10,
                          paste(vimp_obj$s, collapse = ", "),
                          paste(c(vimp_obj$s[1:10], "..."), collapse = ", "))
            vimp_est$group <- paste("s = ", print_s, sep = "")
        }
        ## plot by ordered vimp measure
        vimp_plot <- vimp_est %>%
            arrange(desc(est)) %>%
            filter(row_number() <= num_plot) %>%
            ggplot(aes(x = est, y = group[order(est)])) +
            geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
            geom_point(size = point_size) +
            ggtitle(title) +
            xlab(x_lab) +
            ylab("Feature group") +
            scale_x_continuous(breaks = round(seq(x_lim[1], x_lim[2], 0.1), 1),
                           labels = as.character(round(seq(x_lim[1], x_lim[2], 0.1), 1)),
                           limits = x_lim) +
            theme(legend.position = lgnd_pos,
                  text = element_text(size = main_font_size),
                  axis.title = element_text(size = main_font_size))

        return(vimp_plot)
    } else {
        print("No variable importance estimates to plot.")
    }
}

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
# vimp_plot_name <- function(vimp_obj) {
#     row_nm <- vimp_obj$mat$print_name
#     tmp_string_init <- strsplit(row_nm, "_", fixed = TRUE)[[1]]
#     tmp_string <- paste0(tmp_string_init[-length(tmp_string_init)], collapse = "_")
#     if (grepl("iip", tmp_string)) {
#         plot_nm <- "IIP"
#     } else if (grepl("pc.ic50", tmp_string)) {
#         plot_nm <- "IC-50"
#     } else if (grepl("pc.ic80", tmp_string)) {
#         plot_nm <- "IC-80"
#     } else if (grepl("dichotomous.1", tmp_string)) {
#         plot_nm <- "Estimated sensitivity"
#     } else {
#         plot_nm <- "Multiple Sensitivity"
#     }
#     return(plot_nm)
# }
#
# vimp_nice_rownames <- function(vimp_obj) {
#     row_nm <- vimp_obj$mat$print_name
#     return(unlist(lapply(strsplit(row_nm, "_", fixed = TRUE), function(x) tail(x, n = 1))))
# }
