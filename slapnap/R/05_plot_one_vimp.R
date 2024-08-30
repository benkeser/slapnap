#!/usr/bin/env Rscript

## plot vimp for a single outcome/group combination
#' @param vimp_obj the variable importance object
#' @param title the title of the plot
#' @param x_lim the x-axis limits
#' @param x_lab the x-axis label
#' @param lgnd_pos the legend position
#' @param point_size the point size
#' @param main_font_size the size of text
#' @param cv whether or not this is cv importance
#' @param num_plot the number of features to plot (in descending order)
#' @importFrom grid grob
#' @importFrom forcats fct_reorder
#' @importFrom tibble add_column
#' @export
plot_one_vimp <- function(vimp_obj, title = "Variable importance", x_lim = c(0, ifelse(!is.null(vimp_obj), ifelse(max(vimp_obj$mat$ciu) > 1, max(vimp_obj$mat$ciu) + 0.2, 1), 0)), x_lab = expression(paste(R^2)), lgnd_pos = c(0.1, 0.3), cv = FALSE, grp = TRUE, num_plot = 50, text_size = 9, threshold = 0.05, opts) {
    text_pos <- x_lim[2] - 0.05
    signif_pos <- x_lim[2] - 0.025
    ylab_txt <- ifelse(grp, "Feature group", "Feature")
    ylab_angle <- ifelse(grp, 15, 0)
    if (!is.null(vimp_obj)) {
        ## get the variable importances
        if (!is.null(vimp_obj$mat)) {
            vimp_est <- vimp_obj$mat
            nice_rownames <- vimp_nice_rownames(vimp_obj, cv = cv)
            if (grp) {
                vimp_est$group <- vimp_nice_group_names(nice_rownames)
            } else {
                vimp_est$group <- vimp_nice_ind_names(nice_rownames)
            }
        } else {
            vimp_est <- cbind(est = vimp_obj$est, se = vimp_obj$se, cil = vimp_obj$cil, ciu = vimp_obj$ciu)
            print_s <- ifelse(length(vimp_obj$s) <= 10,
                          paste(vimp_obj$s, collapse = ", "),
                          paste(c(vimp_obj$s[1:10], "..."), collapse = ", "))
            vimp_est$group <- paste("s = ", print_s, sep = "")
        }
        tmp_cis <- round(vimp_est[, c("cil", "ciu")], 3)
        text_cis <- apply(tmp_cis, 1, function(x) paste0("[", x[1], ", ", x[2], "]"))
        vimp_est <- add_column(vimp_est, text_ci = text_cis)
        vimp_est <- add_column(vimp_est, signif_p = ifelse(vimp_est$p_value < threshold, "*", ""))
        dim_plot <- ifelse(num_plot > 0, min(num_plot, dim(vimp_est)[1]), dim(vimp_est)[1])
        ## plot by ordered vimp measure
        vimp_plot <- vimp_est %>%
            arrange(desc(est)) %>%
            mutate(ord_group = fct_reorder(group, est)) %>%
            filter(row_number() <= num_plot) %>%
            ggplot(aes(x = est, y = ord_group)) +
            geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
            geom_point() +
            geom_text(aes(x = rep(text_pos, dim_plot), label = text_ci), hjust = "inward", color = "blue") +
            geom_text(aes(x = rep(signif_pos, dim_plot), label = signif_p), hjust = "inward", color = "blue") +
            ggtitle(title) +
            xlab(x_lab) +
            ylab(ylab_txt) +
            scale_x_continuous(breaks = round(seq(x_lim[1], x_lim[2], 0.1), 1),
                           labels = as.character(round(seq(x_lim[1], x_lim[2], 0.1), 1)),
                           limits = x_lim) +
            theme_bw() +
            theme(legend.position = lgnd_pos, legend.text = element_text(size = text_size), axis.text = element_text(size = text_size), plot.title.position = "plot", axis.text.y = element_text(angle = ylab_angle, vjust = 1, hjust = 1))
        return(vimp_plot)
    } else {
        return(grob(NULL))
    }
}
