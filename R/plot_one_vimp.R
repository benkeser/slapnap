#!/usr/bin/env Rscript

## plot vimp for a single outcome/group combination
## @param vimp_obj the variable importance object
## @param title the title of the plot
## @param x_lim the x-axis limits
## @param x_lab the x-axis label
## @param lgnd_pos the legend position
## @param point_size the point size
## @param main_font_size the size of text
plot_one_vimp <- function(vimp_obj, title = "Variable importance", x_lim = c(0, 1), x_lab = expression(paste(R^2)), lgnd_pos = c(0.1, 0.3), point_size = 5, main_font_size = 20) {
    ## get the variable importances
    if (!is.null(vimp_obj$mat)) {
        vimp_est <- vimp_obj$mat
        vimp_est$group <- rownames(vimp_obj$mat)
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
}