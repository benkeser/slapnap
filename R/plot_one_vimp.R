#!/usr/bin/env Rscript

## plot vimp for a single outcome/group combination
## @param vimp_est the estimated importance (along with CIs, etc.)
## @param x_lim the x-axis limits
## @param x_lab the x-axis label
## @param lgnd_pos the legend position
## @param point_size the point size
## @param main_font_size the size of text
plot_one_vimp <- function(vimp_est, x_lim = c(0, 1), x_lab = expression(paste(R^2)), lgnd_pos = c(), point_size = 5, main_font_size = 20) {
    ## plot by ordered vimp measure
    vimp_plot <- vimp_est %>%
        arrange(desc(measure)) %>%
        ggplot(aes(x = measure, y = group[order(measure)])) +
        geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
        geom_point(size = point_size) +
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