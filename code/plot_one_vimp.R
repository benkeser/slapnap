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
        vimp_est$group <- vimp_nice_rownames(vimp_obj)
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

vimp_plot_name <- function(vimp_obj) {
    row_nm <- rownames(vimp_obj$mat)[1]
    tmp_string <- strsplit(row_nm, "_", fixed = TRUE)[[1]][1]
    if (tmp_string == "iip") {
        return ("IIP")
    } else if (tmp_string == "pc.ic50") {
        return ("IC-50")
    } else if (tmp_string == "pc.ic80") {
        return ("IC-80")
    } else if (tmp_string == "dichotomous.1") {
        return ("Estimated sensitivity")
    } else {
        return ("Multiple sensitivity")
    }
}

vimp_nice_rownames <- function(vimp_obj) {
    row_nm <- rownames(vimp_obj$mat)
    return(unlist(lapply(strsplit(row_nm, "_", fixed = TRUE), function(x) tail(x, n = 1))))
}

### HERE'S THE CODE I TOOK FROM THE RMD. PUTTING HERE FOR THE TIME BEING
# # Variable importance

# ```{r load_vimp_objects}
# library("dplyr")
# library("ggplot2")
# library("cowplot")
# source("/home/lib/plot_one_vimp.R")
# continuous_outcome_vimp <- readRDS("/home/slfits/continuous_outcome_vimp.rds")
# continuous_outcome_cv_vimp <- readRDS("/home/slfits/continuous_outcome_cv_vimp.rds")
# # only load continuous ones, for now
# # binary_outcome_vimp <- readRDS("/home/slfits/binary_outcome_vimp.rds")
# # binary_outcome_cv_vimp <- readRDS("/home/slfits/binary_outcome_cv_vimp.rds")
# ```

# Figure \@ref(fig:continuous_cv_vimp) shows variable importance estimates for predicting continuous outcomes.

# ```{r plot-continuous-vimp, fig.cap = "Variable importance estimates for continuous outcomes"}
# x_lab <- expression(paste(R^2))
# x_lim <- c(0, 1)
# ## create a plot for each continuous outcome
# continuous_outcome_vimp_plots <- lapply(continuous_outcome_cv_vimp, plot_one_vimp, title = "Cross-validated variable importance")
# plot_grid(plotlist = continuous_outcome_vimp_plots)
# ```
