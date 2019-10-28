## make vimp table for executive summary
make_vimp_executive_summary_table <- function(..., threshold = 0.05, outcome_names = NULL, cv = FALSE) {
    ## capture everything entered before threshold
    L <- list(...)
    ## combine together as a tibble for conditional, marginal, individual
    output_tib_cond <- L[[1]]$conditional$mat
    output_tib_marg <- L[[1]]$marginal$mat
    output_tib_indi <- L[[1]]$individual$mat
    output_tib_cond <- dplyr::bind_cols(output_tib_cond, outcome = rep(outcome_names[1], dim(output_tib_cond)[1]), group = vimp_nice_rownames(L[[1]]$conditional, cv = cv))
    output_tib_marg <- dplyr::bind_cols(output_tib_marg, outcome = rep(outcome_names[1], dim(output_tib_marg)[1]), group = vimp_nice_rownames(L[[1]]$marginal, cv = cv))
    output_tib_indi <- dplyr::bind_cols(output_tib_indi, outcome = rep(outcome_names[1], dim(output_tib_indi)[1]), group = vimp_nice_rownames(L[[1]]$individual, cv = cv))
    if (length(L) > 1) {
        for (i in 2:length(L)) {
            output_tib_cond <- dplyr::bind_rows(output_tib_cond, dplyr::bind_cols(L[[i]]$conditional$mat, outcome = rep(outcome_names[i], dim(L[[i]]$conditional$mat)[1]), group = vimp_nice_rownames(L[[i]]$conditional, cv = cv)))
            output_tib_marg <- dplyr::bind_rows(output_tib_marg, dplyr::bind_cols(L[[i]]$marginal$mat, outcome = rep(outcome_names[i], dim(L[[i]]$marginal$mat)[1]), group = vimp_nice_rownames(L[[i]]$marginal, cv = cv)))
            if (!is.null(L[[i]]$individual)) {
                output_tib_indi <- dplyr::bind_rows(output_tib_indi, dplyr::bind_cols(L[[i]]$individual$mat, outcome = rep(outcome_names[i], dim(L[[i]]$individual$mat)[1]), group = vimp_nice_rownames(L[[1]]$individual, cv = cv)))
            }
        }
    }
    ## determine the ones that are significantly associated, make a tibble
    signif_cond <- output_tib_cond %>%
        filter(p_value <= threshold) %>%
        mutate(import_type = "Group: conditional")
    signif_marg <- output_tib_marg %>%
        filter(p_value <= threshold) %>%
        mutate(import_type = "Group: marginal")
    signif_indi <- output_tib_indi %>%
        filter(p_value <= threshold) %>%
        mutate(import_type = "Individual: marginal")
    output_tib <- dplyr::bind_rows(signif_cond, signif_marg, signif_indi)
    return(output_tib)
}
