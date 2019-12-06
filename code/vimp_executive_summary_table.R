## make vimp table for executive summary
make_vimp_executive_summary_table <- function(..., threshold = 0.05, outcome_names = NULL, cv = FALSE, indi = FALSE) {
    ## capture everything entered before threshold
    L <- list(...)
    ## combine together as a tibble for conditional, marginal, individual
    output_tib_cond <- L[[1]]$conditional$mat
    output_tib_marg <- L[[1]]$marginal$mat
    output_tib_indi <- L[[1]]$individual$mat
    output_tib_cond <- dplyr::bind_cols(output_tib_cond, outcome = rep(outcome_names[1], dim(output_tib_cond)[1]), group = vimp_nice_rownames(L[[1]]$conditional, cv = cv))
    output_tib_marg <- dplyr::bind_cols(output_tib_marg, outcome = rep(outcome_names[1], dim(output_tib_marg)[1]), group = vimp_nice_rownames(L[[1]]$marginal, cv = cv))
    if (indi) {
        output_tib_indi <- dplyr::bind_cols(output_tib_indi, outcome = rep(outcome_names[1], dim(output_tib_indi)[1]), group = vimp_nice_rownames(L[[1]]$individual, cv = cv))
    }

    if (length(L) > 1) {
        for (i in 2:length(L)) {
            output_tib_cond <- dplyr::bind_rows(output_tib_cond, dplyr::bind_cols(L[[i]]$conditional$mat, outcome = rep(outcome_names[i], dim(L[[i]]$conditional$mat)[1]), group = vimp_nice_rownames(L[[i]]$conditional, cv = cv)))
            output_tib_marg <- dplyr::bind_rows(output_tib_marg, dplyr::bind_cols(L[[i]]$marginal$mat, outcome = rep(outcome_names[i], dim(L[[i]]$marginal$mat)[1]), group = vimp_nice_rownames(L[[i]]$marginal, cv = cv)))
            if (!is.null(L[[i]]$individual) & indi) {
                output_tib_indi <- dplyr::bind_rows(output_tib_indi, dplyr::bind_cols(L[[i]]$individual$mat, outcome = rep(outcome_names[i], dim(L[[i]]$individual$mat)[1]), group = vimp_nice_rownames(L[[1]]$individual, cv = cv)))
            }
        }
    }
    ## determine the rankings of each group for each outcome
    cond_tib_with_ranking <- output_tib_cond %>%
        group_by(outcome) %>%
        mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_group_names(group)) %>% 
        ungroup()
    marg_tib_with_ranking <- output_tib_marg %>%
        group_by(outcome) %>%
        mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_group_names(group)) %>% 
        ungroup()
    ## get tables with variable groups in the rows, outcome in the columns, and ranks in the cells
    ## order by mean rank
    cond_summary_tib <- cond_tib_with_ranking %>%
        select(nice_outcome_name, nice_group_name, rank) %>%
        group_by(nice_group_name) %>%
        spread(key = nice_outcome_name, value = rank) %>% 
        ungroup() %>% 
        mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>% 
        arrange(mn_rank)
    marg_summary_tib <- marg_tib_with_ranking %>%
        select(nice_outcome_name, nice_group_name, rank) %>%
        group_by(nice_group_name) %>%
        spread(key = nice_outcome_name, value = rank) %>% 
        ungroup() %>% 
        mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>% 
        arrange(mn_rank)
    if (indi) {
        indi_tib_with_ranking <- output_tib_indi %>%
            group_by(outcome) %>%
            mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_group_names(group))
        indi_summary_tib <- indi_tib_with_ranking %>%
            select(nice_outcome_name, nice_group_name, rank) %>%
            group_by(nice_group_name) %>%
            spread(key = nice_outcome_name, value = rank) %>% 
            ungroup() %>% 
            mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>% 
            arrange(mn_rank)
        output_lst <- list(cond = cond_summary_tib, marg = marg_summary_tib, indi = indi_summary_tib)
    } else {
        output_lst <- list(cond = cond_summary_tib, marg = marg_summary_tib, indi = NA)
    }
    return(output_lst)
}
