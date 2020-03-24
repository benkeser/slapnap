## make vimp table for executive summary
make_vimp_executive_summary_table <- function(..., threshold = 0.05, outcome_names = NULL, cv = FALSE, opts) {
    ## capture everything entered before threshold
    L <- list(...)
    ## combine together as a tibble for conditional, marginal, individual
    output_tib_grp_cond <- L[[1]]$grp_conditional$mat
    output_tib_grp_marg <- L[[1]]$grp_marginal$mat
    output_tib_ind_cond <- L[[1]]$ind_conditional$mat
    output_tib_ind_marg <- L[[1]]$ind_marginal$mat
    if (length(L[[1]]) == 0) {

    } else {
        if ("cond" %in% opts$importance_grp) {
            output_tib_grp_cond <- dplyr::bind_cols(output_tib_grp_cond, outcome = rep(outcome_names[1], dim(output_tib_grp_cond)[1]), group = vimp_nice_rownames(L[[1]]$grp_conditional, cv = cv))
        }
        if ("marg" %in% opts$importance_grp) {
            output_tib_grp_marg <- dplyr::bind_cols(output_tib_grp_marg, outcome = rep(outcome_names[1], dim(output_tib_grp_marg)[1]), group = vimp_nice_rownames(L[[1]]$grp_marginal, cv = cv))
        }
        if ("cond" %in% opts$importance_ind) {
            output_tib_ind_cond <- dplyr::bind_cols(output_tib_ind_cond, outcome = rep(outcome_names[1], dim(output_tib_ind_cond)[1]), group = vimp_nice_rownames(L[[1]]$ind_conditional, cv = cv))
        }
        if ("marg" %in% opts$importance_ind) {
            output_tib_ind_marg <- dplyr::bind_cols(output_tib_ind_marg, outcome = rep(outcome_names[1], dim(output_tib_ind_marg)[1]), group = vimp_nice_rownames(L[[1]]$ind_marginal, cv = cv))
        }
    }
    if (length(L) > 1) {
        for (i in 2:length(L)) {
            if (length(L[[i]]) == 0) {

            } else {
                if ("cond" %in% opts$importance_grp) {
                    output_tib_grp_cond <- dplyr::bind_rows(output_tib_grp_cond, dplyr::bind_cols(L[[i]]$grp_conditional$mat, outcome = rep(outcome_names[i], dim(L[[i]]$grp_conditional$mat)[1]), group = vimp_nice_rownames(L[[i]]$grp_conditional, cv = cv)))
                }
                if ("marg" %in% opts$importance_grp) {
                    output_tib_grp_marg <- dplyr::bind_rows(output_tib_grp_marg, dplyr::bind_cols(L[[i]]$grp_marginal$mat, outcome = rep(outcome_names[i], dim(L[[i]]$grp_marginal$mat)[1]), group = vimp_nice_rownames(L[[i]]$grp_marginal, cv = cv)))
                }
                if ("cond" %in% opts$importance_ind) {
                    output_tib_ind_cond <- dplyr::bind_rows(output_tib_ind_cond, dplyr::bind_cols(L[[i]]$ind_conditional$mat, outcome = rep(outcome_names[i], dim(L[[i]]$ind_conditional$mat)[1]), group = vimp_nice_rownames(L[[i]]$ind_conditional, cv = cv)))
                }
                if ("marg" %in% opts$importance_ind) {
                    output_tib_ind_marg <- dplyr::bind_rows(output_tib_ind_marg, dplyr::bind_cols(L[[i]]$ind_marginal$mat, outcome = rep(outcome_names[i], dim(L[[i]]$ind_marginal$mat)[1]), group = vimp_nice_rownames(L[[i]]$ind_marginal, cv = cv)))
                }
            }
        }
    }
    ## determine the rankings of each group for each outcome
    if ("cond" %in% opts$importance_grp) {
        grp_cond_summary_tib <- output_tib_grp_cond %>%
            group_by(outcome) %>%
            mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_group_names(group)) %>%
            ungroup() %>%
            select(nice_outcome_name, nice_group_name, rank) %>%
            group_by(nice_group_name) %>%
            spread(key = nice_outcome_name, value = rank) %>%
            ungroup() %>%
            mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>%
            arrange(mn_rank)
    }
    if ("marg" %in% opts$importance_grp) {
        grp_marg_summary_tib <- output_tib_grp_marg %>%
            group_by(outcome) %>%
            mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_group_names(group)) %>%
            ungroup() %>%
            select(nice_outcome_name, nice_group_name, rank) %>%
            group_by(nice_group_name) %>%
            spread(key = nice_outcome_name, value = rank) %>%
            ungroup() %>%
            mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>%
            arrange(mn_rank)
    }
    if ("cond" %in% opts$importance_ind) {
        ind_cond_summary_tib <- output_tib_ind_cond %>%
            group_by(outcome) %>%
            mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_ind_names(group)) %>%
            ungroup() %>%
            select(nice_outcome_name, nice_group_name, rank) %>%
            group_by(nice_group_name) %>%
            spread(key = nice_outcome_name, value = rank) %>%
            ungroup() %>%
            mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>%
            arrange(mn_rank)
    }
    if ("marg" %in% opts$importance_ind) {
        ind_marg_summary_tib <- output_tib_ind_marg %>%
            group_by(outcome) %>%
            mutate(rank = row_number(), nice_outcome_name = vimp_plot_name(outcome), nice_group_name = vimp_nice_ind_names(group)) %>%
            ungroup() %>%
            select(nice_outcome_name, nice_group_name, rank) %>%
            group_by(nice_group_name) %>%
            spread(key = nice_outcome_name, value = rank) %>%
            ungroup() %>%
            mutate(mn_rank = select(., -matches("nice_group_name")) %>% rowMeans(.)) %>%
            arrange(mn_rank)
    }
    output_lst <- list(
        grp_cond = switch("cond" %in% opts$importance_grp, grp_cond_summary_tib, NULL),
        grp_marg = switch("marg" %in% opts$importance_grp, grp_marg_summary_tib, NULL),
        ind_cond = switch("cond" %in% opts$importance_ind, ind_cond_summary_tib, NULL),
        ind_marg = switch("marg" %in% opts$importance_ind, ind_marg_summary_tib, NULL)
    )
    return(output_lst)
}
