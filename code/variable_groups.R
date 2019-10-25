## define variable groups for importance

## function to set a variable group
get_variable_group <- function(pred_names, hxb2_sites) {
    ## get only the sites from the pred_names vector
    aa_positions <- unlist(lapply(strsplit(pred_names, ".", fixed = TRUE), function(x) x[2]))

    ## get the ones matching the input sites
    site_character_vars <- pred_names[aa_positions %in% hxb2_sites]
    return(site_character_vars)
}

get_variable_groups <- function(data, pred_names) {
    ## set up sites
    ## CD4bs is feature group 2 from VRC01 paper plus additional sites identified by Adam
    gp120_cd4bs <- unique(c(c(124, 125, 126, 127, 196, 198, 279, 280, 281, 282, 283, 365, 366, 367, 368, 369, 370, 374, 425, 426, 427, 428, 429, 430, 431, 432, 455, 456, 457, 458, 459, 460, 461, 469, 471, 472, 473, 474, 475, 476, 477), c(197, 209, 279, 326, 369, 119, 120, 182, 204, 206, 207, 274, 304, 318, 369, 471), c(62, 64, 66, 207), c(61, 64, 197, 276, 362, 363, 386, 392, 462, 463)))
    gp120_v2_v2g_v2apex <- unique(c(157:196, c(121, 123, 124, 127, 197, 202, 203, 312, 315)))
    gp120_v3_v3g <- unique(c(296:334), c(380, 406, 408, 415, 419, 428, 441, 443, 471), c(156, 137))
    gp41_mper <- c(656:684, 609)
    all_glycosylationvars <- pred_nms[grepl("sequons", pred_nms)]
    all_cysteinesvars <- pred_nms[grepl("cysteine", pred_nms)]
    all_viralgeometryvars <- pred_nms[grepl("length", pred_nms)]

    ## get all variable groups
    aa_gp120_cd4bs_vars <- get_variable_group(pred_names, gp120_cd4bs)
    aa_gp120_v2_vars <- get_variable_group(pred_names, gp120_v2_v2g_v2apex)
    aa_gp120_v3_vars <- get_variable_group(pred_names, gp120_v3_v3g)
    aa_gp41_mper_vars <- get_variable_group(pred_names, gp41_mper)
    aa_glyco_vars <- all_glycosylationvars[all_glycosylationvars %in% colnames(data)]
    aa_cysteine_vars <- all_cysteinesvars[all_cysteinesvars %in% colnames(data)]
    aa_geometry_vars <- all_viralgeometryvars[all_viralgeometryvars %in% colnames(data)]
    return(list(gp120_cd4bs = aa_gp120_cd4bs_vars, gp120_v2 = aa_gp120_v2_vars,
                gp120_v3 = aa_gp120_v3_vars, gp41_mper = aa_gp41_mper_vars,
                glyco = aa_glyco_vars, cysteines = aa_cysteine_vars,
                geometry = aa_geometry_vars))
}
