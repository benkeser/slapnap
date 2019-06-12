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
    CD4_binding_sites <- c(124, 125, 126, 127, 196, 198, 279, 280, 281, 282, 283, 365, 366, 367, 368, 369, 370, 374, 425, 426, 427, 428, 429, 430, 431, 432, 455, 456, 457, 458, 459, 460, 461, 469, 471, 472, 473, 474, 475, 476, 477)
    PNG_sites <- c(130, 139, 143, 156, 187, 197, 241, 262, 289, 339, 355, 363, 406, 408, 410, 442, 448, 460, 462)
    all_glycosylationvars <- c('sequons.total.subset','sequons.total.env','sequons.total.gp120','sequons.total.v5','sequons.total.loop.d','sequons.total.loop.e','sequons.total.vrc01','sequons.total.cd4','sequons.total.sj.fence','sequons.total.sj.trimer')
    all_cysteinesvars <- c('cysteines.total.env','cysteines.total.gp120','cysteines.total.v5','cysteines.total.loop.d','cysteines.total.loop.e','cysteines.total.vrc01','cysteines.total.cd4')
    all_viralgeometryvars <- c('length.env','length.gp120','length.v5','length.v5.outliers','length.loop.d','length.loop.e','length.loop.e.outliers')
    all_stericbulkvars <- c('taylor.small.total.v5','taylor.small.total.loop.d','taylor.small.total.cd4')

    ## get all variable groups
    aa_cd4bs_vars <- get_variable_group(pred_names, CD4_binding_sites)
    aa_pngs_vars <- get_variable_group(pred_names, PNG_sites)
    aa_glyco_vars <- all_glycosylationvars[all_glycosylationvars %in% colnames(data)]
    aa_cysteine_vars <- all_cysteinesvars[all_cysteinesvars %in% colnames(data)]
    aa_geometry_vars <- all_viralgeometryvars[all_viralgeometryvars %in% colnames(data)]
    aa_stericbulk_vars <- all_stericbulkvars[all_stericbulkvars %in% colnames(data)]
    return(list(cd4bs = aa_cd4bs_vars, pngs = aa_pngs_vars,
                glyco = aa_glyco_vars, cysteines = aa_cysteine_vars,
                geometry = aa_geometry_vars, stericbulk = aa_stericbulk_vars))
}
