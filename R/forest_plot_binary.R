#####################################################
# Figures for ic50_censored outcome
#####################################################

# assume that when the code gets called there
# is an object called dat in memory and the relevant outcomes to 
# plot are: log10.pc.ic50, log10.pc.ic80, iip, dichotomous.1,
# dichotomous.2

# there will also be SuperLearner and CV.SuperLearner objects
# in memory called fit_ic50, fit_ic80, fit_iip, 
# fit_dichotomous1, fit_dichotomous2


# ic50_censored: cross-validation forest plots top 5 plus SL 
rm(fit)
# Load the cvma fit on data set 1 along with the dataset
load(file = "fit_cens_set1_v11_newest.RData")
first_row = create_forest_plots(fit, Y.cens, X, "dichotomous", "Dataset 1", "ic50.censored", "CV-AUC", "A", 0.5, 1)

rm(fit)
load(file = "fit_cens_set2_v11_newest.RData")
second_row = create_forest_plots(fit, Y2.cens, X2, "dichotomous", "Dataset 2", "ic50.censored", "CV-AUC", "B", 0.5, 1)

pdf("figures/ic50cens_forest_top5plusSL_crossVal.pdf", width = 20, height=15)
plot_grid(first_row, second_row, labels=c('', ''), ncol=1)
grid.text("Algorithm", x = unit(0.12, "npc"), y = unit(c(0.95,0.45), "npc"), gp=gpar(fontsize=25, col="black", fontface="bold"))
grid.text("Screen", x = unit(0.26, "npc"), y = unit(c(0.95,0.45), "npc"), gp=gpar(fontsize=25, col="black", fontface="bold"))
grid.text("AUC (95% CI)", x = unit(0.43, "npc"), y = unit(c(0.95,0.45), "npc"), gp=gpar(fontsize=25, col="black", fontface="bold"))
dev.off()


# ic50_censored: cross-validation CV-ROC plots AND predicted probability plots
rm(fit)
load(file = "fit_cens_set1_v11_newest.RData")
tbl_forPredProb_set1 = create_CVROC_plots("Dataset 1", fit, SL.library, Y.cens, X, "figures/ic50cens_CVROC_top3_crossVal.pdf")
tab = dataprep(fit)
p1 = predicted_Probability_plot_crossVal(as.data.table(tab), tbl_forPredProb_set1)

rm(fit)
load(file = "fit_cens_set2_v11_newest.RData")
tbl_forPredProb_set2 = create_CVROC_plots("Dataset 2", fit, SL.library, Y2.cens, X2, "figures/ic50cens_CVROC_top3_crossVal.pdf")
tab2 = dataprep(fit)
p2 = predicted_Probability_plot_crossVal(as.data.table(tab2), tbl_forPredProb_set2)

pdf("figures/ic50cens_predProb_top3_crossVal.pdf", width = 13, height=15)
require(cowplot)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, align = 'v', label_size=30)
dev.off()

# TO DO: I couldn't tell what plot this made

# ic50_censored: cross-validation Forest plots with all results 
rm(fit)
load(file = "fit_cens_set1_v11_newest.RData")
tab2 = dataprep(fit)
forestplot_allmodels_dichotomousOutcome(tab2, "Dataset 1", "figures/cens_set1_forestplot_CVAUCs_crossval.pdf", 0.3, 1)

rm(fit)
load(file = "fit_cens_set2_v11_newest.RData")
tab2 = dataprep(fit)
forestplot_allmodels_dichotomousOutcome(tab2, "Dataset 2", "figures/cens_set2_forestplot_CVAUCs_crossval.pdf", 0.3, 1)
