#!/usr/bin/env -S Rscript --vanilla

#-------------------------------------------------#
# Create figure to show distribution of outcomes
#-------------------------------------------------#
# assume that when the code gets called there
# is an object called dat in memory and the relevant outcomes to
# plot are: log10.pc.ic50, log10.pc.ic80, iip, dichotomous.1,
# dichotomous.2


# TO DO:
# need to modify to get rid of Dataset 1/2 structure
# possible to align colors between continuous and dichotomous plots?

dat1 = read.csv(file = "~/data/data1.csv")
dat1$Dataset = "Dataset 1"
dat2 = read.csv(file = "~/data/data2.csv")
dat2$Dataset = "Dataset 2"

dat = rbind(dat1, dat2) %>% select(Dataset, ic50.geometric.mean.imputed.log10, ic80.geometric.mean.imputed.log10, neutralization.slope, ic50.censored, binding.dichotomous.sens.resis)

# Only Continuous outcomes: IC50 and IC80
dat_long = gather(dat[,c("Dataset", "ic50.geometric.mean.imputed.log10", "ic80.geometric.mean.imputed.log10")], outcome, value, ic50.geometric.mean.imputed.log10:ic80.geometric.mean.imputed.log10, factor_key=TRUE)
dat_long$outcome = ifelse(dat_long$outcome=="ic50.geometric.mean.imputed.log10", "log[10]*' IC'[50]",
                          ifelse(dat_long$outcome=="ic80.geometric.mean.imputed.log10", "log[10]*' IC'[80]", dat_long$outcome))
# create proportional stacked bars
ntext <- dat_long[!is.na(dat_long$value),] %>% group_by(Dataset, outcome) %>% dplyr::summarise(n = n())
dat_long = merge(dat_long, ntext, by=c("Dataset", "outcome"))
mid=0
high=0.5
low=-0.5

set.seed(20180430)
p1 = ggplot(dat_long, aes(x=Dataset, y=value, color=value)) + facet_grid(.~outcome, labeller = label_parsed) +
  geom_jitter(size=3, width = 0.25) + xlab("")  +
  scale_y_continuous(limits=c(-2,3.2)) +
  geom_boxplot(outlier.colour=NA, fill="gray70",alpha=0.2) +
  scale_color_gradient2(midpoint=mid, low="blue", mid="gray90", high="red", guide = "colourbar", space="Lab") +
  theme(legend.position="top", legend.title = element_blank(), strip.text = element_text(face="bold", size=18), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18), axis.text.x = element_text(size=17)) +
  ylab("Value") +
  geom_text(aes(label = paste0("n = ", n), y = 2.8), size=7, vjust = 1.5, color = "black")

# Only Continuous outcomes: Slope
dat_long = gather(dat[,c("Dataset", "neutralization.slope")], outcome, value, neutralization.slope, factor_key=TRUE)
dat_long$outcome = ifelse(dat_long$outcome=="neutralization.slope", "Neutralization*' Slope'", dat_long$outcome)
# create proportional stacked bars
ntext <- dat_long[!is.na(dat_long$value),] %>% group_by(Dataset, outcome) %>% dplyr::summarise(n = n())
dat_long = merge(dat_long, ntext, by=c("Dataset", "outcome"))
mid=1
high=quantile(dat_long[!is.na(dat_long$value),]$value, 0.75)
low=quantile(dat_long[!is.na(dat_long$value),]$value, 0.25)

set.seed(20180430)
p2 = ggplot(dat_long, aes(x=Dataset, y=value, color=value)) + facet_grid(.~outcome, labeller = label_parsed) +
  geom_jitter(size=3, width = 0.25) + xlab("")  +
  scale_y_continuous(limits=c(-2,3.2)) +
  geom_boxplot(outlier.colour=NA, fill="gray70",alpha=0.2) +
  scale_color_gradient2(midpoint=mid, low="red", mid="gray90", high="blue", guide = "colourbar", space="Lab") +
  theme(legend.position="top", legend.title = element_blank(), strip.text = element_text(face="bold", size=18), axis.text.x = element_text(size=17), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), plot.margin = unit( c(0,4,0,0) , "mm" )) +
  ylab("Value") + geom_text(aes(label = paste0("n = ", n), y = 2.8), size=7, vjust = 1.5, color = "black")

# Only Binary outcomes
dat_long = gather(dat[,c("Dataset", "ic50.censored", "binding.dichotomous.sens.resis")], outcome, value, ic50.censored:binding.dichotomous.sens.resis, factor_key=TRUE)
dat_long$outcome = ifelse(dat_long$outcome=="ic50.censored", "IC[50]*' Censored'",
                          ifelse(dat_long$outcome=="binding.dichotomous.sens.resis", "Sensitive/Resistant*' Only'", dat_long$outcome))
# create proportional stacked bars
proportion <- dat_long %>%
  group_by(Dataset, outcome, value) %>%
  dplyr::summarise(n = n()) %>% dplyr::mutate(freq = n * 100/ sum(n))

proportion$value = ifelse(is.na(proportion$value), as.integer(2), proportion$value)
proportion$freqPlace = ifelse(proportion$value==0 & proportion$outcome=="IC[50]*' Censored'", 84,
                              ifelse(proportion$value==1, 100,
                                     ifelse(proportion$value==0 & proportion$outcome=="Sensitive/Resistant*' Only'", 60, 83)))
proportion$value = factor(proportion$value, levels = c("1", "2", "0"))

set.seed(20180430)
p3 = ggplot(proportion, aes(factor(Dataset), freq, fill = value)) + facet_grid(~outcome, labeller = label_parsed) +
  geom_bar(stat = "identity", color = "grey40") +
  scale_fill_manual(values = c("red", "gray50", "blue"), labels=c("Resistant", "NA", "Sensitive")) +
  labs(fill = "value") +
  geom_text(aes(label = paste0(round(freq, 2), "%"), y = freqPlace), size=7, vjust = 1.5, color = "white") +
  geom_text(aes(label = paste0("n = ", n), y = freqPlace-7), size=7, vjust = 1.5, color = "white") +
  ylab("Percentage") + xlab("") +
  theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size=18), strip.text = element_text(face="bold", size=18), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18), axis.text.x = element_text(size=17)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

pdf("figures/distribution of outcomes.pdf", width = 20, height=6)
plot_grid(p3, p1, p2, labels=c('', '', ''), align="h", rel_widths=c(1,1,0.5), nrow=1)
grid.text("Resistant", x = unit(0.525, "npc"), y = unit(c(0.95), "npc"), gp=gpar(fontsize=18, col="black", fontface="bold"))
grid.text("Sensitive", x = unit(0.415, "npc"), y = unit(c(0.95), "npc"), gp=gpar(fontsize=18, col="black", fontface="bold"))
grid.text(".", x = unit(0.82, "npc"), y = unit(c(1.92), "npc"), gp=gpar(fontsize=1250, col="white", fontface="bold")) # Add white patch to hide colourbar
dev.off()