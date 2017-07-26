#Input
log10_ks_cutoff <- -0.575

#Show
qplot(log10(syntelogs.complete$ks), geom="histogram", binwidth=.01)  + geom_vline(mapping = NULL, data = NULL, xintercept = log10_ks_cutoff, na.rm = FALSE, show.legend = NA)

#Output
ks_cutoff <- 10^log10_ks_cutoff
