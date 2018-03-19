plotExpressionGrid <- function(homeologs.pairs, maize.expression.sample.avg) {
  data <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    inner_join(maize.expression.sample.avg, by=c("Maize1"="geneID")) %>%
    inner_join(maize.expression.sample.avg, by=c("Maize2"="geneID", "sample"="sample")) %>%
    subset(!is.na(FPKM_avg.x) & !is.na(FPKM_avg.y))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"
  # data <- subset(data, sample %in% c("Silk", "Embryo 20 DAP", "Leaf Zone 2 (Stomatal)"))

  plot <-
    ggplot(data, aes(x=Value_maize1, y=Value_maize2)) +
    geom_point(alpha=.1, aes(color = sample)) +
    geom_abline(mapping = null, data = null, slope = 1, intercept = 0) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
    scale_x_log10() + scale_y_log10() +
    labs(
      title = paste0("Scatterplot of Gene Expression for\nExpressed Retained Duplicates"),
      x = "Expression for Maize1 Gene (log scale)",
      y = "Expression for Maize2 Gene (log scale)"
    ) +
    theme(legend.position = "none") + facet_wrap(~sample, nrow=8)
  return(plot)
}

plotAbundanceGrid <- function(homeologs.pairs, maize.protein.abundance.sample.avg) {
  data <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    inner_join(maize.protein.abundance.sample.avg, by=c("Maize1"="geneID")) %>%
    inner_join(maize.protein.abundance.sample.avg, by=c("Maize2"="geneID", "sample"="sample")) %>%
    subset(!is.na(dNSAF_avg.x) & !is.na(dNSAF_avg.y))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"
  # data <- subset(data, sample %in% c("Silk", "Embryo 20 DAP", "Leaf Zone 2 (Stomatal)"))

  plot <-
    ggplot(data, aes(x=Value_maize1, y=Value_maize2)) +
    geom_point(alpha=.1, aes(color = sample)) +
    geom_abline(mapping = null, data = null, slope = 1, intercept = 0) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
    scale_x_log10() + scale_y_log10() +
    labs(
      title = paste0("Scatterplot of Protein Abundance for\nExpressed Retained Duplicates"),
      x = "Abundance for Maize1 Gene (log scale)",
      y = "Abundance for Maize2 Gene (log scale)"
    ) +
    theme(legend.position = "none") + facet_wrap(~sample, nrow=8)
  return(plot)
}

plotExpressionGrid(homeologs.pairs = homeologs.pairs, maize.expression.sample.avg = maize.expression.sample.avg)
plotAbundanceGrid(homeologs.pairs = homeologs.pairs, maize.protein.abundance.sample.avg = maize.protein.abundance.sample.avg)
