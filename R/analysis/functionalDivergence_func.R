#--------------------------------------------------------------------------------------------------#
# A function that can catch all alternating gene dominances 2/22/2018
#--------------------------------------------------------------------------------------------------#

plotRetainedDuplicatesExpression <- function(retained.duplicates, maize.expression.sample.avg, index) {
  test.index <- index
  test.expression <- maize.expression.sample.avg
  test.expression$FPKM_avg[is.na(test.expression$FPKM_avg)] <- 0

  test.subset <-
    retained.duplicates[test.index,] %>%
    inner_join(test.expression, by=c("Maize1"="geneID")) %>%
    rename(FPKM_avg1=FPKM_avg) %>%
    inner_join(test.expression, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(FPKM_avg2=FPKM_avg)

  test.subset.cor <-
    test.subset %>%
    select(sample, FPKM_avg1, FPKM_avg2)

  test.gene1 <- retained.duplicates[test.index,2]
  test.gene2 <- retained.duplicates[test.index,3]
  test.cor <- cor(x=test.subset.cor$FPKM_avg1, y=test.subset.cor$FPKM_avg2)

  test.subset.melt <-
    test.subset %>%
    select(sample, FPKM_avg1, FPKM_avg2) %>%
    gather(key = subgenome, "FPKM_avg", 2:3)

  test.subset.melt <- sortOnWalleyTissues(df = test.subset.melt, sampleCol = "sample")

  g <- ggplot(data=test.subset.melt, aes(x=sample, y=FPKM_avg, group=subgenome, color=subgenome)) + geom_line() + geom_point() +
    annotate("rect", xmin = .5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "darkgoldenrod4") +
    annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Green") +
    annotate("rect", xmin = 8.5, xmax = 12.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "DarkGreen") +
    annotate("rect", xmin = 12.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Yellow") +
    annotate("rect", xmin = 20.5, xmax = 23.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "darkorange") +
    labs(title="FPKM across tissues for putative retained duplicate in B73", y="FPKM Average(Bioreps)", x="Sample", subtitle = paste0(test.gene1, " and ", test.gene2,": Pearson correlation = ", test.cor)) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(g)
}

plotRetainedDuplicatesAbundance <- function(retained.duplicates, maize.protein.abundance.sample.avg, index) {
  test.index <- index
  test.abundance <- maize.protein.abundance.sample.avg
  test.abundance$dNSAF_avg[is.na(test.abundance$dNSAF_avg)] <- 0

  test.subset <-
    retained.duplicates[test.index,] %>%
    inner_join(test.abundance, by=c("Maize1"="geneID")) %>%
    rename(dNSAF_avg1=dNSAF_avg) %>%
    inner_join(test.abundance, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(dNSAF_avg2=dNSAF_avg)

  test.subset.cor <-
    test.subset %>%
    select(sample, dNSAF_avg1, dNSAF_avg2)

  test.gene1 <- retained.duplicates[test.index,2]
  test.gene2 <- retained.duplicates[test.index,3]
  test.cor <- cor(x=test.subset.cor$dNSAF_avg1, y=test.subset.cor$dNSAF_avg2)

  test.subset.melt <-
    test.subset %>%
    select(sample, dNSAF_avg1, dNSAF_avg2) %>%
    gather(key = subgenome, "dNSAF_avg", 2:3)

  test.subset.melt <- sortOnWalleyTissues(df = test.subset.melt, sampleCol = "sample")

  g <- ggplot(data=test.subset.melt, aes(x=sample, y=dNSAF_avg, group=subgenome, color=subgenome)) + geom_line() + geom_point() +
    annotate("rect", xmin = .5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "darkgoldenrod4") +
    annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Green") +
    annotate("rect", xmin = 8.5, xmax = 12.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "DarkGreen") +
    annotate("rect", xmin = 12.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Yellow") +
    annotate("rect", xmin = 20.5, xmax = 23.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "darkorange") +
    labs(title="dNSAF across tissues for putative retained duplicate in B73", y="dNSAF Average(Bioreps)", x="Sample", subtitle = paste0(test.gene1, " and ", test.gene2,": Pearson correlation = ", test.cor)) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  return(g)
}

sortOnWalleyTissues <- function(df, sampleCol) {
  df[sampleCol] <- factor(df[[sampleCol]], c("Primary Root 5 Days",
                                             "Root - Cortex 5 Days",
                                             "Root - Elongation Zone 5 Days",
                                             "Root - Meristem Zone 5 Days",
                                             "Secondary Root 7-8 Days",
                                             "6-7 internode",
                                             "7-8 internode",
                                             "Vegetative Meristem 16-19 Day",
                                             "Leaf Zone 1 (Symmetrical)",
                                             "Leaf Zone 2 (Stomatal)",
                                             "Leaf Zone 3 (Growth)",
                                             "Mature Leaf 8",
                                             "Ear Primordium 2-4 mm",
                                             "Ear Primordium 6-8 mm",
                                             "Embryo 20 DAP",
                                             "Embryo 38 DAP",
                                             "Endosperm 12 DAP",
                                             "Endosperm Crown 27 DAP",
                                             "Germinatin Kernels 2 DAI",
                                             "Pericarp/Aleurone 27 DAP",
                                             "B73 Mature Pollen",
                                             "Female Spikelet Collected on day as silk",
                                             "Silk"))
  return(df)
}
