#==================================================================================================#
## getTypeSort
#--------------------------------------------------------------------------------------------------#
getTypeSort <- function(homeologs.pairs, maize.walley.v4mapped.expression, foldCutoff, correlationCutoff) {
  ## Build dataframes necessary for sorting types
  # Retained duplicates are gene pairs where one is in maize1 and the other in maize2
  retained.duplicates <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  # Expression data for all genes in retained duplicates list as an avg of biological replicates
  data <-
    maize.walley.v4mapped.expression %>%
    select(geneID, sample, FPKM_avg) %>%
    subset(geneID %in% retained.duplicates$Maize1 | geneID %in% retained.duplicates$Maize2)

  # NA's mess with some of the aggregating later on, make them zeros
  data$FPKM_avg[is.na(data$FPKM_avg)] <- 0

  # A dead gene is a gene with less than 1 FPKM in all tissues. Find these first so the "zero" FPKMs don't mess with foldchanges calculations later on.
  dead.genes <-
    data %>%
    group_by(geneID) %>%
    summarise(max=max(FPKM_avg)) %>%
    ungroup() %>%
    subset(max < 1) %>%
    select(geneID)

  retained.duplicates.2dead <- retained.duplicates[retained.duplicates$Maize1 %in% dead.genes$geneID & retained.duplicates$Maize2 %in% dead.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.2dead)
  retained.duplicates.1dead <- retained.duplicates[retained.duplicates$Maize1 %in% dead.genes$geneID | retained.duplicates$Maize2 %in% dead.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.1dead)

  # Find genes with reciprocal dominance. (BED genes in the Pophaly paper)
  bed.genes <-
    retained.duplicates %>%
    inner_join(data, by=c("Maize1"="geneID")) %>%
    rename(FPKM_avg1=FPKM_avg) %>%
    inner_join(data, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(FPKM_avg2=FPKM_avg) %>%
    select(Maize1,Maize2,sample,FPKM_avg1,FPKM_avg2) %>%
    mutate(fold1=FPKM_avg1/FPKM_avg2, fold2=FPKM_avg2/FPKM_avg1) %>%
    mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., is.infinite(.), 0))) %>%
    group_by(Maize1,Maize2) %>%
    summarise(maxFoldChange1=max(fold1), maxFoldChange2=max(fold2)) %>%
    ungroup() %>%
    subset(maxFoldChange1 > foldCutoff & maxFoldChange2 > foldCutoff) %>%
    select(Maize1, Maize2) %>%
    gather(sub, geneID) %>%
    select(geneID)

  retained.duplicates.bed <- retained.duplicates[retained.duplicates$Maize1 %in% bed.genes$geneID & retained.duplicates$Maize2 %in% bed.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.bed)

  # Find genes with high correlation.  Make sure the "dead" genes and reciprocal dominance genes are out first, since one easy way to get high correlation
  # is to have 2 dead genes, and an instance of reciprocal dominance (requires 2 tissues at min) might occur in otherwise highly correlated genes
  # for all other tissues
  true.genes <-
    retained.duplicates %>%
    inner_join(data, by=c("Maize1"="geneID")) %>%
    rename(FPKM_avg1=FPKM_avg) %>%
    inner_join(data, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(FPKM_avg2=FPKM_avg) %>%
    select(Maize1,Maize2,sample,FPKM_avg1,FPKM_avg2) %>%
    group_by(Maize1,Maize2) %>%
    summarise(correlation=cor(FPKM_avg1, FPKM_avg2, method = "pearson")) %>%
    ungroup() %>%
    subset(correlation > correlationCutoff) %>%
    select(Maize1, Maize2) %>%
    gather(sub, geneID) %>%
    select(geneID)

  retained.duplicates.true <- retained.duplicates[retained.duplicates$Maize1 %in% true.genes$geneID & retained.duplicates$Maize2 %in% true.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.true)

  ## Pack for export from this function
  retained.duplicates$type <- "ambiguous"
  retained.duplicates.1dead$type <- "1dead"
  retained.duplicates.2dead$type <- "2dead"
  retained.duplicates.bed$type <- "bed"
  retained.duplicates.true$type <- "true"

  retained.duplicates <-
    retained.duplicates %>%
    bind_rows(retained.duplicates.1dead) %>%
    bind_rows(retained.duplicates.2dead) %>%
    bind_rows(retained.duplicates.bed) %>%
    bind_rows(retained.duplicates.true)

  return(retained.duplicates)
}

#==================================================================================================#
## getTypeSortProteins
#--------------------------------------------------------------------------------------------------#
getTypeSortProteins <- function(homeologs.pairs, maize.walley.abundance.v4, foldCutoff, correlationCutoff) {
  ## Build dataframes necessary for sorting types
  # Retained duplicates are gene pairs where one is in maize1 and the other in maize2
  retained.duplicates <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  # Abundance data for all proteins in retained duplicates list as an avg of biological replicates
  data <-
    maize.walley.abundance.v4 %>%
    select(geneID, sample, dNSAF_avg) %>%
    subset(geneID %in% retained.duplicates$Maize1 | geneID %in% retained.duplicates$Maize2)

  # NA's mess with some of the aggregating later on, make them zeros
  data$dNSAF_avg[is.na(data$dNSAF_avg)] <- 0

  # A dead gene is a gene with less than 1 FPKM in all tissues. Find these first so the "zero" FPKMs don't mess with foldchanges calculations later on.
  dead.genes <-
    data %>%
    group_by(geneID) %>%
    summarise(max=max(dNSAF_avg)) %>%
    ungroup() %>%
    subset(max < 1) %>%
    select(geneID)

  retained.duplicates.2dead <- retained.duplicates[retained.duplicates$Maize1 %in% dead.genes$geneID & retained.duplicates$Maize2 %in% dead.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.2dead)
  retained.duplicates.1dead <- retained.duplicates[retained.duplicates$Maize1 %in% dead.genes$geneID | retained.duplicates$Maize2 %in% dead.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.1dead)

  # Find genes with reciprocal dominance. (BED genes in the Pophaly paper)
  bed.genes <-
    retained.duplicates %>%
    inner_join(data, by=c("Maize1"="geneID")) %>%
    rename(dNSAF_avg1=dNSAF_avg) %>%
    inner_join(data, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(dNSAF_avg2=dNSAF_avg) %>%
    select(Maize1,Maize2,sample,dNSAF_avg1,dNSAF_avg2) %>%
    mutate(fold1=dNSAF_avg1/dNSAF_avg2, fold2=dNSAF_avg2/dNSAF_avg1) %>%
    mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., is.infinite(.), 0))) %>%
    group_by(Maize1,Maize2) %>%
    summarise(maxFoldChange1=max(fold1), maxFoldChange2=max(fold2)) %>%
    ungroup() %>%
    subset(maxFoldChange1 > foldCutoff & maxFoldChange2 > foldCutoff) %>%
    select(Maize1, Maize2) %>%
    gather(sub, geneID) %>%
    select(geneID)

  retained.duplicates.bed <- retained.duplicates[retained.duplicates$Maize1 %in% bed.genes$geneID & retained.duplicates$Maize2 %in% bed.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.bed)

  # Find genes with high correlation.  Make sure the "dead" genes and reciprocal dominance genes are out first, since one easy way to get high correlation
  # is to have 2 dead genes, and an instance of reciprocal dominance (requires 2 tissues at min) might occur in otherwise highly correlated genes
  # for all other tissues
  true.genes <-
    retained.duplicates %>%
    inner_join(data, by=c("Maize1"="geneID")) %>%
    rename(dNSAF_avg1=dNSAF_avg) %>%
    inner_join(data, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(dNSAF_avg2=dNSAF_avg) %>%
    select(Maize1,Maize2,sample,dNSAF_avg1,dNSAF_avg2) %>%
    group_by(Maize1,Maize2) %>%
    summarise(correlation=cor(dNSAF_avg1, dNSAF_avg2, method = "pearson")) %>%
    ungroup() %>%
    subset(correlation > correlationCutoff) %>%
    select(Maize1, Maize2) %>%
    gather(sub, geneID) %>%
    select(geneID)

  retained.duplicates.true <- retained.duplicates[retained.duplicates$Maize1 %in% true.genes$geneID & retained.duplicates$Maize2 %in% true.genes$geneID,]
  retained.duplicates <- setdiff(retained.duplicates, retained.duplicates.true)

  ## Pack for export from this function
  retained.duplicates$type <- "ambiguous"
  retained.duplicates.1dead$type <- "1dead"
  retained.duplicates.2dead$type <- "2dead"
  retained.duplicates.bed$type <- "bed"
  retained.duplicates.true$type <- "true"

  retained.duplicates <-
    retained.duplicates %>%
    bind_rows(retained.duplicates.1dead) %>%
    bind_rows(retained.duplicates.2dead) %>%
    bind_rows(retained.duplicates.bed) %>%
    bind_rows(retained.duplicates.true)

  return(retained.duplicates)
}

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
