####################################################################################################
## Project:         SubGenomes project
## Script purpose:  Create functions to be sourced and used for this project
##
## Input:
## Output:
## Date: 2017-09-14
## Author: Jesse R. Walsh
####################################################################################################

## Slicing data
#--------------------------------------------------------------------------------------------------#
getExpressionByExperiment <- function(maize.expression.all, homeologs.pairs, experiment) {
  if (experiment <= 1 | experiment >= 70) {
    return(NA)
  }
  expressedGenes <- maize.expression.all[,c(1,experiment)]
  homeolog.expression <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    inner_join(expressedGenes, by=c("Maize1"="geneID")) %>%
    inner_join(expressedGenes, by=c("Maize2"="geneID"))

    names(homeolog.expression)[3] <- "FPKM_maize1"
    names(homeolog.expression)[4] <- "FPKM_maize2"

    return(homeolog.expression)
}

## experimentList = list of column positions from maize.expression.all (the dataset used currently has data in columns 2:69, but only about 10 display well at a time)
## Factor = how many times greater expression has to be to be considered "dominant" -> (i.e. is A >= factor*B)
## omitNA = if FALSE, NA values in the expression dataset for either homeolog result in a "Not Expressed",
##          if TRUE, NA values are converted to zero's, so any expression for one homeolog will automatically dominate the other and NA in both results in a "Neither"
##          this setting does not stop NAs in the case where the gene pair doesn't exist in the expression set
graphGenePairExpressionsByExperiement <- function(maize.expression.all, homeologs.pairs, experimentList, factor, includeNA) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  nPairs <- nrow(genePairs)

  for (i in experimentList) {
    data <- getExpressionByExperiment(maize.expression.all, homeologs.pairs, i)
    if (includeNA) {
      data$FPKM_maize1[is.na(data$FPKM_maize1)] <- 0
      data$FPKM_maize2[is.na(data$FPKM_maize2)] <- 0
    }
    data$dominance <- NA
    data$dominance[is.na(data$FPKM_maize1) | is.na(data$FPKM_maize2)] <- "Not Expressed"
    data$dominance[data$FPKM_maize1 >= factor*data$FPKM_maize2] <- "Maize1"
    data$dominance[data$FPKM_maize2 >= factor*data$FPKM_maize1] <- "Maize2"
    data$dominance[data$FPKM_maize1 <= factor*data$FPKM_maize2 & data$FPKM_maize2 <= factor*data$FPKM_maize1] <- "Neither"

    data <-
      data %>%
      select(Maize1, Maize2, dominance)
    names(data)[names(data)=="dominance"] <- names(maize.expression.all)[i]

    genePairs <-
      genePairs %>%
      left_join(data, by = c("Maize1" = "Maize1", "Maize2" = "Maize2"))
  }

  data <- genePairs
  data <- data[,c(-1,-2)]
  freq <- table(col(data), as.matrix(data))
  Experiment <- names(data)
  data <- data.frame(cbind(freq), Experiment)
  data.m <-gather(data, Variable, Frequency, -Experiment)

  plot <-
    ggplot(data.m, aes(Experiment, Frequency/nPairs*100)) +
    geom_bar(aes(fill=Variable), position="dodge", stat="identity") +
    labs(y = paste0("Percent of gene pairs (n=",nPairs,")"),
         x = "Experiment",
         title = "Which homeolog dominates expression?"
    )

  return(plot)
}


## Utility Functions
#--------------------------------------------------------------------------------------------------#
saveMyObjects <- function() {
  rm(list = ls())
  source("~/git/SubGenomes/R/loadData.R")
  source("~/git/SubGenomes/R/cleanData.R")
  save.image("~/git/SubGenomes/Data/SavedObjects/loadedData.RData")
}

loadMyObjects <- function() {
  load("~/git/SubGenomes/Data/SavedObjects/loadedData.RData", .GlobalEnv)
}
