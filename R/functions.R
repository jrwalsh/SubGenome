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
# getExpressionByExperiment <- function(maize.expression.all, homeologs.pairs, experiment) {
#   if (experiment <= 1 | experiment >= 70) {
#     return(NA)
#   }
#   expressedGenes <- maize.expression.all[,c(1,experiment)]
#   homeolog.expression <-
#     homeologs.pairs %>%
#     subset(Maize1 != "" & Maize2 != "") %>%
#     select(Maize1, Maize2) %>%
#     inner_join(expressedGenes, by=c("Maize1"="geneID")) %>%
#     inner_join(expressedGenes, by=c("Maize2"="geneID"))
#
#     names(homeolog.expression)[3] <- "FPKM_maize1"
#     names(homeolog.expression)[4] <- "FPKM_maize2"
#
#     return(homeolog.expression)
# }

## maize.data = long-form expression/protein data
## Factor = how many times greater expression has to be to be considered "dominant" -> (i.e. is A >= factor*B)
## omitNA = if FALSE, NA values in the expression dataset for either homeolog result in a "Not Expressed",
##          if TRUE, NA values are converted to zero's, so any expression for one homeolog will automatically dominate the other and NA in both results in a "Neither"
##          this setting does not stop NAs in the case where the gene pair doesn't exist in the expression set
graphGenePairDominanceByExperiment <- function(maize.data, homeologs.pairs, factor, includeNA) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  data <-
    genePairs %>%
    inner_join(maize.data, by=c("Maize1"="geneID")) %>%
    inner_join(maize.data, by=c("Maize2"="geneID", "Sample"="Sample"))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"

  ## We only count gene pairs where both genes were at least tested in the data set
  nPairs <-
    data %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    nrow()

  if (includeNA) {
    data$Value_maize1[is.na(data$Value_maize1)] <- 0
    data$Value_maize2[is.na(data$Value_maize2)] <- 0
  }
  data$dominance <- NA
  data$dominance[is.na(data$Value_maize1) | is.na(data$Value_maize2)] <- "Not Expressed"
  data$dominance[data$Value_maize1 >= factor*data$Value_maize2] <- "Maize1"
  data$dominance[data$Value_maize2 >= factor*data$Value_maize1] <- "Maize2"
  data$dominance[data$Value_maize1 <= factor*data$Value_maize2 & data$Value_maize2 <= factor*data$Value_maize1] <- "Neither"

  data <-
    data %>%
    select(Sample, dominance) %>%
    group_by(Sample, dominance) %>%
    count()

  plot <-
    ggplot(data, aes(Sample, n/nPairs*100)) +
    geom_bar(aes(fill=dominance), position="dodge", stat="identity") +
    labs(y = paste0("Percent of gene pairs (n=",nPairs,")"),
         x = "Experiment",
         title = paste0("Which homeolog dominates expression? (factor=",factor,")")
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(plot)
}


## Utility Functions
#--------------------------------------------------------------------------------------------------#
saveMyObjects <- function() {
  rm(list = ls())
  rm(params)
  source("~/git/SubGenomes/R/loadData.R")
  source("~/git/SubGenomes/R/cleanData.R")
  save.image("~/git/SubGenomes/Data/SavedObjects/loadedData.RData")
}

loadMyObjects <- function() {
  load("~/git/SubGenomes/Data/SavedObjects/loadedData.RData", .GlobalEnv)
}
