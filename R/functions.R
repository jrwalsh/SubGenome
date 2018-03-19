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
## maize.data = long-form expression/protein data
## Factor = how many times greater expression has to be to be considered "dominant" -> (i.e. is A >= factor*B)
## omitNA = if FALSE, NA values in the expression dataset for either homeolog result in a "Not Expressed",
##          if TRUE, NA values are converted to zero's, so any expression for one homeolog will automatically dominate the other and NA in both results in a "Neither"
##          this setting does not stop NAs in the case where the gene pair doesn't exist in the expression set
dataGenePairDominanceByExperiment <- function(maize.data, homeologs.pairs, factor, includeNA) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  data <-
    genePairs %>%
    inner_join(maize.data, by=c("Maize1"="geneID")) %>%
    inner_join(maize.data, by=c("Maize2"="geneID", "sample"="sample"))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"

  if (includeNA) {
    data$Value_maize1[is.na(data$Value_maize1)] <- 0
    data$Value_maize2[is.na(data$Value_maize2)] <- 0
  }
  data$dominance <- NA
  data$dominance[is.na(data$Value_maize1) | is.na(data$Value_maize2)] <- "Not Expressed"
  data$dominance[!is.na(data$Value_maize1) & !is.na(data$Value_maize2) & data$Value_maize1 >= factor*data$Value_maize2] <- "Maize1"
  data$dominance[!is.na(data$Value_maize1) & !is.na(data$Value_maize2) & data$Value_maize2 >= factor*data$Value_maize1] <- "Maize2"
  data$dominance[!is.na(data$Value_maize1) & !is.na(data$Value_maize2) & data$Value_maize1 <= factor*data$Value_maize2 & data$Value_maize2 <= factor*data$Value_maize1] <- "Neither"

  data <-
    data %>%
    select(sample, dominance) %>%
    group_by(sample, dominance) %>%
    count()

  return(data)
}

graphGenePairDominanceByExperiment <- function(maize.data, homeologs.pairs, factor, includeNA) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  data <-
    genePairs %>%
    inner_join(maize.data, by=c("Maize1"="geneID")) %>%
    inner_join(maize.data, by=c("Maize2"="geneID", "sample"="sample"))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"

  ## We only count gene pairs where both genes were at least tested in the data set
  nPairs <-
    data %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    nrow()

  data <- dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, factor, includeNA)
  plot <-
    ggplot(data, aes(sample, n/nPairs*100)) +
    geom_bar(aes(fill=dominance), position="dodge", stat="identity") +
    labs(y = paste0("Percent of gene pairs (n=",nPairs,")"),
         x = "Sample",
         title = paste0("Which homeolog has greater expression? (factor=",factor,")")
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(plot)
}

## Deliberately ignore all pairs that have a "Not Expressed" designation
## Full_join needed to prevent missing data for one Sample/Factor from hiding data from others, need to make this missing data "0"
graphGenePairDominanceByExperimentAcrossFactors <- function(maize.data, homeologs.pairs, includeNA) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  data <-
    genePairs %>%
    inner_join(maize.data, by=c("Maize1"="geneID")) %>%
    inner_join(maize.data, by=c("Maize2"="geneID", "sample"="sample"))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"

  ## We only count gene pairs where both genes were at least tested in the data set
  nPairs <-
    data %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    nrow()

  data <-
    dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 1, includeNA) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 2, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 3, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 4, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 5, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 6, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 7, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 8, includeNA), by = c("sample", "dominance")) %>%
    full_join(dataGenePairDominanceByExperiment(maize.data, homeologs.pairs, 9, includeNA), by = c("sample", "dominance"))
  colnames(data)[c(-1,-2)] <- c("1","2","3","4","5","6","7","8","9")#,"10")
  data <-
    data %>%
    gather(Factor, n, c(-1,-2))
  data[is.na(data)] <- 0
  data <- data[!data$dominance=="Not Expressed",]

  plot <-
    ggplot(data, aes(x=Factor, y=n/nPairs*100, group=dominance, color=dominance)) +
    geom_line() + facet_wrap(~sample, nrow=6) +
    labs(
      title = paste0("Percent of homeologs from each subgenome\nwhich are overexpressed compared to their pair\nby a given factor.")
    )
  return(plot)
}

## geneTranscript.counts = # of isoforms per gene
## Factor = how many times greater expression has to be to be considered "dominant" -> (i.e. is A >= factor*B)
graphGenePairDominanceByIsoForms <- function(geneTranscript.counts, homeologs.pairs, factor) {
  data <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct() %>%
    inner_join(geneTranscript.counts, by=c("Maize1"="gene")) %>%
    inner_join(geneTranscript.counts, by=c("Maize2"="gene"))
  names(data)[3] <- "Value_maize1"
  names(data)[4] <- "Value_maize2"
  data$dominance <- NA
  data$dominance[is.na(data$Value_maize1) | is.na(data$Value_maize2)] <- "Not Data"
  data$dominance[data$Value_maize1 >= factor*data$Value_maize2] <- "Maize1"
  data$dominance[data$Value_maize2 >= factor*data$Value_maize1] <- "Maize2"
  data$dominance[data$Value_maize1 <= factor*data$Value_maize2 & data$Value_maize2 <= factor*data$Value_maize1] <- "Neither"

  data <-
    data %>%
    count(dominance)

  plot <-
    ggplot(data, aes(dominance, n)) +
    geom_bar(aes(fill=dominance), position="dodge", stat="identity") +
    labs(y = paste0("Number of pairs dominated (counts)"),
         x = "Dominance",
         title = paste0("Which homeolog has greater isoform counts? (factor=",factor,")")
    ) +
    geom_text(aes(label=n), vjust=-0.3, size=3.5)

  return(plot)
}


homeologFoldChanges <- function(maize.data, homeologs.pairs) {
  genePairs <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>%
    select(Maize1, Maize2) %>%
    distinct()

  data <-
    genePairs %>%
    inner_join(maize.data, by=c("Maize1"="geneID")) %>%
    inner_join(maize.data, by=c("Maize2"="geneID", "sample"="sample"))
  names(data)[4] <- "Value_maize1"
  names(data)[5] <- "Value_maize2"

  ## Remove NAs and calculate foldchange
  data <- subset(data, !is.na(data$Value_maize1) & !is.na(data$Value_maize2))
  data$foldChange <- log2(data$Value_maize1) - log2(data$Value_maize2)
  data <- arrange(data, desc(foldChange))
  return(data)
}


## Utility Functions
#--------------------------------------------------------------------------------------------------#
saveMyObjects <- function() {
  rm(list = ls())
  # rm(params) # For knitr params
  source("./R/loadData.R")
  source("./R/cleanData.R")
  save.image("./data/cleanedData.RData")
}

loadMyObjects <- function() {
  load("./data/cleanedData.RData", .GlobalEnv)
}

getGenePageURL <- function(geneID) {
  base <- "https://maizegdb.org/gene_center/gene/"
  return(paste0(base, geneID))
}

getGBrowsePageURL <- function(geneID) {
  base <- "https://www.maizegdb.org/gbrowse/maize_v4/?name="
  return(paste0(base, geneID))
}


## Example: apply(gene.foldchange[,1], 2, function(x) { wrapInHREF(getGenePageURL(x), x) })
wrapInHREF <- function(url, text) {
  return(paste0("<a href=\"", url, "\">", text, "</a>"))
}
