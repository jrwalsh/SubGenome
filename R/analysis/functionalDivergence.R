####################################################################################################
## Project:
## Script purpose:
##
## Input:
## Output:
## Date: 2018-02-23
## Author: Jesse R. Walsh
####################################################################################################

#--------------------------------------------------------------------------------------------------#
# A function that can catch all alternating gene dominances 2/22/2018
#--------------------------------------------------------------------------------------------------#
# Type 1: subA gene always beats subB gene
# Implies dieing subB gene
# Type 2: subB gene always beats subA gene
# Implies dieing subA gene
# Type 3: Dominance switches, but both genes express
# *Type 4: Dominance switches, but one gene is turned off in each case

library(tidyr)
library(dplyr)
library(ggplot2)
source("./R/scripts/functionalDivergence_func.R")

maize.expression.sample.avg
expressedPairs
retained.duplicates <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "")

which(retained.duplicates=="Zm00001d033174", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034368", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034429", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034517", arr.ind = TRUE)
index <- 120
test.cutoff <- 2
plotRetainedDuplicatesExpression(retained.duplicates, maize.expression.sample.avg, index)
plotRetainedDuplicatesAbundance(retained.duplicates, maize.protein.abundance.sample.avg, index)
index
index <- index + 1

getType3 <- function() {

}

getRetainedDuplicatedStats <- function(retained.duplicates, maize.expression.sample.avg) {
  retained.duplicates.stats <- cbind(retained.duplicates$Maize1, retained.duplicates$Maize2)
  retained.duplicates.stats <- unique(retained.duplicates.stats)
  retained.duplicates.stats

  test.expression <- maize.expression.sample.avg
  test.expression$FPKM_avg[is.na(test.expression$FPKM_avg)] <- 0

  test.expression %>%
    group_by(geneID) %>%
    summarise(min = min(FPKM_avg), max = max(FPKM_avg), avg = mean(FPKM_avg), var = var(FPKM_avg))
}

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
