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

index <- 360
test.cutoff <- 2
getType3(retained.duplicates, maize.expression.sample.avg, index)
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

# annotate("text", x = Inf, y = 10, label = test.cor) +
# annotate("rect", xmin = .5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Green") +
# annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
# annotate("rect", xmin = 3.5, xmax = 9.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
# annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
# annotate("rect", xmin = 10.5, xmax = 11.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
# annotate("rect", xmin = 11.5, xmax = 15.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "DarkGreen") +
# annotate("rect", xmin = 15.5, xmax = 16.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
# annotate("rect", xmin = 16.5, xmax = 21.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Brown") +
# annotate("rect", xmin = 21.5, xmax = 22.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
# annotate("rect", xmin = 22.5, xmax = 23.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Green") +

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
