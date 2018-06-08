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

# retained.duplicates <-
#   homeologs.pairs %>%
#   subset(Maize1 != "" & Maize2 != "")
#
# data <-
#   maize.walley.v4mapped.expression %>%
#   select(geneID, sample, FPKM_avg) %>%
#   subset(geneID %in% retained.duplicates$Maize1 | geneID %in% retained.duplicates$Maize2)
#
# data$FPKM_avg[is.na(data$FPKM_avg)] <- 0

# expressedPairs

which(retained.duplicates=="Zm00001d033174", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034368", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034429", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d034517", arr.ind = TRUE)
index <- 120
test.cutoff <- 2
plotRetainedDuplicatesExpression(retained.duplicates, data, index)
plotRetainedDuplicatesAbundance(retained.duplicates, maize.walley.abundance.v4, index)
index
index <- index + 1

gene.pair.types <- getTypeSort(homeologs.pairs, maize.walley.v4mapped.expression, 2, .95)
summary(as.factor(getTypeSort(homeologs.pairs, maize.walley.v4mapped.expression, 2, .95)$type))
protein.pair.types <- getTypeSortProteins(homeologs.pairs, maize.walley.abundance.v4, 2, .95)
summary(as.factor(getTypeSortProteins(homeologs.pairs, maize.walley.abundance.v4, 2, .95)$type))
intersect(gene.pair.types, protein.pair.types)
intersect(gene.pair.types[!gene.pair.types$type %in% c("ambiguous"),], protein.pair.types[!protein.pair.types$type %in% c("ambiguous"),])
#dead expression and still dead abundance (17/297)
protein.pair.types[protein.pair.types$Maize1 %in% gene.pair.types$Maize1[gene.pair.types$type %in% c("1dead","2dead")],] %>% subset(type %in% c("1dead","2dead"))

index <- which(retained.duplicates=="Zm00001d034896", arr.ind = TRUE)[1]
plotRetainedDuplicatesExpression(retained.duplicates, data, index)

## Direct method, but shouldn't we avg replicates first?
# data <- maize.walley.v4mapped.expression.replicate
# data$rowMax <- apply(data[,-1], 1, function(x) max(x, na.rm = TRUE))
# data$rowMax[is.infinite(data$rowMax)] <- 0
# data <-
#   data %>%
#   select(geneID, rowMax) %>%
#   subset(geneID %in% retained.duplicates$Maize1 | geneID %in% retained.duplicates$Maize2)
# data %>%
#   subset(rowMax < 1)

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
test.expression.temp <- maize.walley.v4mapped.expression
test.expression.temp$FPKM_avg[is.na(test.expression.temp$FPKM_avg)] <- 0
# Give me the very highly correlated pairs (306/3571 > .95, 597/3571 > .9, 916/3571 > .85)
  retained.duplicates %>%
  inner_join(test.expression.temp, by=c("Maize1"="geneID")) %>%
  rename(FPKM_avg1=FPKM_avg) %>%
  inner_join(test.expression.temp, by=c("Maize2"="geneID", "sample"="sample")) %>%
  rename(FPKM_avg2=FPKM_avg) %>%
  select(Maize1,Maize2,sample,FPKM_avg1,FPKM_avg2) %>%
  group_by(Maize1,Maize2) %>%
  summarise(correlation=cor(FPKM_avg1, FPKM_avg2), mean1=mean(FPKM_avg1), mean2=mean(FPKM_avg2)) %>%
  subset(correlation > .95 & (mean1 > 5 | mean2 > 5))

# Give me the dead genes (310/3571 have 1 or both < 1, 500/3571 < 2, 625/3571 < 3, 324/3571 one is > 10 and other < 3 + 129/3571 where both < 3)
  retained.duplicates %>%
  inner_join(test.expression.temp, by=c("Maize1"="geneID")) %>%
  rename(FPKM_avg1=FPKM_avg) %>%
  inner_join(test.expression.temp, by=c("Maize2"="geneID", "sample"="sample")) %>%
  rename(FPKM_avg2=FPKM_avg) %>%
  select(Maize1,Maize2,sample,FPKM_avg1,FPKM_avg2) %>%
  group_by(Maize1,Maize2) %>%
  summarise(max1=max(FPKM_avg1), max2=max(FPKM_avg2), sd1=sd(FPKM_avg1), sd2=sd(FPKM_avg2), mean1=mean(FPKM_avg1), mean2=mean(FPKM_avg2)) %>%
  subset(max2 < 1 & max1 < 1)

which(retained.duplicates=="Zm00001d001887", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d001961", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d002011", arr.ind = TRUE)
index <- 2928
index <- 2908
index <- 2900
plotRetainedDuplicatesExpression(retained.duplicates, test.expression.temp, index)
plotRetainedDuplicatesAbundance(retained.duplicates, maize.protein.abundance.sample.avg, index)

# Give me BED genes (1464/3571 with 2way foldchange of 2, 884/3571 with 2way foldchange of 4, 687/3571 with 2way foldchange of 10)
# Give me BED genes "avg must be above 3" (1,070/3571 with 2way foldchange of 2, 434/3571 with 2way foldchange of 4, 227/3571 with 2way foldchange of 10)
myTemp <-
  retained.duplicates %>%
  inner_join(test.expression.temp, by=c("Maize1"="geneID")) %>%
  rename(FPKM_avg1=FPKM_avg) %>%
  inner_join(test.expression.temp, by=c("Maize2"="geneID", "sample"="sample")) %>%
  rename(FPKM_avg2=FPKM_avg) %>%
  select(Maize1,Maize2,sample,FPKM_avg1,FPKM_avg2) %>%
  mutate(diff=FPKM_avg1-FPKM_avg2, fold1=FPKM_avg1/FPKM_avg2, fold2=FPKM_avg2/FPKM_avg1) %>%
  subset(FPKM_avg1 > 3 | FPKM_avg2 > 3) %>%
  subset(diff > 1 | diff < -1) %>%
  group_by(Maize1,Maize2) %>%
  summarise(dom1=sum(diff>0), dom2=sum(diff<0), maxFoldChange1=max(fold1), maxFoldChange2=max(fold2))

myTemp$maxFoldChange1[is.infinite(myTemp$maxFoldChange1)] <- 100
myTemp$maxFoldChange2[is.infinite(myTemp$maxFoldChange2)] <- 100

  myTemp %>%
  subset((dom1 == 0 | dom2 == 0) & (maxFoldChange1 > 2 | maxFoldChange2 > 2))

which(retained.duplicates=="Zm00001d001839", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d001862", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d001989", arr.ind = TRUE)
which(retained.duplicates=="Zm00001d002721", arr.ind = TRUE)
index <- 2943
index <- 2936
index <- 2906
index <- 2786
plotRetainedDuplicatesExpression(retained.duplicates, test.expression.temp, index)
plotRetainedDuplicatesAbundance(retained.duplicates, maize.protein.abundance.sample.avg, index)

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
