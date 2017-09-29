library(tidyr)
library(dplyr)
####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Create and combine objects to use for analysis and plots.
##
## Input:
## Output:
## Date: 2017-09-14
## Author: Jesse R. Walsh
####################################################################################################
#==================================================================================================#
## parseSubgenomes.R
#--------------------------------------------------------------------------------------------------#
source("~/git/SubGenomes/R/parseSubgenomes.R")

#==================================================================================================#
## addGOTerms.R
#--------------------------------------------------------------------------------------------------#
## GO Analysis -> grouping into useful slices
goAnnotations.sub1 <-
  subgenome %>%
  subset(subgenome == "sub1") %>%
  distinct() %>%
  inner_join(go.maize.clean, by = c("gene2" = "geneID")) %>%
  select(gene2, goTerm, evCode, type) %>%
  distinct()

goAnnotations.sub2 <-
  subgenome %>%
  subset(subgenome == "sub2") %>%
  distinct() %>%
  inner_join(go.maize.clean, by = c("gene2" = "geneID")) %>%
  select(gene2, goTerm, evCode, type) %>%
  distinct()

goAnnotations.sub1.aggr <-
  goAnnotations.sub1 %>%
  select(goTerm) %>%
  group_by(goTerm) %>%
  summarise(n=n())

goAnnotations.sub2.aggr <-
  goAnnotations.sub2 %>%
  select(goTerm) %>%
  group_by(goTerm) %>%
  summarise(n=n())

#==================================================================================================#
## parseSorghumGO.R
#--------------------------------------------------------------------------------------------------#
maizeWithSorghumGO <-
  syntelogs.mutated %>%
  select(gene1, gene2)

maizeWithSorghumGO <-
  maizeWithSorghumGO %>%
  left_join(go.sorghum.clean, by=c("gene1"="sorghumID"))

maizeWithSorghumGO <-
  maizeWithSorghumGO %>%
  filter(!is.na(goTerm)) %>%
  select(gene2, goTerm) %>%
  distinct()

#==================================================================================================#
## maize.expression.all
#--------------------------------------------------------------------------------------------------#
maize.expression.all <- maize.expression.clean

## Convert to v4 ids
maize.expression.all <-
  maize.expression.all %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("tracking_id" = "v3_id")) %>%
  rename(geneID=v4_id)

## Reorder columns
maize.expression.all <- maize.expression.all[,c(70,2:69)]

## Merge duplicate rows by adding FPKM values together
maize.expression.all <-
  maize.expression.all %>%
  group_by(geneID) %>%
  summarise_all(funs(sum))

#==================================================================================================#
## maize.expression.sample.avg
#--------------------------------------------------------------------------------------------------#
maize.expression.sample.avg <- maize.expression.clean

## Convert to v4 ids
maize.expression.sample.avg <-
  maize.expression.sample.avg %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("tracking_id" = "v3_id")) %>%
  rename(geneID=v4_id)

## Reorder columns
maize.expression.sample.avg <- maize.expression.sample.avg[,c(70,2:69)]

## convert to sample names, merge replicates in each sample using mean, and output in long form
maize.expression.sample.avg <-
  maize.expression.sample.avg %>%
  gather("tracking_id", "FPKM",-1) %>%
  left_join(experiment.map, by=c("tracking_id"="tracking_id")) %>%
  select(geneID, Sample, FPKM) %>%
  group_by(geneID, Sample) %>%
  summarise(FPKM_avg=mean(FPKM, na.rm=TRUE)) %>%
  arrange(geneID)

## When all replicates have NA, mean returns NaN.  Convert it back to NA.
maize.expression.sample.avg$FPKM_avg[is.nan(maize.expression.sample.avg$FPKM_avg)] <- NA

#==================================================================================================#
## maize.protein.abundance.sample.avg
#--------------------------------------------------------------------------------------------------#
maize.protein.abundance.sample.avg <- maize.protein.abundance.clean

## Convert to v4 ids
maize.protein.abundance.sample.avg <-
  maize.protein.abundance.sample.avg %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("v3_id" = "v3_id")) %>%
  rename(geneID=v4_id)

## Reorder columns
maize.protein.abundance.sample.avg <- maize.protein.abundance.sample.avg[,c(length(maize.protein.abundance.sample.avg),3:length(maize.protein.abundance.sample.avg)-1)]

## Reduce columns to only the ones comparable to expression data
colsToKeep <- colnames(maize.protein.abundance.sample.avg) %in% c("geneID",experiment.map.proteins$Replicate[!is.na(experiment.map.proteins$ExpressionSampleName)])
maize.protein.abundance.sample.avg <- maize.protein.abundance.sample.avg[,colsToKeep]

## convert to sample names, merge replicates in each sample using mean, and output in long form
maize.protein.abundance.sample.avg <-
  maize.protein.abundance.sample.avg %>%
  gather("Replicate", "dNSAF",-1) %>%
  left_join(experiment.map.proteins, by=c("Replicate"="Replicate")) %>%
  select(geneID, Sample, dNSAF) %>%
  group_by(geneID, Sample) %>%
  summarise(dNSAF_avg=mean(dNSAF, na.rm=TRUE)) %>%
  arrange(geneID)

## When all replicates have NA, mean returns NaN.  Convert it back to NA.
maize.protein.abundance.sample.avg$dNSAF_avg[is.nan(maize.protein.abundance.sample.avg$dNSAF_avg)] <- NA

#==================================================================================================#
## parseExpressionData.R
#--------------------------------------------------------------------------------------------------#
expressedGenes <- data.frame(ID=maize.expression.clean[,1], Means=rowMeans(maize.expression.clean[,-1], na.rm = TRUE))

## Remove rows with NA
expressedGenes <-
  expressedGenes %>%
  rename(geneID = tracking_id, FPKM_mean = Means) %>%
  subset(!is.na(FPKM_mean))

## Convert to v4 ids
expressedGenes <-
  expressedGenes %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, FPKM_mean) %>%
  rename(geneID=v4_id)

## Converting from v3 to v4 will give duplicate values (when gene models are merged, etc.), assume they are they same length so we can add their FPKM together
expressedGenes <-
  expressedGenes %>%
  group_by(geneID) %>%
  summarise(FPKM_mean=sum(FPKM_mean))

## Attach log2(FPKM_mean) values to homeologous pairs
expressedPairs <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  inner_join(expressedGenes, by=c("Maize1"="geneID")) %>%
  rename(FPKM_mean1=FPKM_mean) %>%
  inner_join(expressedGenes, by=c("Maize2"="geneID")) %>%
  rename(FPKM_mean2=FPKM_mean)

#==================================================================================================#
## analyzeGODiffs.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/analyzeGODiffs.R")

#==================================================================================================#
## createTopGO.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/createTopGO.R")

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
