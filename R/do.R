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
