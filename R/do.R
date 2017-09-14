library(tidyr)
library(dplyr)
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
## parseExpressionData.R
#--------------------------------------------------------------------------------------------------#
expressedGenes <- data.frame(ID=maize.expression.clean[,1], Means=rowMeans(maize.expression.clean[,-1], na.rm = TRUE))

## Remove rows with NA
expressedGenes <-
  expressedGenes %>%
  mutate(FPKM_mean = Means) %>%
  select(tracking_id, FPKM_mean) %>%
  rename(geneID = tracking_id) %>%
  subset(!is.na(FPKM_mean))

## Convert to v4 ids (this is a bad conversion table, only has CornCyc Genes in it.... get a better one!)
geneIDConvertTable <- maize.genes.v3_to_v4_map.clean

expressedGenes <-
  expressedGenes %>%
  inner_join(geneIDConvertTable, by=c("geneID" = "v3_id")) %>%
  select(v4_id, FPKM_mean) %>%
  rename(geneID=v4_id)

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
