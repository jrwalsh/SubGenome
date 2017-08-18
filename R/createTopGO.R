library(topGO)
library(dplyr)
library(tidyr)
####################################################################################################
## Project:
## Script purpose:
##
## Input:
## Output:
## Date: 2017-08-18
## Author: Jesse R. Walsh
####################################################################################################

## GOBPTerm, GOMFTerm and GOCCTerm are part of the topGO package
MFterms <- ls(GOMFTerm)
head(MFterms)

## Create a topGO custom geneID2GO annotation object, mapping gene ids to go terms
geneID2GO.temp <-
  goAnnotations.sub1 %>%
  group_by(gene2) %>%
  summarise(GO = paste(`GO Term`, collapse = ", "))
geneID2GO.temp$GO <- as.character(strsplit(geneID2GO.temp$GO, ","))
geneID2GO <- as.character(geneID2GO.temp$GO)
names(geneID2GO) <- geneID2GO.temp$gene2


## Which annotations are MF?
goAnnotations.sub1.temp <- goAnnotations.sub1
goAnnotations.sub1.temp$`GO Term` <- gsub(goAnnotations.sub1.temp$`GO Term`, pattern = "\\|", replacement = "")
goAnnotations.sub1.temp %>%
  subset(`GO Term` %in% MFterms)






# geneID2GO <-
#   goAnnotations.sub1 %>%
#   nest("GO Term", .key="GO") %>%
#   select(gene2, GO)
# geneID2GO <- as.character(geneID2GO$GO)
# setNames(geneID2GO) <- geneID2GO$gene2
