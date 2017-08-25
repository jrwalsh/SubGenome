library(topGO)
library(dplyr)
library(tidyr)
####################################################################################################
## Project:
## Script purpose:
##
## Input: goAnnotations.mutated from the addGOTerms.R script
## Output:
##
## Date: 2017-08-18
## Author: Jesse R. Walsh
####################################################################################################
## GOBPTerm, GOMFTerm and GOCCTerm are part of the topGO package
MFterms <- ls(GOMFTerm)
# head(MFterms)


# goDataFile <- params$goDataFile
# goDataFile <- "C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO from maizecyc.tab"
## Read in GO Annotation data for maize genes
# if (!file.exists(goDataFile)) {
#   goAnnotations.raw <- read_delim(goDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
# } else {
#   stop()
# }


## Create a topGO custom geneID2GO annotation object, mapping gene ids to go terms
geneID2GO.temp <-
  goAnnotations.mutated %>%
  select(geneID, `GO Term`) %>%
  group_by(geneID) %>%
  summarise(GO = paste(`GO Term`, collapse = ", "))
geneID2GO <- strsplit(geneID2GO.temp$GO, ", ")
names(geneID2GO) <- geneID2GO.temp$geneID

geneNames <- names(geneID2GO)
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
# myInterestingGenes <- geneNames
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

topGOdata <- new("topGOdata", description = "getGO", ontology="MF", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)



## Which annotations are MF?
goAnnotations.sub1.temp <- goAnnotations.sub1
goAnnotations.sub1.temp$`GO Term` <- gsub(goAnnotations.sub1.temp$`GO Term`, pattern = "\\|", replacement = "")
goAnnotations.sub1.temp %>%
  subset(`GO Term` %in% MFterms)



# geneID2GO.temp$GO <- as.character(strsplit(geneID2GO.temp$GO, ","))
# geneID2GO <- as.character(geneID2GO.temp$GO)
# names(geneID2GO) <- geneID2GO.temp$gene2

# geneID2GO <-
#   goAnnotations.sub1 %>%
#   nest("GO Term", .key="GO") %>%
#   select(gene2, GO)
# geneID2GO <- as.character(geneID2GO$GO)
# setNames(geneID2GO) <- geneID2GO$gene2
