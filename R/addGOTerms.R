library(readr)
library(dplyr)
####################################################################################################
## Project: Subgenomes project
## Script purpose:
##
## Input:
## Output:
##
## Date: 2017-07-28
## Author: Jesse R. Walsh
####################################################################################################

goDataFile <- params$goDataFile

# goDataFile <- "C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO from maizecyc.tab"

## Read in GO Annotation data for maize genes
goAnnotations.raw <- read_delim(goDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)

## Parse out citation string
goAnnotations.mutated <-
  goAnnotations.raw %>%
  separate(Citation, c("Publication", "EVCode","TimeStamp","Curator"), sep=":", extra="drop")

## Remove |'s from GO Terms
goAnnotations.mutated$`GO Term` <- gsub(goAnnotations.mutated$`GO Term`, pattern = "\\|", replacement = "")

## GO Analysis -> grouping into useful slices
goAnnotations.sub1 <-
  subgenome %>%
  subset(subgenome == "sub1") %>%
  distinct() %>%
  inner_join(goAnnotations.mutated, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term", EVCode) %>%
  distinct()

goAnnotations.sub2 <-
  subgenome %>%
  subset(subgenome == "sub2") %>%
  distinct() %>%
  inner_join(goAnnotations.mutated, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term", EVCode) %>%
  distinct()
