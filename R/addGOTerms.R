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

goAnnotations.sub1.aggr <-
  goAnnotations.sub1 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())

goAnnotations.sub2.aggr <-
  goAnnotations.sub2 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())


## A simple list of genes in subgenome 1
subgenome.sub1 <-
  subgenome %>%
  filter(subgenome == "sub1") %>%
  select(gene2) %>%
  distinct()
## A simple list of genes in subgenome 2
subgenome.sub2 <-
  subgenome %>%
  filter(subgenome == "sub2") %>%
  select(gene2) %>%
  distinct()
