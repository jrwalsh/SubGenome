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

# Read in GO Annotation data for maize genes
goAnnotations <- read_delim("./Data/GO from maizecyc.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

# Merge v3 homeolog gene IDs to the GO term data, each line has a unique Gene -> GO term pair, where genes and go terms can be reused
# This list is only genes that have homeologs and have GO annotation
# !! This is an unnessary step, but might be useful later...
# goAnnotations %>%
#   select(V4_ID, `MaizeCyc2.2 Accession-1`, `GO Term`) %>%
#   group_by(V4_ID) %>%
#   merge(homeologs.genes, by.x = "MaizeCyc2.2 Accession-1", by.y = "gene2") %>%
#   select(`MaizeCyc2.2 Accession-1`, `GO Term`) %>%
#   distinct() -> homeologs.go

# A simple list of (v3) genes in subgenome 1
subgenome %>%
  filter(subgenome == "sub1") %>%
  select(gene2) %>%
  distinct() -> subgenome.sub1

# A simple list of (v3) genes in subgenome 2
subgenome %>%
  filter(subgenome == "sub2") %>%
  select(gene2) %>%
  distinct() -> subgenome.sub2

# Rows from go annotation file annotating a gene from sub1 or sub2
goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub1$gene2,]
goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub2$gene2,]

#
goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub1$gene2,] %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(count(`GO Term`))
