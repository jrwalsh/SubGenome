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

goAnnotations <- read_delim("./Data/GO from maizecyc.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

goAnnotations %>%
  select(V4_ID, `MaizeCyc2.2 Accession-1`, `GO Term`) %>%
  group_by(V4_ID) %>%
  merge(homeologs.genes, by.x = "MaizeCyc2.2 Accession-1", by.y = "gene2") %>%
  select(`MaizeCyc2.2 Accession-1`, `GO Term`) %>%
  distinct() -> homeologs.go

subgenome %>%
  filter(subgenome == "sub1") %>%
  select(gene2) %>%
  distinct() -> subgenome.sub1

subgenome %>%
  filter(subgenome == "sub2") %>%
  select(gene2) %>%
  distinct() -> subgenome.sub2

goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub1$gene2,] %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(count(`GO Term`))

goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub2$gene2,]
