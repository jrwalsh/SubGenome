library(readr)
library(tidyr)
library(dplyr)
####################################################################################################
## Project: Subgenomes project
## Script purpose: This script should label each gene in the input file as either part of the dominant
##                subgenome or the recessive subgenome.
##
## Input: data input fileName, log10_ks_cutoff
## Output:
##        syntelogs.raw         ->
##        syntelogs.mutated     ->
##        homeologs.genes       ->
##        homeologs.block       ->
##        homeologs.chromosome  ->
##        subgenome             ->
##        ks_cutoff             ->
## Date: 2017-07-26
## Author: Jesse R. Walsh
####################################################################################################

inputDataFile <- params$inputDataFile
subGenomeFile <- params$subGenomeFile
log10_ks_cutoff <- params$log10_ks_cutoff

# inputDataFile <- "C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\sorghum_v1_vs_maize_v1.tab"
# subGenomeFile <- "C:\\Users\\Jesse\\Dropbox (Personal)\\deleteme\\outfile_processed_byHand.tab"

## Import the raw data from the parsed SynMap output
syntelogs.raw <- read_delim(inputDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)

## Import "True" sugenomes as reported by Schnable 2011
subgenome.truth <- setNames(data.frame(
  c(1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,7,7,8,8,8,9,9,9,10,10,10),
  c(1,5,9,7,2,3,8,5,4,4,2,2,10,1,6,10,4,1,10,3,6,10,8,5,9,6),
  c("sub1","sub2","sub2","sub1","sub2","sub1","sub2","sub1","sub2",
    "sub1","sub2","sub1","sub2","sub1","sub1","sub1","sub2","sub1",
    "sub1","sub2","sub1","sub1","sub2","sub1","sub1","sub2"),
  stringsAsFactors=FALSE), c("chr1","chr2","subgenome"))

## Add median and geneCount values (aggregated by block/org_chr1/org_chr2) to each row.  Calculate gene sizes.  Add chromosome ids as numbers.
syntelogs.mutated <-
  syntelogs.raw %>%
  group_by(block, org_chr1, org_chr2) %>%
  summarise(median_ks=median(ks, na.rm=TRUE), blockGeneCount=n()) %>%
  left_join(syntelogs.raw, ., by = c("block", "org_chr1", "org_chr2")) %>%
  mutate(gene_length1=abs(stop1-start1)) %>%
  mutate(gene_length2=abs(stop2-start2)) %>%
  mutate(chr1=as.numeric(regmatches(org_chr1, regexpr("\\d*$",org_chr1)))) %>%
  mutate(chr2=as.numeric(regmatches(org_chr2, regexpr("\\d*$",org_chr2))))

## Normal scale of the log10_ks_cutoff
ks_cutoff <- 10^log10_ks_cutoff

## Find which sets of chromosomes with syntelogs should be in group 1 or group 2, where group 1 has larger syntenic blocks
## First group rows with same syntenic block and chromosome, then summarize each block's ks values with median/mean/count
## Following the schnable article, syntenic blocks must have 12 genes "The median synonymous substitution rate of all gene
## pairs in a syntenic block between maize and sorghum can be used to classify syntenic blocks of 12 or more genes unambiguously
## as orthologous or homoeologous, however" and have a median ks value that discriminates for the alpha duplication event.
homeologs.block <-
  syntelogs.mutated %>%
  select(block, chr1, org_chr1, chr2, org_chr2, median_ks, blockGeneCount) %>%
  distinct() %>%
  filter(blockGeneCount >= 12 & median_ks <= ks_cutoff)

homeologs.genes <-
  homeologs.block %>%
  inner_join(syntelogs.mutated, by=c("block", "org_chr1", "org_chr2")) %>%
  select(org_chr1, org_chr2, gene1, gene2) %>%
  distinct()

# Determine which chromosome has largest syntenic block, the syntelogs on this chromosome are subgenome 1
homeologs.chromosome <-
  homeologs.block %>%
  select(block, org_chr1, chr1, org_chr2, chr2, blockGeneCount) %>%
  group_by(org_chr1, chr1, org_chr2, chr2) %>%
  summarise(chromosomeGeneCount=sum(blockGeneCount)) %>%
  arrange(chr1, desc(chromosomeGeneCount)) %>%
  ungroup()

## Add chromosome wide start/stop values, output for greedy algorithm processing in external lanquage
homeologs.chromosomeStopStart <-
  syntelogs.raw %>%
  inner_join(homeologs.block, by=c("block", "org_chr1", "org_chr2")) %>%
  mutate(low = pmin(start1, stop1), high = pmax(start1, stop1)) %>%
  group_by(org_chr1, org_chr2) %>%
  summarise(chrLow = min(low,high), chrHigh = max(low,high)) %>%
  left_join(homeologs.chromosome, ., by=c("org_chr1", "org_chr2")) %>%
  select(org_chr1, org_chr2, chromosomeGeneCount, chrLow, chrHigh)

## Currently need to perform the greedy sorting by hand... need to automate this step.
# write.table(homeologs.chromosomeStopStart, "outfile.tab", sep="\t")
subgenome.chromosomes <- read_delim(subGenomeFile, "\t", escape_double = FALSE, trim_ws = TRUE)

subgenome <-
  subgenome.chromosomes %>%
  select(org_chr1, org_chr2, subgenome) %>%
  left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2"))

## Also, since we know this comparison is 1 sorghum gene = 2 maize genes, lets group the homeologous pairs
homeologs.pairs <-
  subgenome %>%
  select(gene1, gene2, subgenome) %>%
  distinct() %>%
  group_by(gene1) %>%
  mutate(ind = row_number()) %>%
  spread(subgenome, gene2) %>%
  select(gene1, sub1, sub2) %>%
  group_by(gene1) %>%
  summarise(Maize1=trimws(toString(na.omit(sub1))), Maize2=trimws(toString(na.omit(sub2))))

## Only keep genes where there is a duplicate still on subgenome 1 and subgenome 2
subgenome.homeologs <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "")
