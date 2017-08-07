library(readr)
library(tidyr)
library(dplyr)
####################################################################################################
## Project: Subgenomes project
## Script purpose: This script should label each gene in the input file as either part of the dominant
##                subgenome or the recessive subgenome.
##
## Input: data input fileName, log10_ks_cutoff
## Output: not really meant to have defined output, since several intermediary objects can be of use
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
log10_ks_cutoff <- params$log10_ks_cutoff

## Import the raw data from the parsed SynMap output
syntelogs.raw <- read_delim(inputDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)

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
# First group rows with same syntenic block and chromosome, then summarize each block's ks values with median/mean/count
# Following the schnable article, syntenic blocks must have 12 genes "The median synonymous substitution rate of all gene pairs in a syntenic block between maize and sorghum can be used to classify syntenic blocks of 12 or more genes unambiguously as orthologous or homoeologous, however" and have a median ks value that discriminates for the alpha duplication event.
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

# Reconstructed ancestral chromosomes
# ! This isn't quite right, since it assumes the single largest homologous.chromosome is the sub1 genome, but sometimes (half the time?)
# ! the sub1 genome is split over multiple chromosomes
subgenome <-
  match(unique(homeologs.chromosome$chr1), homeologs.chromosome$chr1) %>%
  slice(homeologs.chromosome, .) %>%
  select(org_chr1, org_chr2) %>%
  mutate(subgenome="sub1") %>%
  left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2"))

subgenome$subgenome[is.na(subgenome$subgenome)] <- "sub2"

## Add chromosome wide start/stop values, output for greedy algorithm processing in external lanquage
# ! this is a better method than above chunk
# homeologs.chromosomeStopStart <-
#   syntelogs.raw %>%
#   inner_join(homeologs.block, by=c("block", "org_chr1", "org_chr2")) %>%
#   mutate(low = pmin(start1, stop1), high = pmax(start1, stop1)) %>%
#   group_by(org_chr1, org_chr2) %>%
#   summarise(chrLow = min(low,high), chrHigh = max(low,high)) %>%
#   left_join(homeologs.chromosome, ., by=c("org_chr1", "org_chr2")) %>%
#   select(org_chr1, org_chr2, chromosomeGeneCount, chrLow, chrHigh)
# write.table(homeologs.chromosomeStopStart, "outfile.tab", sep="\t")
# subgenome.chromosomes <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\deleteme\\outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
subgenome <-
  subgenome.chromosomes %>%
  select(org_chr1, org_chr2, subgenome) %>%
  left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2"))

# Also, since we know this comparison is 1 sorghum gene = 2 maize genes, lets get the maize genes associated with each other
# using the format
# Maize1  Maize2  Sorghum
# subgenome %>%
#   select(gene1, gene2, subgenome) %>%
#   spread(gene2, subgenome)
subgenome.done <-
  subgenome %>%
  select(gene1, gene2, subgenome) %>%
  distinct() %>%
  group_by(gene1) %>%
  mutate(ind = row_number()) %>%
  spread(subgenome, gene2) %>%
  select(gene1, sub1, sub2) %>%
  group_by(gene1) %>%
  summarise(Maize1=trimws(toString(na.omit(sub1))), Maize2=trimws(toString(na.omit(sub2))))

# Only keep genes where there is a duplicate still on subgenome 1 and subgenome 2
subgenome.homeologs <-
  subgenome.done %>%
  subset(Maize1 != "" & Maize2 != "")




####################################################################################################
## Function greedOpt ->
##        Use a greedy algorithm to sort chromosomes into sugenome 1 and subgenome 2, where
##        subgenome 1 is the larger (by gene count) region that may include multiple org2 chromosomes,
##        but may not overlap when projected to the org1 chromosome.
##
##  homeologs.chromosome -> must be sorted by gene count
####################################################################################################
greedOpt <- function(homeologs.chromosome, syntelogs.raw){
  # match(unique(homeologs.chromosome$chr1), homeologs.chromosome$chr1) %>%
  # slice(homeologs.chromosome, .) %>%
  # select(org_chr1, org_chr2) %>%
  # mutate(subgenome="sub1") %>%
  # left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2"))

  chromosomes <-
    homeologs.chromosome %>%
    select(org_chr2,chromosomeGeneCount) %>%
    distinct()

  homeologs.chromosomeStopStart <-
    syntelogs.raw %>%
    inner_join(homeologs.block, by=c("block", "org_chr1", "org_chr2")) %>%
    mutate(low = pmin(start1, stop1), high = pmax(start1, stop1)) %>%
    group_by(org_chr1, org_chr2) %>%
    summarise(chrLow = min(low,high), chrHigh = max(low,high)) %>%
    left_join(homeologs.chromosome, ., by=c("org_chr1", "org_chr2")) %>%
    select(org_chr1, org_chr2, chromosomeGeneCount, chrLow, chrHigh)

  x <- 1
  while(x < nrow(chromosomes)) {
    currentChromosome <-
      homeologs.chromosomeStopStart %>%
      inner_join(chromosomes[x,], by="org_chr1")

    for (row in 1:nrow(currentChromosome)) {
      count <- currentChromosome[row, "chromosomeGeneCount"]
      low <- currentChromosome[row, "chrLow"]
      high <- currentChromosome[row, "chrHigh"]

    }
    x <- x+1;
  }

  # N           <- ncol(X)
  # weights     <- rep(0L, N)
  # pred        <- 0 * X
  # sum.weights <- 0L
  #
  # while(sum.weights < iter) {
  #
  #   sum.weights   <- sum.weights + 1L
  #   pred          <- (pred + X) * (1L / sum.weights)
  #   errors        <- sqrt(colSums((pred - Y) ^ 2L))
  #   best          <- which.min(errors)
  #   weights[best] <- weights[best] + 1L
  #   pred          <- pred[, best] * sum.weights
  # }
  # return(weights / sum.weights)
}
# greedOpt(homeologs.chromosome)
