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
syntelogs.raw %>%
  group_by(block, org_chr1, org_chr2) %>%
  summarise(median_ks=median(ks, na.rm=TRUE), blockGeneCount=n()) %>%
  left_join(syntelogs.raw, ., by = c("block", "org_chr1", "org_chr2")) %>%
  mutate(gene_length1=abs(stop1-start1)) %>%
  mutate(gene_length2=abs(stop2-start2)) %>%
  mutate(chr1=as.numeric(regmatches(org_chr1, regexpr("\\d*$",org_chr1)))) %>%
  mutate(chr2=as.numeric(regmatches(org_chr2, regexpr("\\d*$",org_chr2)))) -> syntelogs.mutated

## Normal scale of the log10_ks_cutoff
ks_cutoff <- 10^log10_ks_cutoff

## Find which sets of chromosomes with syntelogs should be in group 1 or group 2, where group 1 has larger syntenic blocks
# First group rows with same syntenic block and chromosome, then summarize each block's ks values with median/mean/count
# Following the schnable article, syntenic blocks must have 12 genes "The median synonymous substitution rate of all gene pairs in a syntenic block between maize and sorghum can be used to classify syntenic blocks of 12 or more genes unambiguously as orthologous or homoeologous, however" and have a median ks value that discriminates for the alpha duplication event.
syntelogs.mutated %>%
  select(block, chr1, org_chr1, chr2, org_chr2, median_ks, blockGeneCount) %>%
  distinct() %>%
  filter(blockGeneCount >= 12 & median_ks <= ks_cutoff) -> homeologs.block

homeologs.block %>%
  inner_join(syntelogs.mutated, by=c("block", "org_chr1", "org_chr2")) %>%
  select(org_chr1, org_chr2, gene1, gene2) %>%
  distinct() -> homeologs.genes

# Determine which chromosome has largest syntenic block, the syntelogs on this chromosome are subgenome 1
homeologs.block %>%
  select(block, org_chr1, chr1, org_chr2, chr2, blockGeneCount) %>%
  group_by(org_chr1, chr1, org_chr2, chr2) %>%
  summarise(chromosomeGeneCount=sum(blockGeneCount)) %>%
  arrange(chr1, desc(chromosomeGeneCount)) %>%
  ungroup() -> homeologs.chromosome

# Reconstructed ancestral chromosomes
# ! This isn't quite right, since it assumes the single largest homologous.chromosome is the sub1 genome, but sometimes (half the time?)
# ! the sub1 genome is split over multiple chromosomes
match(unique(homeologs.chromosome$chr1), homeologs.chromosome$chr1) %>%
  slice(homeologs.chromosome, .) %>%
  select(org_chr1, org_chr2) %>%
  mutate(subgenome="sub1") %>%
  left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2")) -> subgenome
subgenome$subgenome[is.na(subgenome$subgenome)] <- "sub2"

# Also, since we know this comparison is 1 sorghum gene = 2 maize genes, lets get the maize genes associated with each other
# subgenome %>%
#   select(gene1, gene2, subgenome) %>%
#   spread(gene2, subgenome)
