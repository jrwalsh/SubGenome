library(GenomicFeatures)
library(tidyr)
library(dplyr)
####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Create and combine objects to use for analysis and plots.
##
## Input:
## Output:
## Date: 2017-09-14
## Author: Jesse R. Walsh
####################################################################################################
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

goAnnotations.sub1.aggr <-
  goAnnotations.sub1 %>%
  select(goTerm) %>%
  group_by(goTerm) %>%
  summarise(n=n())

goAnnotations.sub2.aggr <-
  goAnnotations.sub2 %>%
  select(goTerm) %>%
  group_by(goTerm) %>%
  summarise(n=n())

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

## Convert to v4 ids
expressedGenes <-
  expressedGenes %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, FPKM_mean) %>%
  rename(geneID=v4_id)

## Attach log2(FPKM_mean) values to homeologous pairs
expressedPairs <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  inner_join(expressedGenes, by=c("Maize1"="geneID")) %>%
  rename(FPKM_mean1=FPKM_mean) %>%
  inner_join(expressedGenes, by=c("Maize2"="geneID")) %>%
  rename(FPKM_mean2=FPKM_mean)

#==================================================================================================#
## analyzeGODiffs.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/analyzeGODiffs.R")

#==================================================================================================#
## createTopGO.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/createTopGO.R")

#==================================================================================================#
## txdb -> geneTranscript.map
#--------------------------------------------------------------------------------------------------#
## Only work with chromosomes, ignore unplaced contigs
seqlevels(txdb) <- c("1","2","3","4","5","6","7","8","9","10")

## Get gene/transcript names
geneTranscript.map <- data.frame(transcripts(txdb)$tx_name)
# GRList <- exonsBy(txdb, by = "tx")
# tx_ids <- names(GRList)
# head(select(txdb, keys=tx_ids, columns=c("GENEID","TXNAME"), keytype="TXID"))

## Clean geneTranscript.map
geneTranscript.map <-
  geneTranscript.map %>%
  rename(transcript=transcripts.txdb..tx_name)
geneTranscript.map$transcript <- sub("transcript:", "", geneTranscript.map$transcript)
geneTranscript.map$gene <- sub("(Zm[0-9]{5}d[0-9]{6}).*", "\\1", geneTranscript.map$transcript)
geneTranscript.map <- geneTranscript.map[!startsWith(geneTranscript.map$transcript, "MI"),]
geneTranscript.counts <-
  geneTranscript.map %>%
  select(gene) %>%
  group_by(gene) %>%
  summarise(n=n())

#--------------------------------------------------------------------------------------------------#
detach("package:GenomicFeatures", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
