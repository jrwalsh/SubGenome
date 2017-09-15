source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

library(readr)
library(GenomicFeatures)
library(tidyr)
library(dplyr)

gene_transcripts_v4 <- read_delim("~/git/SubGenomes/Data/gene_transcripts_v4.tab",
                                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# View(gene_transcripts_v4)

gene_transcripts_v4_gff <- read_delim("~/git/SubGenomes/Data/gene_transcripts_v4_gff.tab",
                                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

unique(gene_transcripts_v4_gff)

intersect(gene_transcripts_v4, gene_transcripts_v4_gff)
setdiff(gene_transcripts_v4, gene_transcripts_v4_gff)
setdiff(gene_transcripts_v4_gff, gene_transcripts_v4)
nrow(gene_transcripts_v4) - nrow(gene_transcripts_v4_gff)




#==================================================================================================#
## geneTranscript.map and geneTranscript.counts
#--------------------------------------------------------------------------------------------------#
## Load the data
gff_file <- system.file("extdata", "GFF3_files", "Zea_mays.AGPv4.32.gff3.gz", package="GenomicFeatures")

## Create the txdb object
txdb <- makeTxDbFromGFF("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\Zea_mays.AGPv4.32.gff3.gz", format="gff3")

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
geneTranscript.map <- geneTranscript.map[!startsWith(geneTranscript.map$transcript, "MI"),]
geneTranscript.map$gene <- sub("(Zm[0-9]{5}d[0-9]{6}).*", "\\1", geneTranscript.map$transcript)
geneTranscript.counts <-
  geneTranscript.map %>%
  select(gene) %>%
  group_by(gene) %>%
  summarise(n=n())

