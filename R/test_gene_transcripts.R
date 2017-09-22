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
# data <- geneTranscript.counts[geneTranscript.counts$gene %in% subgenome$gene2[subgenome$subgenome=="sub1"],]
data <- geneTranscript.counts %>% inner_join(subgenome, by=c("gene"="gene2")) %>% select(gene, n, subgenome) %>% filter(subgenome=="sub1" | subgenome=="sub2")
ggplot(data, aes(n,color=subgenome)) + geom_density() + scale_x_log10()
ggplot(data, aes(n,color=subgenome)) + geom_histogram(bins = 50)

# data <- geneTranscript.counts[geneTranscript.counts$gene %in% subgenome$gene2[subgenome$subgenome=="sub2"],]
# ggplot(data, aes(n)) + geom_bar()

#==================================================================================================#
##
#--------------------------------------------------------------------------------------------------#
expressedGenes
expressedPairs
getExpressionByExperiment(maize.expression.all, homeologs.pairs, 2)

data <- getExpressionByExperiment(maize.expression.all, homeologs.pairs, 12)
data <- data[!is.na(data$FPKM_maize1) & !is.na(data$FPKM_maize2),]
# data$maize1_Dom <- ifelse(data$FPKM_maize1 > 1*data$FPKM_maize2, TRUE, FALSE)
# nrow(data[data$maize1_Dom,]) - nrow(data[!data$maize1_Dom,])
# data$maize1_Dom <- ifelse(data$FPKM_maize1 > 2*data$FPKM_maize2, TRUE, FALSE)
# nrow(data[data$maize1_Dom,]) - nrow(data[!data$maize1_Dom,])
data$maize1_Dom <- ifelse(data$FPKM_maize1/data$FPKM_maize2 >= 4, TRUE, FALSE)
data$maize2_Dom <- ifelse(data$FPKM_maize2/data$FPKM_maize1 >= 4, TRUE, FALSE)
# data$maize1_Dom <- ifelse((data$FPKM_maize1-data$FPKM_maize2)/data$FPKM_maize2 >= 4, TRUE, FALSE)
# data$maize2_Dom <- ifelse((data$FPKM_maize2-data$FPKM_maize1)/data$FPKM_maize1 >= 4, TRUE, FALSE)
nrow(data[data$maize1_Dom,]) / nrow(data[data$maize1_Dom == FALSE & data$maize2_Dom == FALSE,]) * 100
nrow(data[data$maize2_Dom,]) / nrow(data[data$maize1_Dom == FALSE & data$maize2_Dom == FALSE,]) * 100
nrow(data[data$maize1_Dom == FALSE & data$maize2_Dom == FALSE,])









