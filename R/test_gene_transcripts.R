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

## Plot frequency of isoform counts for each gene
data <- expressedGenes %>% inner_join(subgenome, by=c("geneID"="gene2")) %>% select(geneID, FPKM_mean, subgenome) %>% filter(subgenome=="sub1" | subgenome=="sub2")
ggplot(data, aes(FPKM_mean,color=subgenome)) +
  geom_density() +
  scale_x_log10() +
  labs(y = "Density (FPKM)",
       x = "log10(FPKM)",
       title = "Comparison of FPKM Expression between\nGenes in Subgenome 1 and Subgenome 2"
  )

ggplot(data, aes(subgenome,log10(FPKM_mean))) +
  geom_boxplot() +
  labs(y = "log10(FPKM)",
       x = "Subgenome",
       title = "Comparison of FPKM Expression between\nGenes in Subgenome 1 and Subgenome 2"
  )

##### Any way to compare when sub1>sub2 vs when sub2>sub1
data <- expressedPairs %>% mutate(diff=FPKM_mean1-FPKM_mean2)
ggplot(data, aes(diff)) +
  geom_density() +
  scale_x_log10() +
  labs(y = "Density (FPKM)",
       x = "log10(FPKM)",
       title = "Comparison of FPKM Expression between\nGenes in Subgenome 1 and Subgenome 2"
  )

