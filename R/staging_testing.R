# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")

# library(readr)
# library(GenomicFeatures)
library(tidyr)
library(dplyr)


cmdBlast <- function(query, subject){
  paste0("blastn -query <(", query, ") -subject <(", subject, ")")
}

cmdNeedle <- function(query, subject){
  paste0("needle <(echo ", query, ") <(echo ", subject, ") -gapopen 10 -gapextend 0.5 -outfile >(cat)")
}

getSequence <- function(geneID){

}



# ---

graphGenePairDominanceByIsoForms(geneTranscript.counts, homeologs.pairs, 4)

data <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  distinct() %>%
  inner_join(geneTranscript.counts, by=c("Maize1"="gene")) %>%
  inner_join(geneTranscript.counts, by=c("Maize2"="gene"))
names(data)[3] <- "Value_maize1"
names(data)[4] <- "Value_maize2"
data <- subset(data, !is.na(data$Value_maize1) & !is.na(data$Value_maize2))

nSingleExonGenes.total <- nrow(geneTranscript.counts[geneTranscript.counts$n == 1,])
nSingleExonGenes.inSet <- unique(c(data$Maize1[data$Value_maize1 == 1], data$Maize2[data$Value_maize2 == 1]))

data <- subset(data, data$Value_maize1 > 1 & data$Value_maize2 > 1) ## Ignore single-exon genes
data$foldChange <- log2(data$Value_maize1) - log2(data$Value_maize2)
data <- arrange(data, desc(foldChange))
data[1:10,]
