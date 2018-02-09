# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")

# library(readr)
# library(GenomicFeatures)
library(fitdistrplus)
library(tidyr)
library(dplyr)
library(ggplot2)
select <- dplyr::select

#--------------------------------------------------------------------------------------------------#
# Local blast searchs
#--------------------------------------------------------------------------------------------------#
cmdBlast <- function(query, subject) {
  paste0("blastn -query <(", query, ") -subject <(", subject, ")")
}

cmdNeedle <- function(query, subject) {
  return(system(paste("perl ./Perl/runNeedle.pl", query, subject)))
}

getSequence <- function(geneID) {

}

parseNeedleOutput <- function() {

}


system("perl ./Perl/runNeedle.pl acacgcacca acgcacgca")

query <- "ATGTCGGGGCGCGGCAAGGGCGGGAAGGGCCTGGGCAAGGGCGGAGCGAAGCGCCACCGGAAGGTGCTCCGCGACAACATCCAGGGGATCACGAAGCCGGCGATCCGGCGGCTGGCGCGCAGGGGCGGCGTGAAGCGCATCTCCGGGCTCATCTACGAGGAGACCCGCGGCGTCCTCAAGATCTTCCTGGAGAACGTGATCCGCGACGCCGTCACCTACACCGAGCACGCGCGCCGCAAGACCGTCACAGCCATGGACGTCGTCTACGCGCTCAAGCGCCAGGGCCGCACCCTCTACGGCTTCGGCGGCTAGGCTATGCCGGTCGCGCTCGCGCCTCGCTGTCGTCGCTGCGGTTCTTGTTTCTGCAATTCGCAAGGTGTGGCTTGACGGGGAAAATGCTGGTTAATGTAGCACTAGATCAAGCATGTTGTGCGTGTCTGAAGAAATGCCAATTGTGACTTGCTGTTTGTTCAAACATGTTCCGTTTCTTCACATATCAAC"
subject <- "ATGGAGATCCTCGAATCGACCCTGTTGGGCGAGTTCATCGGCTTCATCAAGGGGAACTGGTCAGCGCACTCGCGGGTTGACCAGCGCCGGCGCCGCCTGCGCCAGCTGGTCTCGAAGGTGCGCATGGTGGTGGACGCCGCCGAGGGCAATGCGGGGGCGGCCGTGCGGGACGAGTCGTTCTCGGCCTGGCTGCAGGTGCTCCGATCCGAGGCGCTGCGGGGGCAGGAGGTGCTCGACGCCGCGGGCCGCGCCGCGGCCGTCGCGAGCTCCGCCCGCCGCTTCCTCGCGGGCTTCAGGGCGCTCTTCGTCTGCAGCGACGAGGTGGACCGCCTCACGGAGGCCGTCGAGGAGCTCGAGCGCCTGGCGGGGCCCGGCGGCGACCTCGACATGTTCGTCAGGGTCCTCAGCCTGGACGCCGCCCGGACCGCCGCCACGCAGGAGGACATGGACGTCGACGGGCGCCCGGCACCCGGCGCTCGCCACCGACGGGAGGGGAGCGGGTGCGCGGGCTCCGTCGCGGCCCTCCTCCCCTCGCCCGGCGCCAAGAGGAAGCGCGCGTGCGGCGGCGGCTCGACGTCGCGTGGCGAGGACGACGCGGTGCAGCCGCCGAAGCGGCGGCCCCTGGTGTGGATGCGGTCGCACCGGTGGCCCCCGTCCGGCTTCGGCAGGGTCTCTTCCGCACCTCGTGAGCCGCCGCCGGCGGCGTTCCCTCGCTCGCGTCGCGCCCAGACGGTGGCATTGGCGATGTCTAGGATCCGACGCCGCATCGGAAAGCCAACGACGACGAGGCACCGCCGGGAGCCCAGCCTCGGGCAGCATTTCTCGCGGATTACGTTGTAG"
cmdNeedle(query, subject)

#--------------------------------------------------------------------------------------------------#
# GO
#--------------------------------------------------------------------------------------------------#
df1 <- MaizeGO %>%
  subset(evCode=="EXP" & source=="MaizeCyc") %>%
  select(geneID, goTerm) %>%
  distinct()

df2 <- MaizeGO %>%
  subset(evCode=="EXP") %>%
  select(geneID, goTerm) %>%
  distinct()

#--------------------------------------------------------------------------------------------------#
# V4 data
#--------------------------------------------------------------------------------------------------#
maize.expression.sample.avg.clean
data("maize.walley.v4mapped.expression", package = "MaizeOmics")

df <-
  maize.expression.sample.avg.clean %>%
  rename(FPKM_avg_v3 = FPKM_avg) %>%
  inner_join(maize.walley.v4mapped.expression, by=c("geneID"="geneID","sample"="sample")) %>%
  rename(FPKM_avg_v4 = FPKM_avg)

df$FPKM_avg_v3[is.na(df$FPKM_avg_v3)] <- 0
df$FPKM_avg_v4[is.na(df$FPKM_avg_v4)] <- 0

df <-
  df %>%
  subset(!is.na(FPKM_avg_v3) & !is.na(FPKM_avg_v4))

cor(df$FPKM_avg_v3, df$FPKM_avg_v4, method = c("pearson", "kendall", "spearman"))
cor(df$FPKM_avg_v3, df$FPKM_avg_v4, method = "pearson")
cor(df$FPKM_avg_v3, df$FPKM_avg_v4, method = "spearman")
t.test(df$FPKM_avg_v3, df$FPKM_avg_v4, alternative = "less", paired = TRUE)
cor.test(df$FPKM_avg_v3, df$FPKM_avg_v4, method = "spearman")

ggplot(df, aes(x=log2(FPKM_avg_v4), y=log2(FPKM_avg_v3))) +
  geom_point() +
  geom_smooth(method = lm) +
  geom_abline(mapping = null, data = null, slope = 1, intercept = 0, col="red")

#--------------------------------------------------------------------------------------------------#
# Interactive Exploration
#--------------------------------------------------------------------------------------------------#
library(rggobi)
g <- ggobi(data)

# Explore Kaeppler
data <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  distinct() %>%
  inner_join(maize.kaeppler.expression.sample.avg, by=c("Maize1"="geneID")) %>%
  inner_join(maize.kaeppler.expression.sample.avg, by=c("Maize2"="geneID", "Sample"="Sample")) %>%
  subset(!is.na(FPKM_avg.x) & !is.na(FPKM_avg.y))
names(data)[4] <- "FPKM_maize1"
names(data)[5] <- "FPKM_maize2"

## Add paired abundance data
data$foldChange_expr <- log2(data$FPKM_maize1) - log2(data$FPKM_maize2)

#--------------------------------------------------------------------------------------------------#
detach("package:fitdistrplus", unload=TRUE)
