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
# CSHL Syntelogs vs. My SynMap results
#--------------------------------------------------------------------------------------------------#
library(readr)
CSHLPairs <- read_csv("~/Dropbox/CSHLPairs.csv")
CSHLPairs <- rename(CSHLPairs, gene1=ort2, gene2=ort1)
MyPairs <- syntelogs.mutated %>%
  select(gene1, gene2)

intersect(CSHLPairs[,1], MyPairs[,1])
intersect(unique(CSHLPairs[,1]), unique(MyPairs[,1]))

View(CSHLPairs)
View(MyPairs)

#--------------------------------------------------------------------------------------------------#
# A function that can catch all alternating gene dominances
#--------------------------------------------------------------------------------------------------#
# Type 1: subA gene always beats subB gene
  # Implies dieing subB gene
# Type 2: subB gene always beats subA gene
  # Implies dieing subA gene
# Type 3: Dominance switches, but both genes express
# *Type 4: Dominance switches, but one gene is turned off in each case

maize.expression.sample.avg
expressedPairs
retained.duplicates <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "")

getRetainedDuplicatedStats <- function(retained.duplicates, maize.expression.sample.avg) {
  retained.duplicates.stats <- cbind(retained.duplicates$Maize1, retained.duplicates$Maize2)
  retained.duplicates.stats <- unique(retained.duplicates.stats)
  retained.duplicates.stats

  test.expression <- maize.expression.sample.avg
  test.expression$FPKM_avg[is.na(test.expression$FPKM_avg)] <- 0

  test.expression %>%
    group_by(geneID) %>%
    summarise(min = min(FPKM_avg), max = max(FPKM_avg), avg = mean(FPKM_avg), var = var(FPKM_avg))
}

getType3 <- function(retained.duplicates, maize.expression.sample.avg) {
  test.index <- 0
  test.cutoff <- 2
  test.expression <- maize.expression.sample.avg
  test.expression$FPKM_avg[is.na(test.expression$FPKM_avg)] <- 0

  test.subset <-
    retained.duplicates[test.index,] %>%
    inner_join(test.expression, by=c("Maize1"="geneID")) %>%
    rename(FPKM_avg1=FPKM_avg) %>%
    inner_join(test.expression, by=c("Maize2"="geneID", "sample"="sample")) %>%
    rename(FPKM_avg2=FPKM_avg)

  test.subset.cor <-
    test.subset %>%
    select(sample, FPKM_avg1, FPKM_avg2)

  test.gene1 <- retained.duplicates[test.index,2]
  test.gene2 <- retained.duplicates[test.index,3]
  test.cor <- cor(x=test.subset.cor$FPKM_avg1, y=test.subset.cor$FPKM_avg2)

  test.subset.melt <-
    test.subset %>%
    select(sample, FPKM_avg1, FPKM_avg2) %>%
    gather(key = subgenome, "FPKM_avg", 2:3)

  ggplot(data=test.subset.melt, aes(x=sample, y=FPKM_avg, group=subgenome, color=subgenome)) + geom_line() + geom_point() +
    # annotate("text", x = Inf, y = 10, label = test.cor) +
    annotate("rect", xmin = .5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Green") +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
    annotate("rect", xmin = 3.5, xmax = 9.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
    annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
    annotate("rect", xmin = 10.5, xmax = 11.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
    annotate("rect", xmin = 11.5, xmax = 15.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "DarkGreen") +
    annotate("rect", xmin = 15.5, xmax = 16.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Yellow") +
    annotate("rect", xmin = 16.5, xmax = 21.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Brown") +
    annotate("rect", xmin = 21.5, xmax = 22.5, ymin = -Inf, ymax = Inf, alpha = .2, fill = "Orange") +
    annotate("rect", xmin = 22.5, xmax = 23.5, ymin = -Inf, ymax = Inf, alpha = .1, fill = "Green") +
    labs(title="FPKM across tissues for putative retained duplicate in B73", y="FPKM Average(Bioreps)", x="Sample", subtitle = paste0(test.gene1, " and ", test.gene2,": Pearson correlation = ", test.cor)) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

  test.index
  test.index <- test.index + 1

  return()
}


#--------------------------------------------------------------------------------------------------#
detach("package:fitdistrplus", unload=TRUE)
