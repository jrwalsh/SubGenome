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
data <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  distinct() %>%
  inner_join(maize.walley.v4mapped.expression.data, by=c("Maize1"="geneID")) %>%
  inner_join(maize.walley.v4mapped.expression.data, by=c("Maize2"="geneID", "sample"="sample")) %>%
  subset(!is.na(FPKM_avg.x) & !is.na(FPKM_avg.y))
names(data)[4] <- "Value_maize1"
names(data)[5] <- "Value_maize2"

data$m1 <- 0
data$m2 <- 0
data$m1[data$Value_maize1 > data$Value_maize2] <- 1
data$m2[data$Value_maize1 < data$Value_maize2] <- 1
data.temp <- data %>% select(sample, m1, m2) %>% group_by(sample) %>% summarise(vm1=sum(m1), vm2=sum(m2)) %>% ungroup()
data.temp$per <- data.temp$vm1/(data.temp$vm1 + data.temp$vm2)
#--------------------------------------------------------------------------------------------------#
plot <- graphGenePairDominanceByExperiment(maize.walley.v4mapped.expression %>% select(geneID, sample, FPKM_avg), homeologs.pairs, 4, FALSE)
plot(plot)

data <-
  homeologs.pairs %>%
  # subset(Maize1 != "" & Maize2 != "") %>% #only matching
  subset(!(Maize1 %in% c("") & Maize2 %in% c(""))) %>% #only nonmatching
  distinct() %>%
  select(Maize1, Maize2) %>%
  left_join(maize.walley.v4mapped.expression, by=c("Maize1"="geneID")) %>%
  rename(FPKM_avg1=FPKM_avg) %>%
  left_join(maize.walley.v4mapped.expression, by=c("Maize2"="geneID", "sample"="sample")) %>%
  rename(FPKM_avg2=FPKM_avg) %>%
  select(Maize1, Maize2, sample, FPKM_avg1, FPKM_avg2) %>%
  # subset(sample %in% c("B73 Mature Pollen")) %>%
  group_by(sample) %>%
  summarize(Mean1=mean(FPKM_avg1, na.rm=TRUE), Mean2=mean(FPKM_avg2, na.rm=TRUE))
data <-
  data %>%
  gather(subgenome, FPKM, -1)

ggplot(data, aes(x=sample, y=FPKM, fill=subgenome)) + geom_col(position = position_dodge()) + theme(axis.text.x=element_text(angle=90,hjust=1))




data <-
  homeologs.pairs %>%
  subset(Maize2 %in% c("")) %>% #only nonmatching
  select(Maize1) %>%
  left_join(maize.walley.v4mapped.expression, by=c("Maize1"="geneID")) %>%
  select(Maize1, sample, FPKM_avg)
temp <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>% #only matching
  select(Maize1) %>%
  left_join(maize.walley.v4mapped.expression, by=c("Maize1"="geneID")) %>%
  select(Maize1, sample, FPKM_avg)

data$status <- "nonMatch"
temp$status <- "match"
data <-
  data %>%
  bind_rows(temp)

ggplot(data, aes(status,log10(FPKM_avg))) +
  geom_boxplot() +
  labs(y = "log10(FPKM)",
       x = "Subgenome",
       title = "Comparison of FPKM Expression between\nGenes in Subgenome 1 and Subgenome 2"
  )

#--------------------------------------------------------------------------------------------------#
detach("package:fitdistrplus", unload=TRUE)
