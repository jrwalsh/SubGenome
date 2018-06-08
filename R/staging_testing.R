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
subgenome.positions <-
  subgenome %>%
  left_join(B73.v4.gene.positions, by=c("gene2"="geneID")) %>%
  select(Chromosome, Pos_start, Pos_end, subgenome) %>%
  rename(Subgenome=subgenome) %>%
  subset(Subgenome %in% c("sub1", "sub2"))

gene.pair.types <- getTypeSort(homeologs.pairs, maize.walley.v4mapped.expression, 2, .95)
# dominant_homeolog_exp <- data.frame(geneID=append(gene.pair.types$Maize1[gene.pair.types$type %in% c("1dead")], gene.pair.types$Maize2[gene.pair.types$type %in% c("1dead")]))
dominant_homeolog_exp <- data.frame(geneID=append(gene.pair.types$Maize1[gene.pair.types$type %in% c("1dead") & gene.pair.types$Maize1 %in% dead.genes$geneID],
                                                  gene.pair.types$Maize2[gene.pair.types$type %in% c("1dead") & gene.pair.types$Maize2 %in% dead.genes$geneID]))
dominant_homeolog_exp$Subgenome <- "Dead"
dominant_homeolog_exp <-
  dominant_homeolog_exp %>%
  select(geneID, Subgenome) %>%
  left_join(B73.v4.gene.positions, by=c("geneID"="geneID")) %>%
  select(Chromosome, Pos_start, Pos_end, Subgenome)

df1 <- B73.v4.centromere.positions
df1$x <- -.45
df1$xend <- .45
df2 <- subgenome.positions
df2$x <- -.45
df2$xend <- -.05
df3 <- dominant_homeolog_exp
df3$x <- .05
df3$xend <- .45
segment_data <- bind_rows(df2,df3)
segment_data <- segment_data[!is.na(segment_data$Pos_start),]
segment_data <-
  segment_data %>%
  ungroup()
segment_data$Subgenome[segment_data$Subgenome %in% c("sub1")] <- "Maize1"
segment_data$Subgenome[segment_data$Subgenome %in% c("sub2")] <- "Maize2"
# segment_data$Subgenome[segment_data$Subgenome %in% c("dom")] <- "Dominant Homeolog"
segment_data <- segment_data %>% subset(Subgenome %in% c("Maize1", "Maize2", "Dead"))
centromere_data <- df1

## Make the plot
chromosome_plot_exp <-
  ggplot(data=B73.v4.chr.size) +
  geom_bar(aes(B73.v4.chr.size$chromosome, B73.v4.chr.size$size), stat="identity", fill="grey70") +
  # geom_segment(data=segment_data, aes(x=Chromosome+x, xend=Chromosome+xend, y=Pos_start, yend=Pos_end, colour=Subgenome), size=.001) +
  geom_rect(data=segment_data, inherit.aes = F, aes(xmin=Chromosome+x, xmax=Chromosome+xend, ymin=Pos_start, ymax=Pos_end, color=Subgenome, fill=Subgenome)) +
  geom_rect(data=centromere_data, inherit.aes = F, aes(xmin=Chromosome+x, xmax=Chromosome+xend, ymin=Pos_start, ymax=Pos_end)) +
  scale_fill_manual(values=c("green", "red", "blue", "black")) +
  scale_color_manual(values=c("green", "red", "blue", "black")) +
  scale_x_discrete(limits=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10")) +
  scale_y_continuous(labels=format_si())
chromosome_plot_exp  +
  labs(
    # title="Subgenome Locations by Chromosome",
    x="Chromosome",
    y="Position on Chromosome"
  )

#--------------------------------------------------------------------------------------------------#
detach("package:fitdistrplus", unload=TRUE)
