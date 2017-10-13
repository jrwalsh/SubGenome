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
# Correlation between datasets
#--------------------------------------------------------------------------------------------------#
## Start with the paired expression data
data <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  distinct() %>%
  inner_join(maize.expression.sample.avg, by=c("Maize1"="geneID")) %>%
  inner_join(maize.expression.sample.avg, by=c("Maize2"="geneID", "Sample"="Sample")) %>%
  subset(!is.na(FPKM_avg.x) & !is.na(FPKM_avg.y))
names(data)[4] <- "FPKM_maize1"
names(data)[5] <- "FPKM_maize2"

## Add paired abundance data
data <-
  data %>%
  inner_join(maize.protein.abundance.sample.avg, by=c("Maize1"="geneID", "Sample"="Sample")) %>%
  inner_join(maize.protein.abundance.sample.avg, by=c("Maize2"="geneID", "Sample"="Sample")) %>%
  subset(!is.na(dNSAF_avg.x) & !is.na(dNSAF_avg.y))
names(data)[6] <- "dNSAF_maize1"
names(data)[7] <- "dNSAF_maize2"

## Fit test -> This data is not normal!  Log-normal seems like a best fit.
descdist(data$FPKM_maize1, discrete = FALSE)
descdist(data$FPKM_maize2, discrete = FALSE)
descdist(data$dNSAF_maize1, discrete = FALSE)
descdist(data$dNSAF_maize2, discrete = FALSE)
fit.lnorm <- fitdist(data$FPKM_maize1, "lnorm", method = "mme")
plot(fit.lnorm)
fit.lnorm$aic
fit <- logspline(stats)

## Correlation test
cor(data$FPKM_maize1, data$FPKM_maize2, method = "pearson")
cor(data$dNSAF_maize1, data$dNSAF_maize2, method = "pearson")
cor(data$FPKM_maize1, data$dNSAF_maize1, method = "pearson")

cor(data$FPKM_maize1, data$FPKM_maize2, method = "spearman")
cor(data$dNSAF_maize1, data$dNSAF_maize2, method = "spearman")
cor(data$FPKM_maize1, data$dNSAF_maize1, method = "spearman")

fit <- lm(data=data, FPKM_maize1 ~ dNSAF_maize1)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)


data$foldChange_expr <- log2(data$FPKM_maize1) - log2(data$FPKM_maize2)
data$foldChange_abun <- log2(data$dNSAF_maize1) - log2(data$dNSAF_maize2)
ggplot(data, aes(foldChange_expr, foldChange_abun)) +
  geom_point(alpha=.1) +
  xlim(-10,10) + ylim(-10,10) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  # geom_abline(mapping = null, data = null, slope = 1, intercept = 0) +
  # scale_x_log10() + scale_y_log10() +
  # facet_wrap(~Sample, nrow=6) +
  labs(
    title = paste0("Scatterplot of Expression in Maize1 vs. Abundance in Maize1\nLog-Log Scale for gene pair(s).\nAlpha=.1, linear abline, loess smoothing"),
    x = "Expression for Maize1 Gene (log scale)",
    y = "Abundance for Maize1 Gene (log scale)"
  ) +
  annotate("text",x=8,y=-8, label=paste("Pearson Cor = ", round(cor(data$foldChange_expr, data$foldChange_abun, method = c("pearson")),2)))

#--------------------------------------------------------------------------------------------------#
detach("package:fitdistrplus", unload=TRUE)
