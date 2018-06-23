####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Create and combine objects to use for analysis and plots.
##
## Input:
## Output:
## Date: 2017-09-14
## Author: Jesse R. Walsh
####################################################################################################
# library(devtools)
# install_github("jrwalsh/MaizeGO@v0.2.0", force = TRUE)
# install_github("jrwalsh/MaizeOmics@v0.2.1", force = TRUE)
# install_github("jrwalsh/MaizeMap@v0.2.1", force = TRUE)
# install.packages("readr")
# install.packages("readxl")
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("VennDiagram")
# install.packages("gridExtra")
## BioConductor Packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
# biocLite("topGO")

source("./R/functions.R")
# saveMyObjects() # only need to run this is there is a change in raw data or the cleaning procedures
loadMyObjects()
source("./R/parseSubgenomes.R")
library(tidyr)
library(dplyr)

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
## parseExpressionData.R
#--------------------------------------------------------------------------------------------------#
## Get means for each gene across all tissues and samples
expressedGenes <- data.frame(ID=maize.walley.v4mapped.expression.replicate[,1], Means=rowMeans(maize.walley.v4mapped.expression.replicate[,-1], na.rm = TRUE))

## Remove rows with NA
expressedGenes <-
  expressedGenes %>%
  rename(FPKM_mean = Means) %>%
  subset(!is.na(FPKM_mean))

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
## parseAbundanceData.R
#--------------------------------------------------------------------------------------------------#
## Get means for each protein across all tissues and samples
maize.walley.abundance.v4.replicate <-
  maize.walley.abundance.v4 %>%
  spread(key=sample, value=dNSAF_avg)
expressedProteins <- data.frame(ID=maize.walley.abundance.v4.replicate[,1], Means=rowMeans(maize.walley.abundance.v4.replicate[,-1], na.rm = TRUE))

## Remove rows with NA
expressedProteins <-
  expressedProteins %>%
  rename(dNSAF_mean = Means) %>%
  subset(!is.na(dNSAF_mean))

## Attach log2(FPKM_mean) values to homeologous pairs
expressedProteinsPairs <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  select(Maize1, Maize2) %>%
  inner_join(expressedProteins, by=c("Maize1"="geneID")) %>%
  rename(dNSAF_mean1=dNSAF_mean) %>%
  inner_join(expressedProteins, by=c("Maize2"="geneID")) %>%
  rename(dNSAF_mean2=dNSAF_mean)

#==================================================================================================#
## Correlation between exp and abundance?
#--------------------------------------------------------------------------------------------------#
## Start with the paired expression data
# data <-
#   homeologs.pairs %>%
#   subset(Maize1 != "" & Maize2 != "") %>%
#   select(Maize1, Maize2) %>%
#   distinct() %>%
#   inner_join(maize.expression.sample.avg, by=c("Maize1"="geneID")) %>%
#   inner_join(maize.expression.sample.avg, by=c("Maize2"="geneID", "Sample"="Sample")) %>%
#   subset(!is.na(FPKM_avg.x) & !is.na(FPKM_avg.y))
# names(data)[4] <- "FPKM_maize1"
# names(data)[5] <- "FPKM_maize2"
#
# ## Add paired abundance data
# data <-
#   data %>%
#   inner_join(maize.protein.abundance.sample.avg, by=c("Maize1"="geneID", "Sample"="Sample")) %>%
#   inner_join(maize.protein.abundance.sample.avg, by=c("Maize2"="geneID", "Sample"="Sample")) %>%
#   subset(!is.na(dNSAF_avg.x) & !is.na(dNSAF_avg.y))
# names(data)[6] <- "dNSAF_maize1"
# names(data)[7] <- "dNSAF_maize2"
#
# ## Fit test -> This data is not normal!  Log-normal seems like a best fit.
# descdist(data$FPKM_maize1, discrete = FALSE)
# descdist(data$FPKM_maize2, discrete = FALSE)
# descdist(data$dNSAF_maize1, discrete = FALSE)
# descdist(data$dNSAF_maize2, discrete = FALSE)
# fit.lnorm <- fitdist(data$FPKM_maize1, "lnorm", method = "mme")
# plot(fit.lnorm)
# fit.lnorm$aic
# fit <- logspline(stats)
#
# ## Correlation test
# cor(data$FPKM_maize1, data$FPKM_maize2, method = "pearson")
# cor(data$dNSAF_maize1, data$dNSAF_maize2, method = "pearson")
# cor(data$FPKM_maize1, data$dNSAF_maize1, method = "pearson")
#
# cor(data$FPKM_maize1, data$FPKM_maize2, method = "spearman")
# cor(data$dNSAF_maize1, data$dNSAF_maize2, method = "spearman")
# cor(data$FPKM_maize1, data$dNSAF_maize1, method = "spearman")
#
# fit <- lm(data=data, FPKM_maize1 ~ dNSAF_maize1)
# summary(fit)
# par(mfrow = c(2, 2))
# plot(fit)

#==================================================================================================#
## analyzeGODiffs.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/analyzeGODiffs.R")

#==================================================================================================#
## createTopGO.R
#--------------------------------------------------------------------------------------------------#
# source("~/git/SubGenomes/R/createTopGO.R")

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
