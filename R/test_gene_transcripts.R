source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")

library(readr)
library(GenomicFeatures)
library(tidyr)
library(dplyr)



# ## Fix experiment names
# mm <- match(names(maize.expression.avg), experiment.map$tracking_id)
# names(maize.expression.avg)[!is.na(mm)] <- as.character(experiment.map$Replicate[na.omit(mm)])
#
#
# maize.expression.sample.avg <- maize.expression.all
# maize.expression.sample.avg <-
#   maize.expression.sample.avg %>%
#   gather("tracking_id", "FPKM",-1) %>%
#   left_join(experiment.map, by=c("tracking_id"="tracking_id")) %>%
#   select(geneID, Sample, FPKM) %>%
#   group_by(geneID, Sample) %>%
#   summarise(FPKM_avg=mean(FPKM, na.rm=TRUE)) %>%
#   arrange(geneID)
# maize.expression.sample.avg$FPKM_avg[is.nan(maize.expression.sample.avg$FPKM_avg)] <- NA
#
# plot <- graphGenePairExpressionsByExperiement(maize.expression.sample.avg, homeologs.pairs, 4, FALSE)
# plot(plot)


library(readxl)
maize.protein.abundance.raw <- read_xlsx("~/git/SubGenomes/Data/ProteinAbundance/aag1125_SupportingFile_Table_S2-1.xlsx", sheet = 3, col_names = TRUE)
maize.protein.abundance.clean <- maize.protein.abundance.raw
maize.protein.abundance.clean <-
  maize.protein.abundance.clean %>%
  select(1:9, 25:36, 41:50, 55:60, 68:85, 90:116, 121:136, 141:144) %>%
  rename(v3_id=X__1)

maize.protein.abundance.clean <-
  maize.protein.abundance.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("v3_id" = "v3_id")) %>%
  rename(geneID=v4_id)

maize.protein.abundance.clean <- maize.protein.abundance.clean[,c(103,2:102)]

## Rejected: no match in expression dataset
names(maize.protein.abundance.clean %>% select(-c(1:9, 25:36, 41:50, 55:60, 68:85, 90:116, 121:136, 141:144)))


