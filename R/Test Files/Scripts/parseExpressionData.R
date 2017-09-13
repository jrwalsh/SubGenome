library(readr)
library(tidyr)
library(dplyr)
# GSE50191_FPKM <- read_delim("~/git/SubGenomes/Data/GSE50191_FPKM.tsv",
#                             "\t",
#                             trim_ws = TRUE)
#
# ## Remove low FPKM values, calculate row means
# GSE50191_FPKM[GSE50191_FPKM < 1] <- NA
expressedGenes <- data.frame(ID=maize.expression.clean[,1], Means=rowMeans(maize.expression.clean[,-1], na.rm = TRUE))

## Remove rows with NA
expressedGenes <-
  expressedGenes %>%
  mutate(FPKM_mean = Means) %>%
  select(tracking_id, FPKM_mean) %>%
  rename(geneID = tracking_id) %>%
  subset(!is.na(FPKM_mean))

## Convert to v4 ids (this is a bad conversion table, only has CornCyc Genes in it.... get a better one!)
geneIDConvertTable <- maize.genes.v3_to_v4_map.clean
  # go.maize.raw %>%
  # select(V4_ID, `MaizeCyc2.2 Accession-1`) %>%
  # rename(v4ID=V4_ID, v3ID=`MaizeCyc2.2 Accession-1`) %>%
  # distinct()

expressedGenes <-
  expressedGenes %>%
  inner_join(geneIDConvertTable, by=c("geneID" = "v3_id")) %>%
  select(v4_id, FPKM_mean) %>%
  rename(geneID=v4_id)
