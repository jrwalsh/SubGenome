library(readr)
GSE50191_FPKM <- read_delim("~/git/SubGenomes/Data/GSE50191_FPKM.tsv",
                            "\t",
                            trim_ws = TRUE)

## Remove low FPKM values, calculate row means
GSE50191_FPKM[GSE50191_FPKM < 1] <- NA
expressedGenes <- data.frame(ID=GSE50191_FPKM[,1], Means=rowMeans(GSE50191_FPKM[,-1], na.rm = TRUE))

## Remove rows with NA
expressedGenes <-
  expressedGenes %>%
  mutate(FPKM_mean = Means) %>%
  select(tracking_id, FPKM_mean) %>%
  rename(geneID = tracking_id) %>%
  subset(!is.na(FPKM_mean))

## Convert to v4 ids (this is a bad conversion table, only has CornCyc Genes in it.... get a better one!)
geneIDConvertTable <-
  goAnnotations.mutated %>%
  select(geneID, `MaizeCyc2.2 Accession-1`) %>%
  rename(v4ID=geneID, v3ID=`MaizeCyc2.2 Accession-1`) %>%
  distinct()

expressedGenes <-
  expressedGenes %>%
  inner_join(geneIDConvertTable, by=c("geneID" = "v3ID")) %>%
  select(v4ID, FPKM_mean) %>%
  rename(geneID=v4ID)
