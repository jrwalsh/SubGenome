library(readr)
GSE50191_FPKM <- read_delim("~/git/SubGenomes/Data/GSE50191_FPKM.tsv",
                            "\t",
                            trim_ws = TRUE)

GSE50191_FPKM[GSE50191_FPKM < 1] <- NA
expressedGenes <- data.frame(ID=GSE50191_FPKM[,1], Means=rowMeans(GSE50191_FPKM[,-1], na.rm = TRUE))
expressedGenes %>%
  mutate(log2FPKM_mean = log2(Means)) %>%
  select(tracking_id, log2FPKM_mean) %>%
  subset(!is.na(log2FPKM_mean))
