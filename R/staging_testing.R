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
