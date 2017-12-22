####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Clean up raw datafiles needed for this project, including renaming column and
##        selecting relevant columns, handling missing data and incomplete rows, and any global
##        conversions or modification to data to make it usable.
##
## Input:
##        maize.genes.v3_to_v4_map.raw
##        maize.expression.raw
##        maize.expression.sample.avg.raw
##        maize.protein.abundance.sample.avg.raw
##        maize.kaeppler.expression.raw
##        maize.kaeppler.expression.sample.avg.raw
##        go.maize.raw
##        go.maize.v3.raw
##        syntelogs.sorghum.v1.maize.v1.raw
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        go.sorghum.raw
##        gene.transcript.map
##
## Output:
##        maize.genes.v3_to_v4_map.clean
##        maize.expression.clean
##        maize.expression.sample.avg.clean
##        maize.protein.abundance.sample.avg.clean
##        maize.kaeppler.expression.clean
##        maize.kaeppler.expression.sample.avg.clean
##        go.maize.clean
##        syntelogs.sorghum.v1.maize.v1.clean
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean
##        go.sorghum.clean
##        geneTranscript.counts
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################
library(tidyr)
library(dplyr)
# startsWith = getFromNamespace("startsWith", "backports") # if R version < 3.3.0

#==================================================================================================#
## maize.genes.v3_to_v4_map.raw
#--------------------------------------------------------------------------------------------------#
maize.genes.v3_to_v4_map.clean <- maize.genes.v3_to_v4_map.raw

#==================================================================================================#
## maize.expression.raw
#--------------------------------------------------------------------------------------------------#
maize.expression.clean <- maize.expression.raw

## Convert to v4 ids
maize.expression.clean <-
  maize.expression.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  rename(geneID=v4_id)

## Reorder columns to get the v4 id in column 1 and drop the v3 id
maize.expression.clean <- maize.expression.clean[,c(70,2:69)]

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding FPKM values together
maize.expression.clean <-
  maize.expression.clean %>%
  group_by(geneID) %>%
  summarise_all(funs(sum))

#==================================================================================================#
## maize.expression.sample.avg.raw
#--------------------------------------------------------------------------------------------------#
maize.expression.sample.avg.clean <- ungroup(maize.expression.sample.avg.raw)

## Convert to v4 ids. Keep v4 id, drop v3 id
maize.expression.sample.avg.clean <-
  maize.expression.sample.avg.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, sample, FPKM_avg) %>%
  rename(geneID=v4_id)

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding FPKM values together
maize.expression.sample.avg.clean <-
  maize.expression.sample.avg.clean %>%
  group_by(geneID, sample) %>%
  summarise_all(funs(sum)) %>%
  arrange(geneID)

#==================================================================================================#
## maize.protein.abundance.raw
#--------------------------------------------------------------------------------------------------#
maize.protein.abundance.sample.avg.clean <- ungroup(maize.protein.abundance.sample.avg.raw)

## Convert to v4 ids
maize.protein.abundance.sample.avg.clean <-
  maize.protein.abundance.sample.avg.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, sample, dNSAF_avg) %>%
  rename(geneID=v4_id)

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding FPKM values together
maize.protein.abundance.sample.avg.clean <-
  maize.protein.abundance.sample.avg.clean %>%
  group_by(geneID, sample) %>%
  summarise_all(funs(sum)) %>%
  arrange(geneID)

#==================================================================================================#
## maize.kaeppler.expression.raw
#--------------------------------------------------------------------------------------------------#
maize.kaeppler.expression.clean <- maize.kaeppler.expression.raw

## Convert to v4 ids
maize.kaeppler.expression.clean <-
  maize.kaeppler.expression.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  rename(geneID=v4_id)

## Reorder columns to get the v4 id in column 1 and drop the v3 id
maize.kaeppler.expression.clean <- maize.kaeppler.expression.clean[,c(81,2:80)]

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding FPKM values together
maize.kaeppler.expression.clean <-
  maize.kaeppler.expression.clean %>%
  group_by(geneID) %>%
  summarise_all(funs(sum))

#==================================================================================================#
## maize.kaeppler.expression.sample.avg.raw
#--------------------------------------------------------------------------------------------------#
maize.kaeppler.expression.sample.avg.clean <- ungroup(maize.kaeppler.expression.sample.avg.raw)

## Convert to v4 ids. Keep v4 id, drop v3 id
maize.kaeppler.expression.sample.avg.clean <-
  maize.kaeppler.expression.sample.avg.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, sample, FPKM_avg) %>%
  rename(geneID=v4_id)

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding FPKM values together
maize.kaeppler.expression.sample.avg.clean <-
  maize.kaeppler.expression.sample.avg.clean %>%
  group_by(geneID, sample) %>%
  summarise_all(funs(sum)) %>%
  arrange(geneID)

#==================================================================================================#
## go.maize.raw
#--------------------------------------------------------------------------------------------------#
go.maize.clean <- go.maize.v3.raw

## Convert to v4 ids, remove duplicates that might happen for merged gene models
go.maize.clean <-
  go.maize.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, goTerm, publication, evCode, curator, source, type) %>%
  rename(geneID=v4_id) %>%
  distinct()

## Merge these go annotations the ones already assigned to v4 ids
go.maize.clean <-
  go.maize.clean %>%
  bind_rows(go.maize.raw) %>%
  subset(!is.na(type))

## TODO
## There are about 1000 dulicates of the geneID/goTerm assignment.  Some from type (both exp and comp) and some from source.  Need a way
## to pick which of the duplicates to keep and which to toss.
# go.maize.clean %>%
#   select(geneID, goTerm) %>%
#   distinct()


#==================================================================================================#
## syntelogs.sorghum.v1.maize.v1.raw
#--------------------------------------------------------------------------------------------------#
## ks and kn values are sometimes NA, this will be allowed
syntelogs.sorghum.v1.maize.v1.clean <- syntelogs.sorghum.v1.maize.v1.raw

#==================================================================================================#
## syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
#--------------------------------------------------------------------------------------------------#
## ks and kn values are sometimes NA, this will be allowed
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean <- syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw

## For v4 imports, need to remove the CDS: and _T00# parts of the gene2 column
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2, pattern = "CDS:", replacement = "")
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2, pattern = "_T\\d\\d\\d", replacement = "")

## For sorghum v3.1, need to change identifier to SORBI_ format
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1, pattern = "Sobic.", replacement = "SORBI_")
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1, pattern = ".\\d.v\\d.\\d", replacement = "")

#==================================================================================================#
## go.sorghum.raw
#--------------------------------------------------------------------------------------------------#
go.sorghum.clean <- go.sorghum.raw

go.sorghum.clean <-
  go.sorghum.clean %>%
  rename(sorghumID=`Gene stable ID`, goTerm=`GO term accession`) %>%
  select(sorghumID, goTerm)

#==================================================================================================#
## geneTranscript.counts
#--------------------------------------------------------------------------------------------------#
## Condense gene:transcript map to a gene:countOfTranscripts map
geneTranscript.counts <-
  gene.transcript.map %>%
  select(gene) %>%
  group_by(gene) %>%
  summarise(n=n())

#--------------------------------------------------------------------------------------------------#
## Clean up raw files
rm(maize.genes.v3_to_v4_map.raw)
rm(maize.expression.raw)
rm(maize.expression.sample.avg.raw)
rm(maize.protein.abundance.sample.avg.raw)
rm(maize.kaeppler.expression.raw)
rm(maize.kaeppler.expression.sample.avg.raw)
rm(go.maize.raw)
rm(go.maize.v3.raw)
rm(syntelogs.sorghum.v1.maize.v1.raw)
rm(syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw)
rm(go.sorghum.raw)

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
