####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Clean up raw datafiles needed for this project, including renaming column and
##        selecting relevant columns, handling missing data and incomplete rows, and any global
##        conversions or modification to data to make it usable.
##
## Input:
##        maize.walley.abundance
##        MaizeGO.B73.Uniprot, MaizeGO.B73.v3, MaizeGO.B73.v4
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        go.sorghum.raw
##        gene.transcript.map
##
## Output:
##        maize.walley.abundance.v4
##        go.maize.clean
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
## maize.genes.v3_to_v4.map
#--------------------------------------------------------------------------------------------------#
maize.genes.v3_to_v4.map <- maize.genes.v3_to_v4.map[!startsWith(maize.genes.v3_to_v4.map$v4_id, "GRMZM"),]  ## shouldn't need this after next update to MaizeMap

#==================================================================================================#
## maize.protein.abundance.raw
#--------------------------------------------------------------------------------------------------#
maize.walley.abundance.v4 <- ungroup(maize.walley.abundance)

## Convert to v4 ids.  Split genes can have their value applied to both new v4 genes (an implicit side affect of this step),
## while merged genes should have values summed together in next step. Inner join so that non-matches are excluded
maize.walley.abundance.v4 <-
  maize.walley.abundance.v4 %>%
  inner_join(maize.genes.v3_to_v4.map, by=c("geneID" = "v3_id")) %>%
  select(v4_id, sample, dNSAF_avg) %>%
  rename(geneID=v4_id)

## In cases where gene models merge, there will be multiple rows with the same v4 gene model.  Merge duplicate rows by adding dNSAF values together
maize.walley.abundance.v4 <-
  maize.walley.abundance.v4 %>%
  group_by(geneID, sample) %>%
  summarise_all(funs(sum)) %>%
  arrange(geneID)

#==================================================================================================#
## MaizeGO.B73.Uniprot, MaizeGO.B73.v3, MaizeGO.B73.v4
#--------------------------------------------------------------------------------------------------#
go.maize.clean <- MaizeGO.B73.v4

## Convert to v4 ids, remove duplicates that might happen for merged gene models
## Ignore uniprot to v3 mappings, only 2 map to v3 and 0 of those map to v4 anyhow.
MaizeGO.B73.v3_to_v4 <-
  MaizeGO.B73.v3 %>%
  inner_join(maize.genes.v3_to_v4.map, by=c("geneID" = "v3_id")) %>%
  select(v4_id, goTerm, publication, evCode, curator, source, type.x) %>%
  rename(geneID=v4_id, type=type.x) %>%
  distinct()

MaizeGO.B73.uniprot_to_v4 <-
  MaizeGO.B73.Uniprot %>%
  inner_join(maize.genes.uniprot_to_v4.map, by=c("geneID" = "UniProt_Acc")) %>%
  select(v4_id, goTerm, publication, evCode, curator, source, type) %>%
  rename(geneID=v4_id) %>%
  distinct()

## Merge these go annotations and select relevant columns
go.maize.clean <-
  bind_rows(MaizeGO.B73.v3_to_v4, MaizeGO.B73.uniprot_to_v4) %>%
  subset(!is.na(type)) %>%
  select(geneID, goTerm, evCode, type) %>%
  distinct()

## TODO
## For any row that is duplicated in gene + go, mark the row as duplicated.  In a second pass, remove all duplicated that is COMP, then drop the duplicated column

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
rm(MaizeGO.B73.v3_to_v4)
rm(MaizeGO.B73.uniprot_to_v4)
rm(syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw)

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
