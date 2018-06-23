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
##        B73.v4.centromere.positions.raw
##
## Output:
##        maize.walley.abundance.v4
##        go.maize.clean
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean
##        go.sorghum.clean
##        geneTranscript.counts
##        B73.v4.centromere.positions
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################
library(tidyr)
library(dplyr)
startsWith = getFromNamespace("startsWith", "backports") # if R version < 3.3.0

#==================================================================================================#
## maize.genes.v3_to_v4.map
#--------------------------------------------------------------------------------------------------#
maize.genes.v3_to_v4.map <- maize.genes.v3_to_v4.map[!startsWith(maize.genes.v3_to_v4.map$v4_id, "GRMZM"),]  ## shouldn't need this after next update to MaizeMap

#==================================================================================================#
## maize.protein.abundance.raw
#--------------------------------------------------------------------------------------------------#
maize.walley.abundance.v4 <- ungroup(maize.walley.abundance)

## Convert to v3 ids. At this time, there is no mapping file to do this. Since the vast majority of v2 ids are the same in v3 even if the model itself
## was improved, we can just treat the v2 ids as v3 ids with a few exceptions. New v3 gene models can be ignored, since they weren't in v2 and wouldn't
## map anyway. Removed gene models can also be ignored, since they won't be in v3 and thus will not be converted to v4 ids in the next step.  There
## are 10 renamed gene models and 9 pairs of merged gene models.  These we will manually address below.
## First, the renamed gene models.
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC147602.5_FG004")] <- "GRMZM6G741210"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC190882.3_FG003")] <- "GRMZM6G961377"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC192244.3_FG001")] <- "GRMZM6G869379"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC194389.3_FG001")] <- "GRMZM6G399977"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC204604.3_FG008")] <- "GRMZM6G220418"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC210529.3_FG004")] <- "GRMZM6G945840"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC232289.2_FG005")] <- "GRMZM6G404540"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC233893.1_FG001")] <- "GRMZM6G310687"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC233910.1_FG005")] <- "GRMZM6G729818"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("AC235534.1_FG001")] <- "GRMZM6G798998"

## Second, merged gene models should be summed together by renaming the model that is lost in the merge and summing the new model name by gene+sample
## The merging process picked the name of the merged gene model based on: geneA + geneB = geneA, so that geneB was merged into geneA and kept the geneA's name.
## Therefore, we need to rename geneB <- geneA, then combine the duplicate rows of geneA.
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G103315")] <- "GRMZM2G000964"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G452386")] <- "GRMZM2G045892"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G518717")] <- "GRMZM2G119720"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G020429")] <- "GRMZM2G142383"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G439578")] <- "GRMZM2G319465"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G117517")] <- "GRMZM2G338693"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM5G864178")] <- "GRMZM5G861997"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM2G143862")] <- "GRMZM5G872800"
maize.walley.abundance.v4$geneID[maize.walley.abundance.v4$geneID %in% c("GRMZM5G823855")] <- "GRMZM5G891969"

maize.walley.abundance.v4 <-
  maize.walley.abundance.v4 %>%
  group_by(geneID, sample) %>%
  summarise_all(funs(sum)) %>%
  arrange(geneID) %>%
  ungroup()

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
  arrange(geneID) %>%
  ungroup()

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

#==================================================================================================#
## B73.v4.centromere.positions.raw
#--------------------------------------------------------------------------------------------------#
B73.v4.centromere.positions <- B73.v4.centromere.positions.raw
B73.v4.centromere.positions <-
  B73.v4.centromere.positions %>%
  rename(Chromosome="X__1", Pos_start="X__6", Pos_end="X__7") %>%
  select(Chromosome, Pos_start, Pos_end)

# Remove the extra region on chromosome 9, that is just confusing
B73.v4.centromere.positions <- B73.v4.centromere.positions[-c(10),]

B73.v4.centromere.positions$Chromosome <- as.integer(B73.v4.centromere.positions$Chromosome)

#--------------------------------------------------------------------------------------------------#
## Clean up raw files
rm(MaizeGO.B73.v3_to_v4)
rm(MaizeGO.B73.uniprot_to_v4)
rm(syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw)

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
