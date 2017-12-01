library(GenomicFeatures)
library(topGO)
library(tidyr)
library(dplyr)
startsWith = getFromNamespace("startsWith", "backports") # if R version < 3.3.0
####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Clean up raw datafiles needed for this project, including renaming column and
##        selecting relevant columns, handling missing data and incomplete rows, and any global
##        conversions or modification to data to make it usable.
##
## Input:
##        syntelogs.sorghum.v1.maize.v1.raw
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        maize.expression.raw
##        maize.kaeppler.expression.raw
##        go.maize.raw
##        go.sorghum.raw
##        maize.genes.v3_to_v4_map.raw
##        txdb
##        maize.protein.abundance.raw
## Output:
##        syntelogs.sorghum.v1.maize.v1.clean
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean
##        maize.expression.clean
##        maize.kaeppler.expression.clean
##        go.maize.clean
##        go.sorghum.clean
##        maize.genes.v3_to_v4_map.clean
##        geneTranscript.map
##        maize.protein.abundance.clean
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################
#==================================================================================================#
## syntelogs.sorghum.v1.maize.v1.raw
#--------------------------------------------------------------------------------------------------#
syntelogs.sorghum.v1.maize.v1.clean <- syntelogs.sorghum.v1.maize.v1.raw

## ks and kn values are sometimes NA, this will be allowed

#==================================================================================================#
## syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
#--------------------------------------------------------------------------------------------------#
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean <- syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw

## For v4 imports, need to remove the CDS: and _T00# parts of the gene2 column
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2, pattern = "CDS:", replacement = "")
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene2, pattern = "_T\\d\\d\\d", replacement = "")

## For sorghum v3.1, need to change identifier to SORBI_ format
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1, pattern = "Sobic.", replacement = "SORBI_")
syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1 <- gsub(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean$gene1, pattern = ".\\d.v\\d.\\d", replacement = "")

## ks and kn values are sometimes NA, this will be allowed

#==================================================================================================#
## maize.expression.raw
#--------------------------------------------------------------------------------------------------#
maize.expression.clean <- maize.expression.raw

## Remove low FPKM values
maize.expression.clean[maize.expression.clean < 1] <- NA

#==================================================================================================#
## maize.protein.abundance.raw
#--------------------------------------------------------------------------------------------------#
maize.protein.abundance.clean <- maize.protein.abundance.raw

## Rename the gene id column
maize.protein.abundance.clean <-
  maize.protein.abundance.clean %>%
  rename(v3_id=X__1)

## Remove low dNSAF values
maize.protein.abundance.clean[maize.protein.abundance.clean < 1] <- NA

#==================================================================================================#
## maize.kaeppler.expression.raw
#--------------------------------------------------------------------------------------------------#
maize.kaeppler.expression.clean <- maize.kaeppler.expression.raw

## Select relevant columns (i.e. get rid of chromosome position columns)
maize.kaeppler.expression.clean <-
  maize.kaeppler.expression.clean %>%
  select(-chromosome, -position_left, -position_right)

## Remove low FPKM values
maize.kaeppler.expression.clean[maize.kaeppler.expression.clean < 1] <- NA

#==================================================================================================#
## go.maize.raw
#--------------------------------------------------------------------------------------------------#
go.maize.clean <- go.maize.raw

## Rename columns
go.maize.clean <-
  go.maize.clean %>%
  # rename("geneID" = "MaizeCyc2.2 Accession-1") %>%
  rename(geneID = V4_ID) %>%
  rename(goTerm = `GO Term`)

## Parse out citation string and select important columns
go.maize.clean <-
  go.maize.clean %>%
  select(geneID, goTerm, Citation) %>%
  separate(Citation, c("publication", "evCode","timeStamp","curator"), sep=":", extra="drop")

## Publication is uncommon, replace blanks with NA's
go.maize.clean$publication[go.maize.clean$publication == ""] <- NA

## Simplify evidence codes, assume missing is computationally derived.
go.maize.clean$evCode[grepl("EV-EXP", go.maize.clean$evCode)] <- "EXP"
go.maize.clean$evCode[grepl("EV-AS", go.maize.clean$evCode)] <- "EXP" #trust author statements
go.maize.clean$evCode[grepl("EV-IC", go.maize.clean$evCode)] <- "EXP" #trust curator inferences
go.maize.clean$evCode[grepl("EV-COMP", go.maize.clean$evCode)] <- "COMP"
go.maize.clean$evCode[go.maize.clean$evCode == ""] <- "COMP"

## Remove |'s from GO Terms
go.maize.clean$goTerm <- gsub(go.maize.clean$goTerm, pattern = "\\|", replacement = "")

## Add type.  GO Terms are either CC, BP, or MF.  Terms without a type, type is set to NA
go.maize.clean$type <- ""
go.maize.clean$type[go.maize.clean$goTerm %in% ls(GOMFTerm)] <- "MF"
go.maize.clean$type[go.maize.clean$goTerm %in% ls(GOBPTerm)] <- "BP"
go.maize.clean$type[go.maize.clean$goTerm %in% ls(GOCCTerm)] <- "CC"
go.maize.clean$type[go.maize.clean$type == ""] <- NA

go.maize.clean <- distinct(go.maize.clean)

#==================================================================================================#
## go.sorghum.raw
#--------------------------------------------------------------------------------------------------#
go.sorghum.clean <- go.sorghum.raw

go.sorghum.clean <-
  go.sorghum.clean %>%
  rename(sorghumID=`Gene stable ID`, goTerm=`GO term accession`) %>%
  select(sorghumID, goTerm)

#==================================================================================================#
## maize.genes.v3_to_v4_map.raw
#--------------------------------------------------------------------------------------------------#
maize.genes.v3_to_v4_map.clean <- maize.genes.v3_to_v4_map.raw

maize.genes.v3_to_v4_map.clean <-
  maize.genes.v3_to_v4_map.clean %>%
  rename(v3_id = `v3 gene ID`, v4_id = `v4 gene ID (if present)`) %>%
  select(v3_id, v4_id)

maize.genes.v3_to_v4_map.clean$v3_id[maize.genes.v3_to_v4_map.clean$v3_id == "na"] <- NA
maize.genes.v3_to_v4_map.clean$v4_id[maize.genes.v3_to_v4_map.clean$v4_id == "na"] <- NA

maize.genes.v3_to_v4_map.clean$v3_id[startsWith(maize.genes.v3_to_v4_map.clean$v3_id, "AF")] <- NA
maize.genes.v3_to_v4_map.clean$v3_id[startsWith(maize.genes.v3_to_v4_map.clean$v3_id, "AY")] <- NA
maize.genes.v3_to_v4_map.clean$v3_id[startsWith(maize.genes.v3_to_v4_map.clean$v3_id, "EF")] <- NA
maize.genes.v3_to_v4_map.clean$v3_id[startsWith(maize.genes.v3_to_v4_map.clean$v3_id, "zma")] <- NA
maize.genes.v3_to_v4_map.clean$v4_id[startsWith(maize.genes.v3_to_v4_map.clean$v4_id, "zma")] <- NA
maize.genes.v3_to_v4_map.clean$v4_id[maize.genes.v3_to_v4_map.clean$v4_id == "not in v4"] <- NA

maize.genes.v3_to_v4_map.clean <-
  maize.genes.v3_to_v4_map.clean %>%
  filter(!is.na(v3_id) & !is.na(v4_id))

## v3 to v4 is 1 to many mapping
# maize.genes.v3_to_v4_map.clean$v4_id[duplicated(maize.genes.v3_to_v4_map.clean$v4_id)]

#==================================================================================================#
## txdb -> geneTranscript.map
#--------------------------------------------------------------------------------------------------#
## Only work with chromosomes, ignore unplaced contigs
seqlevels(txdb) <- c("1","2","3","4","5","6","7","8","9","10")

## Get gene/transcript names
geneTranscript.map <- data.frame(transcripts(txdb)$tx_name)
# GRList <- exonsBy(txdb, by = "tx")
# tx_ids <- names(GRList)
# head(select(txdb, keys=tx_ids, columns=c("GENEID","TXNAME"), keytype="TXID"))

## Clean geneTranscript.map
geneTranscript.map <-
  geneTranscript.map %>%
  rename(transcript=transcripts.txdb..tx_name)
geneTranscript.map$transcript <- sub("transcript:", "", geneTranscript.map$transcript)
geneTranscript.map$gene <- sub("(Zm[0-9]{5}d[0-9]{6}).*", "\\1", geneTranscript.map$transcript)
geneTranscript.map <- geneTranscript.map[!startsWith(geneTranscript.map$transcript, "MI"),]
geneTranscript.counts <-
  geneTranscript.map %>%
  select(gene) %>%
  group_by(gene) %>%
  summarise(n=n())

rm(txdb)

#--------------------------------------------------------------------------------------------------#
detach("package:GenomicFeatures", unload=TRUE)
detach("package:topGO", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
