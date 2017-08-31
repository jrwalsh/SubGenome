library(tidyr)
library(dplyr)
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
##        go.maize.raw
##        go.sorghum.raw
## Output:
##        syntelogs.sorghum.v1.maize.v1.clean
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean
##        maize.expression.clean
##        go.maize.clean
##        go.sorghum.clean
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################
#==================================================================================================#
## syntelogs.sorghum.v1.maize.v1.raw
#--------------------------------------------------------------------------------------------------#
syntelogs.sorghum.v1.maize.v1.clean <- syntelogs.sorghum.v1.maize.v1.raw

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

#==================================================================================================#
## maize.expression.raw
#--------------------------------------------------------------------------------------------------#
maize.expression.clean <- maize.expression.raw

## Remove low FPKM values, calculate row means
maize.expression.clean[maize.expression.clean < 1] <- NA

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

## Remove |'s from GO Terms
go.maize.clean$goTerm <- gsub(go.maize.clean$goTerm, pattern = "\\|", replacement = "")

#==================================================================================================#
## go.sorghum.raw
#--------------------------------------------------------------------------------------------------#
go.sorghum.clean <- go.sorghum.raw

go.sorghum.clean <-
  go.sorghum.clean %>%
  rename(sorghumID=`Gene stable ID`, goTerm=`GO term accession`) %>%
  select(sorghumID, goTerm)

#--------------------------------------------------------------------------------------------------#
detach("package:tidyr", unload=TRUE)
detach("package:dplyr", unload=TRUE)
