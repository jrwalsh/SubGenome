## Carson asked me to do a quick data exploration to see if I turned up any leads that might indicate support for Carolyn's grant to study roots.
## Basically, I would be looking for under-representation in roots systems vs. the whole, ideally related to the "Gold annotation set" we have.
## -- JRW 11/27/17

library(readr)
library(tidyr)
library(dplyr)
#--------------------------------------------------------------------------------------------------#
# Load Data
go.gold.raw <- read_delim("~/git/MaizeGO/Data/maize_v3.gold.gaf", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)

go.gold.clean <- go.gold.raw

## Rename columns and select important columns
go.gold.clean <-
  go.gold.clean %>%
  rename(geneID = db_object_symbol, goTerm = term_accession, evCode = evidence_code) %>%
  select(geneID, goTerm, evCode)

## Convert to v4
go.gold.clean <-
  go.gold.clean %>%
  inner_join(maize.genes.v3_to_v4_map.clean, by=c("geneID" = "v3_id")) %>%
  select(v4_id, goTerm, evCode) %>%
  rename(geneID=v4_id)

## Remove Duplicates
go.gold.clean <- distinct(go.gold.clean)

#--------------------------------------------------------------------------------------------------#
# Compare to Exression
go.gold.clean

maize.expression.roots <-
  maize.expression.sample.avg %>%
  subset(Sample == "Primary Root 5 Days" | Sample == "Root - Cortex 5 Days" | Sample == "Root - Elongation Zone 5 Days" |
           Sample == "Root - Meristem Zone 5 Days" | Sample == "Secondary Root 7-8 Days") %>%
  subset(!is.na(FPKM_avg))

z.gold.genes <- length(unique(go.gold.clean$geneID))
z.genes.expressed.in.roots <- length(unique(maize.expression.roots$geneID))
z.gold.genes.expressed.in.roots <-
  go.gold.clean %>%
  inner_join(maize.expression.roots, by=("geneID"="geneID")) %>%
  select(geneID) %>%
  unique() %>%
  nrow()

maize.expression.roots$gold <- FALSE
maize.expression.roots$gold[maize.expression.roots$geneID %in% go.gold.clean$geneID] <- TRUE
boxplot(log2(FPKM_avg) ~ gold, data = maize.expression.roots, xlab = "Expressed genes part of Gold Set?", ylab = "log2(FPKM)", main = "Expression in root of genes that are part of the Gold set\nversus those that are not")
summary(maize.expression.roots$FPKM_avg[maize.expression.roots$gold == TRUE])
summary(maize.expression.roots$FPKM_avg[maize.expression.roots$gold == FALSE])

maize.expression.not.roots <-
  maize.expression.sample.avg %>%
  subset(!is.na(FPKM_avg)) %>%
  subset(Sample != "Primary Root 5 Days") %>%
  subset(Sample != "Root - Cortex 5 Days") %>%
  subset(Sample != "Root - Elongation Zone 5 Days") %>%
  subset(Sample != "Root - Meristem Zone 5 Days") %>%
  subset(Sample != "Secondary Root 7-8 Days")

go.gold.clean %>%
  inner_join(maize.expression.not.roots , by=("geneID"="geneID")) %>%
  select(geneID) %>%
  unique() %>%
  nrow()

maize.expression.not.roots$gold <- FALSE
maize.expression.not.roots$gold[maize.expression.not.roots$geneID %in% go.gold.clean$geneID] <- TRUE
boxplot(log2(FPKM_avg) ~ gold, data = maize.expression.not.roots, xlab = "Expressed genes part of Gold Set?", ylab = "log2(FPKM)", main = "Expression in non-root of genes that are part of the Gold set\nversus those that are not")
summary(maize.expression.not.roots$FPKM_avg[maize.expression.not.roots$gold == TRUE])
summary(maize.expression.not.roots$FPKM_avg[maize.expression.not.roots$gold == FALSE])
