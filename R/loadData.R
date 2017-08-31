library(readr)
####################################################################################################
## Project:
## Script purpose:
##
## Input:
## Output:
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################

## Import the raw data from the parsed SynMap output
syntelogs.v1.raw <- read_delim("./Data/sorghum_v1_vs_maize_v1.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
syntelogs.v4.raw <- read_delim("./Data/sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

subgenome.chromosomes <- read_delim("./Data/outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## Read in GO Annotation data for maize genes
goAnnotations.raw <- read_delim("./Data/GO from maizecyc.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

sorghum.go <- read_delim("./Data/gramene_sorghumv2_goterms.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

goSlim <- read_delim("./Data/goslim_plant.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

GSE50191_FPKM <- read_delim("./Data/GSE50191_FPKM.tsv", "\t", trim_ws = TRUE)



#--------------------------------------------------------------------------------------------------#
detach(readr)
