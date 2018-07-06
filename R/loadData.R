####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Import raw dataset files relevant to this project
##
## Output:
##        log10_ks_cutoff
##        subgenome.truth
##
##        maize.genes.v3_to_v4.map
##        gene.transcript.map
##        maize.genes.uniprot_to_v4.map
##        maize.walley.v4mapped.expression.replicate
##        maize.walley.v4mapped.expression
##        maize.walley.abundance
##        MaizeGO.B73.Uniprot
##        MaizeGO.B73.v3
##        MaizeGO.B73.v4
##
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        subgenome.assignments
##        chromosome.v4.size
##        gene.positions.v4.chromosome
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################
library(readr)
library(readxl)
library(MaizeGO)
library(MaizeOmics)
library(MaizeMap)

#--------------------------------------------------------------------------------------------------#
## Local parameters
## log10_ks_cutoff, hand picked based on SynMap reported ks/kn data
log10_ks_cutoff <- 0

## This represents the subgenome assignments given in the paper Schnable 2011
subgenome.truth <- setNames(data.frame(
  c(1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,7,7,8,8,8,9,9,9,10,10,10),
  c(1,5,9,7,2,3,8,5,4,4,2,2,10,1,6,10,4,1,10,3,6,10,8,5,9,6),
  c("sub1","sub2","sub2","sub1","sub2","sub1","sub2","sub1","sub2",
    "sub1","sub2","sub1","sub2","sub1","sub1","sub1","sub2","sub1",
    "sub1","sub2","sub1","sub1","sub2","sub1","sub1","sub2"),
  stringsAsFactors=FALSE), c("chr1","chr2","subgenome"))

#--------------------------------------------------------------------------------------------------#
## From Packages
## Mapping data provided by Maggie, based on synteny from SynMap
data("maize.genes.v3_to_v4.map", package = "MaizeMap")

## Load maize GFF data
data("gene.transcript.map", package = "MaizeMap")

## Uniprot Mappings
data("maize.genes.uniprot_to_v4.map", package = "MaizeMap")

## Expression data from the Walley 2016 paper in FPKM for 23 tissues
data("maize.walley.v4mapped.expression.replicate", package = "MaizeOmics")
data("maize.walley.v4mapped.expression", package = "MaizeOmics")

## Protein data from the Walley 2016 paper in dNSAF for 33 tissues
data("maize.walley.abundance", package = "MaizeOmics")

## Read in GO Annotation data for maize genes
data("MaizeGO.B73.Uniprot", package = "MaizeGO")
data("MaizeGO.B73.v3", package = "MaizeGO")
data("MaizeGO.B73.v4", package = "MaizeGO")

#--------------------------------------------------------------------------------------------------#
## From Files
## Import the raw data from the parsed SynMap output
syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw <- read_delim("./data-raw/SynMap/sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## This represents subgenome assignments given using a greedy approach (currently performed by hand)
subgenome.assignments <- read_delim("./data-raw/outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## Data for building chromosome plots
chromosome.v4.size <- read.csv("./data-raw/ChrImages/Zea_mays.AGPv4.36.sizes.csv")
gene.positions.v4.chromosome <- v4_subg <- read.csv("./data-raw/ChrImages/gene_positions.tab", sep = "\t") # this should really be in MaizeMap

B73.v4.chr.size <- read.csv("./data-raw/ChrImages/Zea_mays.AGPv4.36.sizes.csv")
B73.v4.gene.positions <- read_delim("./data-raw/ChrImages/gene_positions.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
B73.v4.centromere.positions.raw <- read_excel("data-raw/ChrImages/centromere_positions.xlsx", col_names = FALSE)

#--------------------------------------------------------------------------------------------------#
detach("package:readr", unload=TRUE)
detach("package:readxl", unload=TRUE)
detach("package:MaizeGO", unload=TRUE)
detach("package:MaizeMap", unload=TRUE)
detach("package:MaizeOmics", unload=TRUE)
