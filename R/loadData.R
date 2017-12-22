####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Import raw dataset files relevant to this project
##
## Output:
##        log10_ks_cutoff
##        subgenome.truth
##
##        maize.genes.v3_to_v4_map.raw
##        maize.expression.raw
##        maize.expression.sample.avg.raw
##        maize.protein.abundance.sample.avg.raw
##        maize.kaeppler.expression.raw
##        maize.kaeppler.expression.sample.avg.raw
##        go.maize.raw
##        go.maize.v3.raw
##        gene.transcript.map
##
##        syntelogs.sorghum.v1.maize.v1.raw
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        go.sorghum.raw
##        go.goSlim.plant
##        subgenome.assignments
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
maize.genes.v3_to_v4_map.raw <- maize.genes.v3_to_v4.map

## Expression data from the Walley 2016 paper in FPKM for 23 tissues
data("maize.walley.expression.replicate", package = "MaizeOmics")
data("maize.walley.expression", package = "MaizeOmics")
maize.expression.raw <- maize.walley.expression.replicate
maize.expression.sample.avg.raw <- maize.walley.expression

## Protein data from the Walley 2016 paper in dNSAF for 33 tissues
data("maize.walley.abundance", package = "MaizeOmics")
maize.protein.abundance.sample.avg.raw <- maize.walley.abundance

## Expression data from the Kaeppler 2015 paper in FPKM for 79 tissues
data("maize.kaeppler.expression.replicate", package = "MaizeOmics")
data("maize.kaeppler.expression", package = "MaizeOmics")
maize.kaeppler.expression.raw <- maize.kaeppler.expression.replicate
maize.kaeppler.expression.sample.avg.raw <- maize.kaeppler.expression

## Read in GO Annotation data for maize genes
data("MaizeGO", package = "MaizeGO")
data("MaizeGO.v3", package = "MaizeGO")
go.maize.raw <- MaizeGO
go.maize.v3.raw <- MaizeGO.v3

## Load maize GFF data
data("gene.transcript.map", package = "MaizeMap")

#--------------------------------------------------------------------------------------------------#
## From Files
## Import the raw data from the parsed SynMap output
syntelogs.sorghum.v1.maize.v1.raw <- read_delim("./data-raw/SynMap/sorghum_v1_vs_maize_v1.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw <- read_delim("./data-raw/SynMap/sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## Read in GO Annotation data for sorghum genes
go.sorghum.raw <- read_delim("./data-raw/GO/gramene_sorghumv2_goterms.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

## Read in list of plant GO-slim subset
go.goSlim.plant <- read_delim("./data-raw/GO/goslim_plant.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## This represents subgenome assignments given using a greedy approach (currently performed by hand)
subgenome.assignments <- read_delim("./data-raw/outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

#--------------------------------------------------------------------------------------------------#
detach("package:readr", unload=TRUE)
detach("package:readxl", unload=TRUE)
detach("package:MaizeGO", unload=TRUE)
detach("package:MaizeMap", unload=TRUE)
detach("package:MaizeOmics", unload=TRUE)
