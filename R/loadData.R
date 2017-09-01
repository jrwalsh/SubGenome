library(readr)
library(readxl)
####################################################################################################
## Project:         Subgenomes project
## Script purpose:  Import raw dataset files relevant to this project
##
## Output:
##        syntelogs.sorghum.v1.maize.v1.raw
##        syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw
##        maize.expression.raw
##        go.maize.raw
##        go.sorghum.raw
##        go.goSlim.plant
##        subgenome.assignments
##        subgenome.truth
##        maize.genes.v3_to_v4_map
##
## Date: 2017-08-25
## Author: Jesse R. Walsh
####################################################################################################

syntelogs.sorghum.v1.maize.v1.raw <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\SynMap\\sorghum_v1_vs_maize_v1.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\SynMap\\sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
maize.expression.raw <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\Expression\\GSE50191_FPKM.tsv", "\t", trim_ws = TRUE)
go.maize.raw <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO\\go_from_ maizecyc.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
go.sorghum.raw <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO\\gramene_sorghumv2_goterms.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
go.goSlim.plant <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO\\goslim_plant.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
subgenome.assignments <- read_delim("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
maize.genes.v3_to_v4_map <- read_xlsx("C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\MaizeGDB_v3_v4.genes.xlsx")

# ## Import the raw data from the parsed SynMap output
# syntelogs.sorghum.v1.maize.v1.raw <- read_delim("./Data/SynMap/sorghum_v1_vs_maize_v1.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
# syntelogs.sorghum.v3.1.maize.v4.and.rejected.raw <- read_delim("./Data/SynMap/sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
#
# ## Expression data from the Walley 2016 paper in FPKM over 68 tissues
# maize.expression.raw <- read_delim("./Data/Expression/GSE50191_FPKM.tsv", "\t", trim_ws = TRUE)
#
# ## Read in GO Annotation data for maize genes
# go.maize.raw <- read_delim("./Data/GO/go_from_ maizecyc.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
#
# ## Read in GO Annotation data for sorghum genes
# go.sorghum.raw <- read_delim("./Data/GO/gramene_sorghumv2_goterms.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#
# ## Read in list of plant GO-slim subset
# go.goSlim.plant <- read_delim("./Data/GO/goslim_plant.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
#
# ## This represents subgenome assignments given using a greedy approach (currently performed by hand)
# subgenome.assignments <- read_delim("./Data/outfile_processed_byHand.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## This represents the subgenome assignments given in the paper Schnable 2011
subgenome.truth <- setNames(data.frame(
  c(1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,7,7,8,8,8,9,9,9,10,10,10),
  c(1,5,9,7,2,3,8,5,4,4,2,2,10,1,6,10,4,1,10,3,6,10,8,5,9,6),
  c("sub1","sub2","sub2","sub1","sub2","sub1","sub2","sub1","sub2",
    "sub1","sub2","sub1","sub2","sub1","sub1","sub1","sub2","sub1",
    "sub1","sub2","sub1","sub1","sub2","sub1","sub1","sub2"),
  stringsAsFactors=FALSE), c("chr1","chr2","subgenome"))

# ## Mapping data provided by Maggie, based on synteny from SynMap
# maize.genes.v3_to_v4_map <- read_xlsx("./Data/MaizeGDB_v3_v4.genes.xlsx")

#--------------------------------------------------------------------------------------------------#
detach("package:readr", unload=TRUE)
