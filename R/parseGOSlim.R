library(readr)
####################################################################################################
## Project:
## Script purpose:
##
## Input:
## Output:
## Date: 2017-08-18
## Author: Jesse R. Walsh
##
## DataSource: http://www.geneontology.org/ontology/subsets/goslim_plant.obo
####################################################################################################

goSlim <- read_delim("../Data/goslim_plant.tab", "\t", escape_double = FALSE, trim_ws = TRUE)


# goSlim$GOSlimTerm <- paste("|",goSlim$GOSlimTerm,"|",sep = "")

# goAssignmentsIn2not1_withSlim <-
#   converted %>%
#   inner_join(goSlim, by=c("GO Term" = "GOSlimTerm")) %>%
#   filter(grepl("EV-EXP", EVCode)) %>%
#   setdiff(., goAnnotations.sub1 %>%
#               inner_join(goSlim, by=c("GO Term" = "GOSlimTerm")) %>%
#               filter(grepl("EV-EXP", EVCode)))
