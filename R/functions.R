####################################################################################################
## Project:         SubGenomes project
## Script purpose:  Create functions to be sourced and used for this project
##
## Input:
## Output:
## Date: 2017-09-14
## Author: Jesse R. Walsh
####################################################################################################

saveMyObjects <- function() {
  rm(list = ls())
  source("~/git/SubGenomes/R/loadData.R")
  source("~/git/SubGenomes/R/cleanData.R")
  save.image("~/git/SubGenomes/Data/SavedObjects/loadedData.RData")
  # save(syntelogs.sorghum.v1.maize.v1.clean, file = "~/git/SubGenomes/Data/SavedObjects/syntelogs.sorghum.v1.maize.v1.clean")
  # save(syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean, file = "~/git/SubGenomes/Data/SavedObjects/syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean")
  # save(maize.expression.clean, file = "~/git/SubGenomes/Data/SavedObjects/maize.expression.clean")
  # save(go.maize.clean, file = "~/git/SubGenomes/Data/SavedObjects/go.maize.clean")
  # save(go.sorghum.clean, file = "~/git/SubGenomes/Data/SavedObjects/go.sorghum.clean")
  # save(maize.genes.v3_to_v4_map.clean, file = "~/git/SubGenomes/Data/SavedObjects/maize.genes.v3_to_v4_map.clean")
  # save(go.goSlim.plant, file = "~/git/SubGenomes/Data/SavedObjects/go.goSlim.plant")
  # save(subgenome.assignments, file = "~/git/SubGenomes/Data/SavedObjects/subgenome.assignments")
  # save(subgenome.truth, file = "~/git/SubGenomes/Data/SavedObjects/subgenome.truth")
  # save(txdb, file = "~/git/SubGenomes/Data/SavedObjects/txdb")
}

loadMyObjects <- function() {
  load("~/git/SubGenomes/Data/SavedObjects/loadedData.RData")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/syntelogs.sorghum.v1.maize.v1.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/syntelogs.sorghum.v3.1.maize.v4.and.rejected.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/maize.expression.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/go.maize.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/go.sorghum.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/maize.genes.v3_to_v4_map.clean")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/go.goSlim.plant")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/subgenome.assignments")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/subgenome.truth")
  # load(file = "~/git/SubGenomes/Data/SavedObjects/txdb")
}
