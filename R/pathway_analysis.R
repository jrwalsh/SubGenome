#--------------------------------------------------------------------------------------------------#
# Pathways
#--------------------------------------------------------------------------------------------------#
getPathwayReactionOverlap <- function(homeologs.pairs) {
  # ## Load
  corncyc.geneids.clean <- MaizeMap::corncyc.gene.map
  corncyc.reactions.clean <- MaizeMap::corncyc.reaction.gene.map
  corncyc.pathways.clean <- MaizeMap::corncyc.pathway.reaction.map

  # ## Clean
  # Merge the pathway/reaction sets
  pathway.reaction.genes <-
    corncyc.pathways.clean %>%
    inner_join(corncyc.reactions.clean, by=c("ReactionID"="ReactionID"))

  ## Do
  # A dataset with homeologous pairs with their frameID's (where available)
  data <-
    homeologs.pairs %>%
    subset(Maize1 != "" & Maize2 != "") %>% #comment this line to get an all subgenome comparison instead of only those with homeologs
    select(Maize1, Maize2) %>%
    distinct() %>%
    left_join(corncyc.geneids.clean, by=c("Maize1"="v4_id")) %>%
    left_join(corncyc.geneids.clean, by=c("Maize2"="v4_id")) %>%
    rename(FrameID1=FrameID.x, FrameID2=FrameID.y)

  # Reactions associated with Maize1 genes
  df1 <-
    data %>%
    subset(!is.na(FrameID1)) %>%
    left_join(pathway.reaction.genes, by=c("FrameID1"="GeneID")) %>%
    subset(!is.na(ReactionID)) %>%
    select(ReactionID) %>%
    distinct()

  # Reactions associated with Maize2 genes
  df2 <-
    data %>%
    subset(!is.na(FrameID2)) %>%
    left_join(pathway.reaction.genes, by=c("FrameID2"="GeneID")) %>%
    subset(!is.na(ReactionID)) %>%
    select(ReactionID) %>%
    distinct()

  nrow(intersect(df1, df2)) #Shared reactions
  nrow(setdiff(df1, df2)) #Unique to M1
  nrow(setdiff(df2, df1)) #Unique to M2


  # Pathways associated with Maize1 genes
  df3 <-
    data %>%
    subset(!is.na(FrameID1)) %>%
    left_join(pathway.reaction.genes, by=c("FrameID1"="GeneID")) %>%
    subset(!is.na(PathwayID)) %>%
    select(PathwayID) %>%
    distinct()

  # Pathways associated with Maize2 genes
  df4 <-
    data %>%
    subset(!is.na(FrameID2)) %>%
    left_join(pathway.reaction.genes, by=c("FrameID2"="GeneID")) %>%
    subset(!is.na(PathwayID)) %>%
    select(PathwayID) %>%
    distinct()

  descriptions <- c("Shared reactions", "Reactions unique to M1", "Reactions unique to M2", "Shared pathways", "Pathways unique to M1", "Pathways unique to M2")
  values <- c(nrow(intersect(df1, df2)), nrow(setdiff(df1, df2)), nrow(setdiff(df2, df1)), nrow(intersect(df3, df4)), nrow(setdiff(df3, df4)), nrow(setdiff(df4, df3)))
  out <- data.frame(Description=descriptions, Value=values)

  return(out)
}
