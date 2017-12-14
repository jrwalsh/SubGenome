#--------------------------------------------------------------------------------------------------#
# Pathways
#--------------------------------------------------------------------------------------------------#
## Load
corncyc.pathways.raw <- read_delim("./Data/Pathways/CornCyc_8.0.1_Pathways.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
corncyc.reactions.raw <- read_delim("./Data/Pathways/CornCyc_8.0.1_Reactions.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
corncyc.geneids.raw <- read_delim("./Data/Pathways/CornCyc_8.0.1_GeneIDs.tab", "\t", escape_double = FALSE, trim_ws = TRUE)

## Clean
corncyc.pathways.clean <- corncyc.pathways.raw

# Remove superpathways, remove pathways with no genes, unnest the reaction list
corncyc.pathways.clean <-
  corncyc.pathways.clean %>%
  subset(is.na(SubPathways)) %>%
  subset(!is.na(Reactions)) %>%
  select(PathwayID, Reactions) %>%
  mutate(Reactions = gsub("\"","",Reactions)) %>%
  mutate(ReactionID = strsplit(as.character(Reactions), " // ")) %>%
  unnest(ReactionID) %>%
  select(PathwayID, ReactionID)

corncyc.reactions.clean <- corncyc.reactions.raw

# Remove superpathways, remove pathways with no genes, unnest the reaction list
corncyc.reactions.clean <-
  corncyc.reactions.clean %>%
  subset(!is.na(Genes)) %>%
  mutate(GeneID = strsplit(as.character(Genes), " // ")) %>%
  unnest(GeneID) %>%
  select(ReactionID, GeneID)

corncyc.geneids.clean <- corncyc.geneids.raw

# Collect only GRMZM and Zm names which can be in either the Accession1 or Accession2 column
corncyc.geneids.clean$v3_id <- NA
corncyc.geneids.clean$v4_id <- NA
corncyc.geneids.clean$v3_id <- ifelse(startsWith(corncyc.geneids.clean$Accession1, "GRMZM"),
                                      corncyc.geneids.clean$Accession1,
                                      ifelse(startsWith(corncyc.geneids.clean$Accession2, "GRMZM"),
                                             corncyc.geneids.clean$Accession2,
                                             NA))
corncyc.geneids.clean$v4_id <- ifelse(startsWith(corncyc.geneids.clean$Accession1, "Zm"),
                                      corncyc.geneids.clean$Accession1,
                                      ifelse(startsWith(corncyc.geneids.clean$Accession2, "Zm"),
                                             corncyc.geneids.clean$Accession2,
                                             NA))
# Remove the transcript/protein suffix.
# For GRMZM, it is _, T or P, 2 numbers, and sometimes a "." and another number.
# For Zm, it is _, T or P, and 3 numbers
corncyc.geneids.clean <-
  corncyc.geneids.clean %>%
  mutate(v3_id=gsub("_[PT]\\d\\d.*\\d*", "", v3_id, perl=TRUE)) %>%
  mutate(v4_id=gsub("_[PT]\\d\\d\\d.\\d", "", v4_id, perl=TRUE)) %>%
  select(FrameID, v4_id) %>%
  subset(!is.na(v4_id))


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

nrow(intersect(df1, df2))
nrow(setdiff(df1, df2))
nrow(setdiff(df2, df1))
# write.table(intersect(df1, df2), "~/Desktop/both.cvs")
# write.table(setdiff(df1, df2), "~/Desktop/maize1.cvs")
# write.table(setdiff(df2, df1), "~/Desktop/maize2.cvs")
