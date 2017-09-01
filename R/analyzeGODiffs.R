library(topGO)
library(tidyr)
library(dplyr)
# MFterms <- ls(GOMFTerm)
# go.terms <- as.list(GOTERM)
####################################################################################################
## Project:
## Script purpose:  Starting with homeologous pairs, attach GO annotations.  Make comparing easier
##        by converting one set of homeologous pair IDs (sub2) to the equivalent ID from sub1.
##        Based on certain criteria (goSlim, ev-exp, molecular function, etc.) pick genes which are
##        annotated differently between homeologs. Then, back up and compare all their GO annotations
##        (using the unfilted list of GO) and print as a human readable file with GO descriptions.
##
## Input:
## Output:
## Date: 2017-09-01
## Author: Jesse R. Walsh
####################################################################################################

## Sub1 genes and their go terms for genes with a homeologous pair
go.sub1 <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  inner_join(., goAnnotations.sub1, by=c("Maize1"="gene2")) %>%
  select(Maize1, goTerm, evCode) %>%
  rename(gene2=Maize1)

## Sub2 genes and their go terms for genes with a homeologous pair, converts to equivalent gene from sub1 using homeologous pairing
go.sub2 <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  inner_join(., goAnnotations.sub2,by=c("Maize2"="gene2"))
go.sub2$Maize2[!is.na(go.sub2$Maize1)] <- as.character(go.sub2$Maize1[!is.na(go.sub2$Maize1)])
go.sub2 <-
  go.sub2 %>%
  select(Maize2, goTerm, evCode) %>%
  rename(gene2=Maize2)

## Filter sub1 and sub2 by desired criteria and look only at what is uniquely assigned to each side
sub1.filtered <-
  go.sub1 %>%
  # inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>%
  filter(grepl("EV-EXP", evCode)) %>%
  select(gene2, goTerm) %>%
  filter(goTerm %in% MFterms) %>%
  distinct()

sub2.filtered <-
  go.sub2 %>%
  # inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>%
  filter(grepl("EV-EXP", evCode)) %>%
  select(gene2, goTerm) %>%
  filter(goTerm %in% MFterms) %>%
  distinct()

sub1.filtered.unique <- setdiff(sub1.filtered, sub2.filtered) %>% select(gene2)
sub2.filtered.unique <- setdiff(sub2.filtered, sub1.filtered) %>% select(gene2)
# intersect(sub1.filtered$gene2, sub2.filtered$gene2)

# ## This provides a lookup table to get all go terms
# sub1_go_annots <-
#   go.sub1 %>%
#   group_by(gene2) %>%
#   summarise(GO = paste(goTerm, collapse = ", "))
# sub2_go_annots <-
#   go.sub2 %>%
#   group_by(gene2) %>%
#   summarise(GO = paste(goTerm, collapse = ", "))


## Function to return object with gene, go terms, and
getGOAnnots <- function(term, go.sub) {
  list <- go.sub[go.sub$gene2 == term,]
  list$goDesc <- lookupGODesc(go.sub$goTerm[go.sub$gene2 == term])
  list$goType <- lookupGOType(go.sub$goTerm[go.sub$gene2 == term])
  list <-
    list %>%
    unite(goTerm, goTerm, goDesc, evCode, goType, sep = "|") %>%
    select(goTerm)
  # list <-
  #   list %>%
  #   group_by(gene2) %>%
  #   summarise(GO = paste(goTerm, collapse = ", ")) %>%
  #   select(GO)

  # go.sub1[go.sub1$gene2 == term,] %>%
  #   mutate(goDesc = lookupGODesc(goTerm))
  #   group_by(gene2) %>%
  #   summary
  #
  # for (row in go.sub1[go.sub1$gene2 == "Zm00001d034918",]) {
  #   print(row$goTerm)
  # }
  # for (x in sub1.filtered.unique$gene2) {
  #   sub1.filtered$goTerm[sub1.filtered$gene2 == x]
  #   sub2.filtered$goTerm[sub1.filtered$gene2 == x]
  #
  #   print(Term(go.terms[[sub1.filtered$goTerm[sub1.filtered$gene2 == x]]]))
  # }
  if (nrow(list) == 0) {
    return(tibble(goTerm=as.character(NA)))
  }
  return(list)
}

lookupGODesc <- function(goTerms) {
  return(as.vector(Term(goTerms)))
  # return(Term(go.terms[[goTerm]]))
}

lookupGOType <- function(goTerms) {
  return(as.vector(Ontology(goTerms)))
  # return(Term(go.terms[[goTerm]]))
}


## Lets report what we found
compareGOAnnots <- tibble(geneID=character(), GO_sub1=character(), GO_sub2=character())
for (gene in sub1.filtered.unique$gene2) {
  a <- getGOAnnots(gene, go.sub1)
  b <- getGOAnnots(gene, go.sub2)
  c <- intersect(a,b)

  # remove duplicates between both sides
  a <- a %>% filter(!(a$goTerm %in% c$goTerm)) %>% rename(GO_sub1 = goTerm)
  b <- b %>% filter(!(b$goTerm %in% c$goTerm)) %>% rename(GO_sub2 = goTerm)

  # add id column to merge by row
  a <- tibble::rowid_to_column(a, "ID")
  b <- tibble::rowid_to_column(b, "ID")
  geneID <- tibble::rowid_to_column(as.data.frame(gene), "ID") %>% rename(geneID=gene)

  compareGOAnnots <- bind_rows(compareGOAnnots, merge(merge(geneID,a,all=TRUE),b,all=TRUE)[-1])
}

compareGOAnnots$geneID[is.na(compareGOAnnots$geneID)] <- ""
compareGOAnnots$GO_sub1[is.na(compareGOAnnots$GO_sub1)] <- ""
compareGOAnnots$GO_sub2[is.na(compareGOAnnots$GO_sub2)] <- ""

write.table(compareGOAnnots, "compareGOAnnots.tab", sep="\t", row.names=FALSE)
