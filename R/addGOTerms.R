library(readr)
library(dplyr)
####################################################################################################
## Project: Subgenomes project
## Script purpose:
##
## Input:
## Output:
##
## Date: 2017-07-28
## Author: Jesse R. Walsh
####################################################################################################

goDataFile <- params$goDataFile

goDataFile <- "C:\\Users\\Jesse\\Dropbox (Personal)\\Link to Subgenome Data\\GO from maizecyc.tab"

## Read in GO Annotation data for maize genes
goAnnotations <- read_delim(goDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
goAnnotations.raw <- read_delim(goDataFile, "\t", escape_double = FALSE, trim_ws = TRUE)

## Parse out citation string
goAnnotations.mutated <-
  goAnnotations.raw %>%
  separate(Citation, c("Publication", "EVCode","TimeStamp","Curator"), sep=":", extra="drop")


## GO Analysis -> grouping into useful slices
goAnnotations.sub1 <-
  subgenome %>%
  subset(subgenome == "sub1") %>%
  distinct() %>%
  inner_join(goAnnotations.mutated, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term", EVCode) %>%
  distinct()

goAnnotations.sub2 <-
  subgenome %>%
  subset(subgenome == "sub2") %>%
  distinct() %>%
  inner_join(goAnnotations.mutated, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term", EVCode) %>%
  distinct()

goAnnotations.sub1.aggr <-
  goAnnotations.sub1 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())

goAnnotations.sub2.aggr <-
  goAnnotations.sub2 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())

#
# setdiff(goAnnotations.sub1.aggr %>% select("GO Term"), goAnnotations.sub2.aggr %>% select("GO Term")) %>% nrow()
# setdiff(goAnnotations.sub2.aggr %>% select("GO Term"), goAnnotations.sub1.aggr %>% select("GO Term"))
#
# full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term", is.na=0) %>%
#   rename(sub1_n=n.x, sub2_n=n.y) %>%
#   subset(is.na(sub2_n)) %>%
#   distinct() %>%
#   nrow()
#
# goAnnotations.diff <-
#   full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term", is.na=0) %>%
#   rename(sub1_n=n.x, sub2_n=n.y) %>%
#   subset(is.na(sub1_n) | is.na(sub2_n))
# goAnnotations.diff.counts <-
#   full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term", is.na=0) %>%
#   rename(sub1_n=n.x, sub2_n=n.y) %>%
#   replace_na(list(sub1_n=0, sub2_n=0)) %>%
#   transmute(`GO Term`=`GO Term`,diff=sub1_n-sub2_n) %>%
#   subset(diff != 0)
# goAnnotations.same <-
#   full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term") %>%
#   rename(sub1_n=n.x, sub2_n=n.y) %>%
#   replace_na(list(sub1_n=0, sub2_n=0)) %>%
#   transmute(`GO Term`=`GO Term`,diff=sub1_n-sub2_n) %>%
#   subset(diff == 0)
















# goAnnotations.mutated <- goAnnotations.raw
# r <- regexpr(":EV.*?:", goAnnotations.mutated$Citation)
# goAnnotations.mutated$EVCode <- rep(NA,length(goAnnotations.mutated$Citation))
# goAnnotations.mutated$EVCode[r!=-1] <- regmatches(goAnnotations.mutated$Citation, r)
# goAnnotations.mutated$EVCode <-
#   goAnnotations.mutated$EVCode %>%
#   gsub(pattern=":", replacement="")


# goAnnotations.mutated <-
#   goAnnotations.raw %>%
#   group_by(block, org_chr1, org_chr2) %>%
#   summarise(median_ks=median(ks, na.rm=TRUE), blockGeneCount=n()) %>%
#   left_join(syntelogs.raw, ., by = c("block", "org_chr1", "org_chr2")) %>%
#   mutate(gene_length1=abs(stop1-start1)) %>%
#   mutate(gene_length2=abs(stop2-start2)) %>%
#   mutate(chr1=as.numeric(regmatches(org_chr1, regexpr("\\d*$",org_chr1))))
#   mutate(chr2=as.numeric(regmatches(org_chr2, regexpr("\\d*$",org_chr2))))


## Merge v3 homeolog gene IDs to the GO term data, each line has a unique Gene -> GO term pair, where genes and go terms can be reused
## This list is only genes that have homeologs and have GO annotation
## !! This is an unnessary step, but might be useful later...
# goAnnotations %>%
#   select(V4_ID, `MaizeCyc2.2 Accession-1`, `GO Term`) %>%
#   group_by(V4_ID) %>%
#   merge(homeologs.genes, by.x = "MaizeCyc2.2 Accession-1", by.y = "gene2") %>%
#   select(`MaizeCyc2.2 Accession-1`, `GO Term`) %>%
#   distinct() -> homeologs.go

# ## A simple list of genes in subgenome 1
subgenome.sub1 <-
  subgenome %>%
  filter(subgenome == "sub1") %>%
  select(gene2) %>%
  distinct()
#
# ## A simple list of genes in subgenome 2
subgenome.sub2 <-
  subgenome %>%
  filter(subgenome == "sub2") %>%
  select(gene2) %>%
  distinct()
#
# ## Rows from go annotation file annotating a gene from sub1 or sub2
# goAnnotations.sub1 <-
#   goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub1$gene2,] %>%
#   select(`MaizeCyc2.2 Accession-1`, `GO Term`)
# goAnnotations.sub2 <-
#   goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub2$gene2,] %>%
#   select(`MaizeCyc2.2 Accession-1`, `GO Term`)

subgenome.sub1.go <-
  subgenome %>%
  subset(subgenome == "sub1") %>%
  distinct() %>%
  inner_join(goAnnotations, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term") %>%
  distinct()

subgenome.sub2.go <-
  subgenome %>%
  subset(subgenome == "sub2") %>%
  distinct() %>%
  inner_join(goAnnotations, by = c("gene2" = "MaizeCyc2.2 Accession-1")) %>%
  select(gene2, "GO Term") %>%
  distinct()


# write.table(goAnnotations.sub1, "goAnnotationsSub1.tab", sep="\t")
# write.table(goAnnotations.sub2, "goAnnotationsSub2.tab", sep="\t")
#
# goAnnotations[goAnnotations$`MaizeCyc2.2 Accession-1` %in% subgenome.sub1$gene2,] %>%
#   select(`GO Term`) %>%
#   group_by(`GO Term`) %>%
#   summarise(count(`GO Term`))


## GO Analysis -> grouping into useful slices
goAnnotations.sub1.aggr <-
  goAnnotations.sub1 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())
goAnnotations.sub2.aggr <-
  goAnnotations.sub2 %>%
  select(`GO Term`) %>%
  group_by(`GO Term`) %>%
  summarise(n=n())
goAnnotations.diff <-
  full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term", is.na=0) %>%
  rename(sub1_n=n.x, sub2_n=n.y) %>%
  subset(is.na(sub1_n) | is.na(sub2_n))
goAnnotations.diff.counts <-
  full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term", is.na=0) %>%
  rename(sub1_n=n.x, sub2_n=n.y) %>%
  replace_na(list(sub1_n=0, sub2_n=0)) %>%
  transmute(`GO Term`=`GO Term`,diff=sub1_n-sub2_n) %>%
  subset(diff != 0)
goAnnotations.same <-
  full_join(goAnnotations.sub1.aggr, goAnnotations.sub2.aggr, by="GO Term") %>%
  rename(sub1_n=n.x, sub2_n=n.y) %>%
  replace_na(list(sub1_n=0, sub2_n=0)) %>%
  transmute(`GO Term`=`GO Term`,diff=sub1_n-sub2_n) %>%
  subset(diff == 0)
# write.table(goAnnotations.diff, "goAnnotationDiffs.tab", sep="\t")

print("Number of GO terms in sub1 not in sub2")
print(nrow(goAnnotations.diff %>% subset(is.na(sub2_n))))
print("Number of GO terms in sub2 not in sub1")
print(nrow(goAnnotations.diff %>% subset(is.na(sub1_n))))

# write.table(goAnnotations.diff %>% subset(is.na(sub1_n)), "out.tab", sep="\t")




### GO Analysis by gene
subgenome.homeologs %>%
  group_by(gene1) %>%
  merge(goAnnotations, by.x = "Maize1", by.y = "MaizeCyc2.2 Accession-1") %>%
  select(gene1, Maize1, `GO Term`, Maize2)
