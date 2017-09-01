library(topGO)
# go.terms <- as.list(GOTERM)
#---------------------------------------#

## Sub1 genes and their go terms for genes with a homeologous pair
sub1 <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  right_join(., goAnnotations.sub1, by=c("Maize1"="gene2"))
sub1$Maize2[!is.na(sub1$Maize1)] <- as.character(sub1$Maize1[!is.na(sub1$Maize1)])
sub1 <-
  sub1 %>%
  select(Maize2, goTerm, evCode) %>%
  rename(gene2=Maize2)

## Sub2 genes and their go terms for genes with a homeologous pair
sub2 <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "") %>%
  right_join(., goAnnotations.sub2,by=c("Maize2"="gene2"))
sub2$Maize2[!is.na(sub2$Maize1)] <- as.character(sub2$Maize1[!is.na(sub2$Maize1)])
sub2 <-
  sub2 %>%
  select(Maize2, goTerm, evCode) %>%
  rename(gene2=Maize2)

## Filter sub1 and sub2 by desired criteria and look only at what is uniquely assigned to each side
sub1.filtered <-
  sub1 %>%
  # inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>%
  filter(grepl("EV-EXP", evCode)) %>%
  select(gene2, goTerm) %>%
  # filter(goTerm %in% MFterms) %>%
  distinct()

sub2.filtered <-
  sub2 %>%
  # inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>%
  filter(grepl("EV-EXP", evCode)) %>%
  select(gene2, goTerm) %>%
  # filter(goTerm %in% MFterms) %>%
  distinct()

sub1.filtered.unique <- setdiff(sub1.filtered, sub2.filtered)
sub2.filtered.unique <- setdiff(sub2.filtered, sub1.filtered)

## Lets report what we found
for (x in sub1.filtered.unique$gene2) {
  sub1.filtered$goTerm[sub1.filtered$gene2 == x]
  sub2.filtered$goTerm[sub1.filtered$gene2 == x]

  Term(go.terms[["GO:0031976"]])
}

intersect(sub1.filtered$gene2, sub2.filtered$gene2)
