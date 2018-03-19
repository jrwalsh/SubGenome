retained.duplicates <-
  homeologs.pairs %>%
  subset(Maize1 != "" & Maize2 != "")

test.goAnnotations.sub1 <-
  goAnnotations.sub1[goAnnotations.sub1$gene2 %in% retained.duplicates$Maize1,]

test.goAnnotations.sub2 <-
  goAnnotations.sub2[goAnnotations.sub2$gene2 %in% retained.duplicates$Maize2,]


# All EV-Codes
grid.newpage()
g <- draw.pairwise.venn(
  area1 = test.goAnnotations.sub1 %>% select(goTerm) %>% distinct() %>% nrow(),
  area2 = test.goAnnotations.sub2 %>% select(goTerm) %>% distinct() %>% nrow(),
  cross.area = intersect(test.goAnnotations.sub1 %>% select(goTerm), test.goAnnotations.sub2 %>% select(goTerm)) %>% nrow(),
  category = c("Maize1", "Maize2"),
  lty = rep("blank", 2),
  fill = c("light blue", "green"),
  alpha = rep(0.5, 2),
  cat.pos = c(0, 180),
  scaled = FALSE,
  euler.d = FALSE,
  sep.dist = 0.3,
  rotation.degree = 0,
  ind = FALSE
)
grid.arrange(gTree(children=g), top="GO Term by Subgenome for Retained Duplicates\n(All EV-Codes)")

# # All EV-Codes, GOSlim only
# grid.newpage()
# g <- draw.pairwise.venn(
#   area1 = goAnnotations.sub1 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% select(goTerm) %>% distinct() %>% nrow(),
#   area2 = goAnnotations.sub2 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% select(goTerm) %>% distinct() %>% nrow(),
#   cross.area = intersect(goAnnotations.sub1 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% select(goTerm),
#                          goAnnotations.sub2 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% select(goTerm)) %>%
#     nrow(),
#   category = c("Maize1", "Maize2"),
#   lty = rep("blank", 2),
#   fill = c("light blue", "green"),
#   alpha = rep(0.5, 2),
#   cat.pos = c(0, 180),
#   scaled = FALSE,
#   euler.d = FALSE,
#   sep.dist = 0.3,
#   rotation.degree = 0,
#   ind = FALSE
# )
# grid.arrange(gTree(children=g), top="GO Term by Subgenome\n(All EV-Codes, GOSlim only)")

# # All EV-Codes, MF terms only
# grid.newpage()
# g <- draw.pairwise.venn(
#   area1 = goAnnotations.sub1 %>% filter(type=="MF") %>% select(goTerm) %>% distinct() %>% nrow(),
#   area2 = goAnnotations.sub2 %>% filter(type=="MF") %>% select(goTerm) %>% distinct() %>% nrow(),
#   cross.area = intersect(goAnnotations.sub1 %>% filter(type=="MF") %>% select(goTerm),
#                          goAnnotations.sub2 %>% filter(type=="MF") %>% select(goTerm)) %>%
#     nrow(),
#   category = c("Maize1", "Maize2"),
#   lty = rep("blank", 2),
#   fill = c("light blue", "green"),
#   alpha = rep(0.5, 2),
#   cat.pos = c(0, 180),
#   scaled = FALSE,
#   euler.d = FALSE,
#   sep.dist = 0.3,
#   rotation.degree = 0,
#   ind = FALSE
# )
# grid.arrange(gTree(children=g), top="GO Term by Subgenome\n(All EV-Codes, MF only)")

#EXP only
grid.newpage()
g <- draw.pairwise.venn(
  area1 = test.goAnnotations.sub1 %>% filter(evCode=="EXP") %>% select(goTerm) %>% distinct() %>% nrow(),
  area2 = test.goAnnotations.sub2 %>% filter(evCode=="EXP") %>% select(goTerm) %>% distinct() %>% nrow(),
  cross.area = intersect(test.goAnnotations.sub1 %>% filter(evCode=="EXP") %>% select(goTerm),
                         test.goAnnotations.sub2 %>% filter(evCode=="EXP") %>% select(goTerm)
  ) %>% nrow(),
  category = c("Maize1", "Maize2"),
  lty = rep("blank", 2),
  fill = c("light blue", "green"),
  alpha = rep(0.5, 2),
  cat.pos = c(0, 180),
  scaled = FALSE,
  euler.d = FALSE,
  sep.dist = 0.3,
  rotation.degree = 0,
  ind = FALSE
)
grid.arrange(gTree(children=g), top="GO Term by Subgenome for Retained Duplicates\n(EV-EXP only)")

# #EXP only, GOSlim only
# grid.newpage()
# g <- draw.pairwise.venn(
#   area1 = goAnnotations.sub1 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% filter(evCode=="EXP") %>% select(goTerm) %>% distinct() %>% nrow(),
#   area2 = goAnnotations.sub2 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% filter(evCode=="EXP") %>% select(goTerm) %>% distinct() %>% nrow(),
#   cross.area = intersect(goAnnotations.sub1 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm"))%>% filter(evCode=="EXP") %>% select(goTerm),
#                          goAnnotations.sub2 %>% inner_join(go.goSlim.plant, by=c(goTerm = "GOSlimTerm")) %>% filter(evCode=="EXP") %>% select(goTerm)
#   ) %>% nrow(),
#   category = c("Maize1", "Maize2"),
#   lty = rep("blank", 2),
#   fill = c("light blue", "green"),
#   alpha = rep(0.5, 2),
#   cat.pos = c(0, 180),
#   scaled = FALSE,
#   euler.d = FALSE,
#   sep.dist = 0.3,
#   rotation.degree = 0,
#   ind = FALSE
# )
# grid.arrange(gTree(children=g), top="GO Term by Subgenome\n(EV-EXP, GOSlim only)")

# #EXP only, MF terms only
# grid.newpage()
# g <- draw.pairwise.venn(
#   area1 = goAnnotations.sub1 %>% filter(evCode=="EXP") %>% filter(type=="MF") %>% select(goTerm) %>% distinct() %>% nrow(),
#   area2 = goAnnotations.sub2 %>% filter(evCode=="EXP") %>% filter(type=="MF") %>% select(goTerm) %>% distinct() %>% nrow(),
#   cross.area = intersect(goAnnotations.sub1 %>% filter(evCode=="EXP") %>% filter(type=="MF") %>% select(goTerm),
#                          goAnnotations.sub2 %>% filter(evCode=="EXP") %>% filter(type=="MF") %>% select(goTerm)
#   ) %>% nrow(),
#   category = c("Maize1", "Maize2"),
#   lty = rep("blank", 2),
#   fill = c("light blue", "green"),
#   alpha = rep(0.5, 2),
#   cat.pos = c(0, 180),
#   scaled = FALSE,
#   euler.d = FALSE,
#   sep.dist = 0.3,
#   rotation.degree = 0,
#   ind = FALSE
# )
# grid.arrange(gTree(children=g), top="GO Term by Subgenome\n(EV-EXP, MF only)")
