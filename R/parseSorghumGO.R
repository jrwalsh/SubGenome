# sorghum.go <- read_delim("./Data/gramene_sorghumv2_goterms.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#
# sorghum.go <-
#   sorghum.go %>%
#   rename(sorghumID=`Gene stable ID`, goTerm=`GO term accession`) %>%
#   select(sorghumID, goTerm)
#
# sorghum.go$sorghumID <- gsub(sorghum.go$sorghumID, pattern = "SORBI_", replacement = "Sobic.")
# sorghum.go$sorghumID <- gsub(sorghum.go$sorghumID, pattern = "$", replacement = ".1.v3.1")

maizeWithSorghumGO <-
  syntelogs.mutated %>%
  select(gene1, gene2)

# maizeWithSorghumGO$gene1 <- gsub(maizeWithSorghumGO$gene1, pattern = ".\\d.v\\d.\\d", replacement = "")

maizeWithSorghumGO <-
  maizeWithSorghumGO %>%
  left_join(go.sorghum.clean, by=c("gene1"="sorghumID"))

maizeWithSorghumGO <-
  maizeWithSorghumGO %>%
  filter(!is.na(goTerm)) %>%
  select(gene2, goTerm) %>%
  distinct()
