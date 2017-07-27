## Maize v4+rejected set vs. Sorghum v3.1 (updated to use up-to-date assemblies and gene models)

```{r parseSynMap2, include=FALSE}
cmd <- paste("perl", "parseKsKnFile.pl", "6807_8082.CDS-CDS.lastz.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks")
system(cmd)
```

```{r importSynMapData2, echo=FALSE}
## Import data and remove incomplete rows
syntelogs.import <- read_delim("../Data/sorghum_v3.1_vs_maize_v4+rejected.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
syntelogs.complete <- syntelogs.import[complete.cases(syntelogs.import),]

sprintf("We have imported %d lines from the parsed SynMap file and removed %d incomplete rows resulting in %d syntenic genes.", nrow(syntelogs.import), nrow(syntelogs.import)- nrow(syntelogs.complete), nrow(syntelogs.complete))
```

We need to select a cutoff ks value to differentiate between the most recent alpha duplication event and the previous beta duplication.

```{r pickKSCutoff2, echo=FALSE}
## Pick a ks value which separates alpha duplication from beta duplication (between the first two peaks)
log10_ks_cutoff <- -0.575
qplot(log10(syntelogs.complete$ks), geom="histogram", binwidth=.01)  + geom_vline(mapping = NULL, data = NULL, xintercept = log10_ks_cutoff, na.rm = FALSE, show.legend = NA)
ks_cutoff <- 10^log10_ks_cutoff
#ks_cutoff <- 0.4234 #include most of the first 2 peaks
```

Find subgenomes

```{r findSubgenomes2, echo=FALSE}
## Find which sets of chromosomes with syntelogs should be in group 1 or group 2, where group 1 has larger syntenic blocks
# First group rows with same syntenic block and chromosome, then summarize each block's ks values with median/mean/count
syntelogs.aggregate <- ddply(syntelogs.complete, ~block+org_chr1+org_chr2, summarise, median=median(ks), mean=mean(ks), count=length(ks))

# Following the schnable article, syntenic blocks must have 12 genes "The median synonymous substitution rate of all gene pairs in a syntenic block between maize and sorghum can be used to classify syntenic blocks of 12 or more genes unambiguously as orthologous or homoeologous, however" and have a median ks value that discriminates for the alpha duplication event.
syntelogs.orthologs <- subset(syntelogs.aggregate, count>=12 & median<=ks_cutoff)

# Determine which chromosome has largest syntenic block, the syntelogs on this chromosome are subgenome 1
orthologs.counts <- ddply(syntelogs.orthologs, ~org_chr1+org_chr2, summarise, total_genes=sum(count))
orthologs.sorted <- orthologs.counts[order(orthologs.counts["org_chr1"], -orthologs.counts["total_genes"]),]

# orthologs.sorted[match(unique(orthologs.sorted$org_chr1), orthologs.sorted$org_chr1),] #first instance of each group (i.e. biggest=sub1)
orthologs.sub1rows <- match(unique(orthologs.sorted$org_chr1), orthologs.sorted$org_chr1) #not exactly right!!!! how to split sub 1 over multiple chromosomes?
orthologs.sub1 <- orthologs.sorted[orthologs.sub1rows,]
orthologs.sub2 <- orthologs.sorted[-orthologs.sub1rows,]
orthologs.sub1$subGenome="sub1"
orthologs.sub2$subGenome="sub2"
```

Write to file

```{r writeSubgenomesToFile2, echo=FALSE}
## Now recover which genes are in these sets
# subset(syntelogs.orthologs, count>=12 & median<=ks_cutoff & org_chr1=="a31607_1" & org_chr2=="b34889_1")
# orthologs.blocks_groups <- syntelogs.orthologs[syntelogs.orthologs$count>=12 & syntelogs.orthologs$median<=ks_cutoff,]
subgenome.groups <- rbind(merge(syntelogs.orthologs, orthologs.sub1, all.x = F), merge(syntelogs.orthologs, orthologs.sub2, all.x = F))
subgenome <- merge(syntelogs.complete, subgenome.groups, all.x = F)
write.table(subgenome, "../Output/recoveredSubGenomes_sorghum_v3.1_vs_maize_v4+rejected.tab", sep="\t")
```

















library(reshape2)
msyntelogs.complete <- melt(syntelogs.complete, id.vars = c("org_chr1", "org_chr2", "block"))
recast(msyntelogs.complete, block + org_chr1 + org_chr2 ~ variable, mean)
