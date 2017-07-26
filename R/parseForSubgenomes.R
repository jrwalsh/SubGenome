#Input: data input fileName, kscutoff
#Output: file with subgenome

library(readr)
syntelogs.import <- read_delim("~/Dropbox/SynMap Subgenomes/sorghum vs maize.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
syntelogs.complete <- syntelogs.import[complete.cases(syntelogs.import),]


library(plyr)
## Find which sets of chromosomes with syntelogs should be in group 1 or group 2
syntelogs.aggregate <- ddply(syntelogs.complete, ~block+org_chr1+org_chr2, summarise, median=median(ks), mean=mean(ks), count=length(ks))
syntelogs.orthologs <- subset(syntelogs.aggregate, count>=12 & median<=0.4234)
# syntelogs.orthologs[order(syntelogs.orthologs["org_chr1"], syntelogs.orthologs["org_chr2"]),]#view sorted orthologs
orthologs.counts <- ddply(syntelogs.orthologs, ~org_chr1+org_chr2, summarise, total_genes=sum(count))
orthologs.sorted <- orthologs.counts[order(orthologs.counts["org_chr1"], -orthologs.counts["total_genes"]),]

# orthologs.sorted[match(unique(orthologs.sorted$org_chr1), orthologs.sorted$org_chr1),] #first instance of each group (i.e. biggest=sub1)
orthologs.sub1rows <- match(unique(orthologs.sorted$org_chr1), orthologs.sorted$org_chr1) #not exactly right!!!! how to split sub 1 over multiple chromosomes?
orthologs.sub1 <- orthologs.sorted[orthologs.sub1rows,]
orthologs.sub2 <- orthologs.sorted[-orthologs.sub1rows,]
orthologs.sub1$subGenome="sub1"
orthologs.sub2$subGenome="sub2"


## Now recover which genes are in these sets
subset(syntelogs.orthologs, count>=12 & median<=0.4234 & org_chr1=="a31607_1" & org_chr2=="b34889_1")
# orthologs.blocks_groups <- syntelogs.orthologs[syntelogs.orthologs$count>=12 & syntelogs.orthologs$median<=0.4234,]
subgenome.groups <- rbind(merge(syntelogs.orthologs, orthologs.sub1, all.x = F), merge(syntelogs.orthologs, orthologs.sub2, all.x = F))
subgenome <- merge(syntelogs.complete, subgenome.groups, all.x = F)
write.table(subgenome, "~/Desktop/", sep="\t")

library(reshape2)
msyntelogs.complete <- melt(syntelogs.complete, id.vars = c("org_chr1", "org_chr2", "block"))
recast(msyntelogs.complete, block + org_chr1 + org_chr2 ~ variable, mean)
