library(readr)
gene_transcripts_v4 <- read_delim("~/git/SubGenomes/Data/gene_transcripts_v4.tab",
                                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
# View(gene_transcripts_v4)

gene_transcripts_v4_gff <- read_delim("~/git/SubGenomes/Data/gene_transcripts_v4_gff.tab",
                                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

unique(gene_transcripts_v4_gff)

intersect(gene_transcripts_v4, gene_transcripts_v4_gff)
setdiff(gene_transcripts_v4, gene_transcripts_v4_gff)
setdiff(gene_transcripts_v4_gff, gene_transcripts_v4)
nrow(gene_transcripts_v4) - nrow(gene_transcripts_v4_gff)


library(GenomicFeatures)
gff_file <- system.file("extdata", "GFF3_files", "~/Downloads/Zea_mays.AGPv4.32.gff3.gz", package="GenomicFeatures")
txdb <- makeTxDbFromGFF(gff_file, format="gff3")
