####################################################################################################
## Project:
## Script purpose:
##
## Input:
## Output:
## Date: 2018-05-17
## Author: Jesse R. Walsh
####################################################################################################

library(ggplot2)
library(tidyr)
library(dplyr)

library(readr)
library(readxl)

v4 <- read.csv("./data-raw/ChrImages/Zea_mays.AGPv4.36.sizes.csv")
v4_subg <- read_csv("data-raw/ChrImages/v4_subg.csv")
positions <- v4_subg <- read.csv("./data-raw/ChrImages/gene_positions.tab", sep = "\t")
OrphanGenes <- read_excel("data-raw/ChrImages/OrphanGenes.xlsx")

v4_subg_jw <-
  subgenome %>%
  left_join(positions, by=c("gene2"="geneID")) %>%
  select(Chromosome, Pos_start, Pos_end, subgenome) %>%
  rename(Subgenome=subgenome) %>%
  subset(Subgenome %in% c("sub1", "sub2"))

v4_dom_jw <- expressedPairs
v4_dom_jw$Subgenome <- "no dominance"
v4_dom_jw$Subgenome[v4_dom_jw$FPKM_mean1 / v4_dom_jw$FPKM_mean2 > 2] <- "dom"
v4_dom_jw$Subgenome[v4_dom_jw$FPKM_mean2 / v4_dom_jw$FPKM_mean1 > 2] <- "dom"
v4_dom_jw$Maize1[v4_dom_jw$FPKM_mean2 / v4_dom_jw$FPKM_mean1 > 2] <- v4_dom_jw$Maize2[v4_dom_jw$FPKM_mean2 / v4_dom_jw$FPKM_mean1 > 2]
df0 <- v4_dom_jw[v4_dom_jw$Subgenome %in% c("no dominance"),]
df0$Maize1 <- df0$Maize2
v4_dom_jw <- bind_rows(v4_dom_jw, df0)
v4_dom_jw <-
  v4_dom_jw %>%
  select(Maize1, Subgenome) %>%
  left_join(positions, by=c("Maize1"="geneID")) %>%
  select(Chromosome, Pos_start, Pos_end, Subgenome)

## 6,773 - 3,875 (no map) -~-  1,018 (no positions); split/merge is fine = 1,896 orphan v4_ids with positions based on Zea_mays.AGPv4.36.chr.gff3
v4_orph_jw <-
  OrphanGenes %>%
  select(gene) %>%
  left_join(maize.genes.v3_to_v4_map.clean, c("gene" = "v3_id"))
v4_orph_na <- v4_orph_jw[is.na(v4_orph_jw$v4_id),]
v4_orph_jw <- v4_orph_jw[!is.na(v4_orph_jw$v4_id),]
v4_orph_jw <-
  v4_orph_jw %>%
  select(v4_id) %>%
  left_join(positions, by=c("v4_id"="geneID"))
v4_orph_jw <- v4_orph_jw[!is.na(v4_orph_jw$Pos_end),]
v4_orph_jw$Subgenome <- "orphan"
v4_orph_jw <-
  v4_orph_jw %>%
  select(Chromosome, Pos_start, Pos_end, Subgenome)

centromere.v4.positions <- read_excel("data-raw/ChrImages/centromere_positions.xlsx", col_names = FALSE)
centromere.v4.positions <-
  centromere.v4.positions %>%
  rename(Chromosome="X__1", Pos_start="X__6", Pos_end="X__7") %>%
  select(Chromosome, Pos_start, Pos_end)
centromere.v4.positions


#--------------------------------------------------------------------------------------------------#
##Maggies
p <-
  ggplot(data=v4, aes(chromosome, size)) +
  geom_bar(stat="identity", fill="grey70")

p + geom_segment(data=v4_subg, aes(x=Chromosome-0.45, xend=Chromosome+0.45, y=Pos_start, yend=Pos_end, colour=Subgenome), size=3)

#--------------------------------------------------------------------------------------------------#
##Jesses
p <-
  ggplot(data=v4, aes(chromosome, size)) +
    geom_bar(stat="identity", fill="grey70")

p + geom_segment(data=v4_subg_jw, aes(x=Chromosome-0.45, xend=Chromosome+0.45, y=Pos_start, yend=Pos_end, colour=Subgenome), size=3)
#--------------------------------------------------------------------------------------------------#
##Jesses dominance
p <-
  ggplot(data=v4, aes(chromosome, size)) +
  geom_bar(stat="identity", fill="grey70")

p + geom_segment(data=v4_dom_jw, aes(x=Chromosome-0.45, xend=Chromosome+0.45, y=Pos_start, yend=Pos_end, colour=Subgenome), size=3)
#--------------------------------------------------------------------------------------------------#
##Multi
# df1 <- v4_subg
# df1$x <- -.45
# df1$xend <- -.17
# df1$Subgenome[df1$Subgenome %in% c("Maize1")] <- "sub1"
# df1$Subgenome[df1$Subgenome %in% c("Maize2")] <- "sub2"
# df1 <- v4_orph_jw
# df1$x <- -.45
# df1$xend <- -.17
# df1 <- centromere.v4.positions
# df1$x <- -.45
# df1$xend <- -.17
# df2 <- v4_subg_jw
# df2$x <- -.15
# df2$xend <- .15
# df3 <- v4_dom_jw
# df3$x <- .17
# df3$xend <- .45
df1 <- centromere.v4.positions
df1$x <- -.45
df1$xend <- .45
df2 <- v4_subg_jw
df2$x <- -.45
df2$xend <- -.05
df3 <- v4_dom_jw
df3$x <- .05
df3$xend <- .45
segment_data <- bind_rows(df2,df3)
segment_data <- segment_data[!is.na(segment_data$Pos_start),]
segment_data <-
  segment_data %>%
  ungroup()
segment_data <- segment_data %>% subset(Subgenome %in% c("sub1", "sub2", "orphan", "sub1dom", "sub2dom", "dom"))
centromere_data <- df1

#--------------------------------------------------------------------------------------------------#
## Make the plot
p <-
  ggplot(data=v4) +
    geom_bar(aes(v4$chromosome, v4$size), stat="identity", fill="grey70") +
    # geom_segment(data=segment_data, aes(x=Chromosome+x, xend=Chromosome+xend, y=Pos_start, yend=Pos_end, colour=Subgenome), size=.001) +
    geom_rect(data=segment_data, inherit.aes = F, aes(xmin=Chromosome+x, xmax=Chromosome+xend, ymin=Pos_start, ymax=Pos_end, color=Subgenome, fill=Subgenome)) +
    geom_rect(data=centromere_data, inherit.aes = F, aes(xmin=Chromosome+x, xmax=Chromosome+xend, ymin=Pos_start, ymax=Pos_end)) +
    scale_fill_manual(values=c("green", "red", "blue", "black")) +
    scale_color_manual(values=c("green", "red", "blue", "black")) +
    scale_x_discrete(limits=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10")) +
    scale_y_continuous(labels=format_si())
p  +
  labs(
  title="Subgenome Locations by Chromosome",
  x="Chromosome",
  y="Position on Chromosome"
)

#--------------------------------------------------------------------------------------------------#
## Try to print
# ggplot2::ggsave("gplot.png", p, width = 10, height=10, units = "in")
ggsave("gplot.png", width=5, height=4, dpi=600)

png(file="gplot.png",width=200,height=160,res=72)
plot(p+labs(
         title="Subgenome Locations by Chromosome",
         x="Chromosome",
         y="Position on Chromosome"))
dev.off()

#--------------------------------------------------------------------------------------------------#
# Add SI units
format_si <- function(...) {
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html

  function(x) {
    limits <- c(1e0, 1e3, 1e6, 1e9, 1e12)
    prefix <- c("b", "kb", "Mb", "Gb", "Tb")
    # prefix <- c("y",   "z",   "a",   "f",   "p",
    #             "n",   "µ",   "m",   " ",   "k",
    #             "M",   "G",   "T",   "P",   "E",
    #             "Z",   "Y")

    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)

    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)

    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...),
          prefix[i])
  }
}
