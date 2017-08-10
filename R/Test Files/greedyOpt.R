####################################################################################################
## Function greedOpt ->
##        Use a greedy algorithm to sort chromosomes into sugenome 1 and subgenome 2, where
##        subgenome 1 is the larger (by gene count) region that may include multiple org2 chromosomes,
##        but may not overlap when projected to the org1 chromosome.
##
##  homeologs.chromosome -> must be sorted by gene count
####################################################################################################
# greedOpt <- function(homeologs.chromosome, syntelogs.raw){
#   # match(unique(homeologs.chromosome$chr1), homeologs.chromosome$chr1) %>%
#   # slice(homeologs.chromosome, .) %>%
#   # select(org_chr1, org_chr2) %>%
#   # mutate(subgenome="sub1") %>%
#   # left_join(homeologs.genes, ., by=c("org_chr1", "org_chr2"))
#
#   chromosomes <-
#     homeologs.chromosome %>%
#     select(org_chr2,chromosomeGeneCount) %>%
#     distinct()
#
#   homeologs.chromosomeStopStart <-
#     syntelogs.raw %>%
#     inner_join(homeologs.block, by=c("block", "org_chr1", "org_chr2")) %>%
#     mutate(low = pmin(start1, stop1), high = pmax(start1, stop1)) %>%
#     group_by(org_chr1, org_chr2) %>%
#     summarise(chrLow = min(low,high), chrHigh = max(low,high)) %>%
#     left_join(homeologs.chromosome, ., by=c("org_chr1", "org_chr2")) %>%
#     select(org_chr1, org_chr2, chromosomeGeneCount, chrLow, chrHigh)
#
#   x <- 1
#   while(x < nrow(chromosomes)) {
#     currentChromosome <-
#       homeologs.chromosomeStopStart %>%
#       inner_join(chromosomes[x,], by="org_chr1")
#
#     for (row in 1:nrow(currentChromosome)) {
#       count <- currentChromosome[row, "chromosomeGeneCount"]
#       low <- currentChromosome[row, "chrLow"]
#       high <- currentChromosome[row, "chrHigh"]
#
#     }
#     x <- x+1;
#   }

# N           <- ncol(X)
# weights     <- rep(0L, N)
# pred        <- 0 * X
# sum.weights <- 0L
#
# while(sum.weights < iter) {
#
#   sum.weights   <- sum.weights + 1L
#   pred          <- (pred + X) * (1L / sum.weights)
#   errors        <- sqrt(colSums((pred - Y) ^ 2L))
#   best          <- which.min(errors)
#   weights[best] <- weights[best] + 1L
#   pred          <- pred[, best] * sum.weights
# }
# return(weights / sum.weights)
# }
# greedOpt(homeologs.chromosome)
