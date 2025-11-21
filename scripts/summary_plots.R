library(ggplot2)

seqtab_file <- snakemake@input[["seqtab"]]
track_out <- snakemake@output[["track"]]
chao_out <- snakemake@output[["chao"]]

seqtab.nochim <- readRDS(seqtab_file)

# basic per-sample stats
reads_nonchim <- rowSums(seqtab.nochim)
observed_asvs <- rowSums(seqtab.nochim > 0)

# compute Chao1 manually, no package
# no vegan::estimateR(counts), compute directly with formula
# Chao1 = Sobs + (F1^2)/(2*F2)
# F1 = singletons, F2 = doubletons
calc_chao1 <- function(counts) {
  f1 <- sum(counts == 1)
  f2 <- sum(counts == 2)
  sobs <- sum(counts > 0)
  if (f2 == 0) {
    # avoid division by zero
    return(sobs + f1 * (f1 - 1) / 2)
  } else {
    return(sobs + (f1^2) / (2 * f2))
  }
}
chao1 <- apply(seqtab.nochim, 1, calc_chao1)

track <- data.frame(
  sample = rownames(seqtab.nochim),
  reads_nonchim = reads_nonchim,
  observed_asvs = observed_asvs,
  chao1 = chao1
)

# save tracking table
write.table(track, file=track_out, sep="\t", quote=FALSE, row.names=FALSE)

# plot richness vs depth
p <- ggplot(track, aes(x=reads_nonchim)) +
  geom_point(aes(y=observed_asvs), size=2) +
  geom_point(aes(y=chao1), size=2, shape=1) +
  labs(
    x="Reads after chimera removal",
    y="Richness",
    title="Observed ASVs and Chao1 vs sequencing depth",
    subtitle="Filled points = Observed ASVs, open circles = Chao1"
  ) +
  theme_minimal()

ggsave(chao_out, p, width=7, height=5)
