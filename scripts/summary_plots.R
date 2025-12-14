library(ggplot2)
library(vegan)
library(tidyverse)

seqtab_file <- snakemake@input[["seqtab"]]

track_out <- snakemake@output[["track"]]
chao_out <- snakemake@output[["chao"]]
asv_wide_out  <- snakemake@output[["asv_wide"]]
asv_long_out  <- snakemake@output[["asv_long"]]
asv_fasta_out <- snakemake@output[["asv_fasta"]]
rare_out      <- snakemake@output[["rarecurve"]]
accum_out     <- snakemake@output[["accum"]]

dir.create(dirname(track_out), recursive = TRUE, showWarnings = FALSE) # check for all dir?

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
p <- ggplot(track, aes(x = reads_nonchim)) +
  geom_point(
    aes(y = observed_asvs, color = "Observed ASVs"),
    size = 2,
    alpha = 0.7
  ) +
  geom_point(
    aes(y = chao1, color = "Chao1"),
    size = 2,
    shape = 1,
    alpha = 0.7
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed ASVs" = "#1f78b4",  # soft blue
      "Chao1"         = "#33a02c"   # soft green
    )
  ) +
  labs(
    x = "Reads after chimera removal",
    y = "Richness",
    title = "Observed ASVs and Chao1 vs sequencing depth",
    subtitle = "Filled points = Observed ASVs, open circles = Chao1"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank()
  )

ggsave(chao_out, p, width = 7, height = 5)

# -------------------------
# 2) ASV tables + sequences
# Original ASV sequences are the column names *before* renaming to ASV IDs
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_len(ncol(seqtab.nochim)))

# wide table uses ASV IDs as column names
seqtab_asv <- seqtab.nochim
colnames(seqtab_asv) <- asv_ids

wide_tab <- as.data.frame(seqtab_asv)
wide_tab$sample <- rownames(wide_tab)
wide_tab <- wide_tab[, c("sample", asv_ids)]

write.table(wide_tab, file = asv_wide_out, sep = "\t", quote = FALSE, row.names = FALSE)

# long table
# (no extra packages needed)
long_list <- lapply(seq_len(nrow(seqtab_asv)), function(i) {
  counts <- seqtab_asv[i, ]
  idx <- which(counts > 0)
  if (length(idx) == 0) return(NULL)
  data.frame(
    sample = rownames(seqtab_asv)[i],
    ASV = asv_ids[idx],
    abundance = as.integer(counts[idx]),
    sequence = asv_seqs[idx],
    stringsAsFactors = FALSE
  )
})
long_tab <- do.call(rbind, long_list)
if (is.null(long_tab)) {
  long_tab <- data.frame(sample=character(), ASV=character(), abundance=integer(), sequence=character())
}
write.table(long_tab, file = asv_long_out, sep = "\t", quote = FALSE, row.names = FALSE)

# FASTA (no Biostrings dependency)
con <- file(asv_fasta_out, open = "w")
for (i in seq_along(asv_seqs)) {
  writeLines(paste0(">", asv_ids[i]), con)
  writeLines(asv_seqs[i], con)
}
close(con)

# -------------------------
# 3) Rarefaction curves

# primer from config for color differentiation
ps <- snakemake@config[["primer_sets"]]

primer_df <- bind_rows(lapply(names(ps), function(set_key) {
  barcs <- ps[[set_key]][["barcodes"]]
  data.frame(
    sample = paste0("barcode", barcs),
    primer = ps[[set_key]][["name"]], # show "515FY-926R" and "DIV4_P5-P7"
    stringsAsFactors = FALSE
  )
}))

# Rarefaction curves (colored by primer set)
rare_list <- lapply(1:nrow(seqtab.nochim), function(i) {
  mat <- seqtab.nochim[i, , drop = FALSE]
  depth <- as.integer(rowSums(mat))

  reads_vec <- seq(1000, depth, by = 5000)
  reads_vec <- reads_vec[reads_vec <= depth]
  
  rare <- vegan::rarefy(mat, sample = reads_vec)
  
  data.frame(
    sample = rownames(mat),
    reads = reads_vec,
    asvs = as.numeric(rare),
    stringsAsFactors = FALSE
  )
})

rare_df <- bind_rows(rare_list)
rare_df <- left_join(rare_df, primer_df, by = "sample")

p_rare <- ggplot(rare_df,
           aes(x = reads, y = asvs, color = primer, group = sample)) +
  geom_line(alpha = 0.7) +
  labs(
    x = "Reads per sample",
    y = "Observed ASVs",
    title = "Rarefaction curves"
  ) +
  theme_bw()

ggsave(rare_out, p_rare, width = 7, height = 5)

# -------------------------
# 4) ASV accumulation across samples plot
# not needed
spec_acc <- vegan::specaccum(seqtab.nochim, method = "random")
pdf(accum_out, width = 7, height = 5)
plot(spec_acc,
     xlab = "Number of samples",
     ylab = "Cumulative ASVs",
     ci.type = "poly",
     ci.col = "grey85",
     col = "#1f78b4",
     lwd = 2)
title("ASV accumulation across samples")
dev.off()
