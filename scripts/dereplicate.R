# dada2
library(dada2)

fn <- snakemake@input[[1]]
out <- snakemake@output[[1]]

drp <- derepFastq(fn)
saveRDS(drp, out)
