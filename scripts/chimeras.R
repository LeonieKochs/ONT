library(dada2)

seqtab_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]

seqtab <- readRDS(seqtab_file)

# with ONT full-length reads, chimera behavior can be a bit different!
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method="consensus",
  multithread=TRUE,
  verbose=TRUE
)

saveRDS(seqtab.nochim, out_file)
