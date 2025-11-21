library(dada2)

dada_files <- snakemake@input
out_file <- snakemake@output[[1]]

# Load denoised dada objects
dadas <- lapply(dada_files, readRDS)

# name list entries by barcode from filename
get_sample_name <- function(x) sub("\\.rds$", "", basename(x))
names(dadas) <- vapply(dada_files, get_sample_name, character(1))

# Build sequence table:
# rows = samples (barcodes)
# cols = ASVs
seqtab <- makeSequenceTable(dadas)

saveRDS(seqtab, out_file)
