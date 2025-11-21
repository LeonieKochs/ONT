library(dada2)

# Snakemake passes a list of derep .rds files
derep_files <- snakemake@input
out_file <- snakemake@output[[1]]

# Load dereplication objects into a list
dereps <- lapply(derep_files, readRDS)

# Name them from filenames (e.g. barcode01, barcode02)
get_sample_name <- function(x) {
  nm <- sub("\\.rds$", "", basename(x))
  nm <- sub("^barcode", "barcode", nm)  # keeps barcodeXX
  nm
}
names(dereps) <- vapply(derep_files, get_sample_name, character(1))

# Learn a single error model across all samples (?)
# DADA2 supports learnErrors() on derep-class objects.
err <- learnErrors(dereps, multithread=TRUE)

# Save model
saveRDS(err, out_file)
