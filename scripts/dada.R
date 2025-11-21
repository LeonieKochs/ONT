library(dada2)

derep_file <- snakemake@input[["derep"]]
err_file   <- snakemake@input[["err"]]
out_file   <- snakemake@output[[1]]

drp <- readRDS(derep_file)
err <- readRDS(err_file)

# tune OMEGA_A, band_size, or use pooling!
dada <- dada(drp, err=err, multithread=TRUE)
saveRDS(dada, out_file)
