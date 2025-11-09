# Snakefile

from snakemake.io import directory
import os

configfile: "config/config.yaml"
config["reads_dir"] = config["reads_dir"].format(USER=os.getenv("USER"))

os.makedirs("logs/slurm", exist_ok=True)

include: "rules/basecalling.smk"
#include: "rules/trimming.smk"

BARCODES_ALL = (
    config["primer_sets"]["set1"]["barcodes"]
    + config["primer_sets"]["set2"]["barcodes"]
)

rule all:
    input:
        # basecalling + summary + demux
        "results/dorado/basecalls.bam",
        "results/dorado/basecall_summary.tsv",
        "results/dorado/demux_fastq",
        #expand("trimmed/barcode{barcode}.fastq.gz", barcode=BARCODES_ALL)


