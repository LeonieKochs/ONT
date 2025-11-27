# Snakefile

from snakemake.io import directory
import os

configfile: "config/config.yaml"
READS_DIR     = config["reads_dir"]
OUTPUT_DIR    = config["output_dir"]
TRIMMED_DIR   = config["trimmed_fastq_dir"]

# If demux_fastq_dir is empty in config, default to OUTPUT_DIR/dorado/demux_fastq
DEMUX_DIR     = config.get("demux_fastq_dir") or os.path.join(OUTPUT_DIR, "dorado", "demux_fastq")

DORADO_BIN    = config["dorado_software"]
DORADO_MODEL  = config["model"]
DORADO_KIT    = config["kit"]
PRIMER_SETS   = config["primer_sets"]

os.makedirs("logs/slurm", exist_ok=True)

BARCODES_ALL = (
    config["primer_sets"]["set1"]["barcodes"]
    + config["primer_sets"]["set2"]["barcodes"]
)

include: "rules/basecalling.smk"
include: "rules/trimming.smk"
include: "rules/dereplicating.smk"

rule all:
    input:
        # basecalling + summary + demux
        #"results/dorado/basecalls.bam",
        #"results/dorado/basecall_summary.tsv",
        #"results/dorado/demux_fastq",
        expand("trimmed/barcode{barcode}.fastq.gz", barcode=BARCODES_ALL),
        expand("dereplicated/barcode{barcode}.fastq.gz", barcode=BARCODES_ALL)


