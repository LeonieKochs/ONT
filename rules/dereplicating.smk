# dereplication dada2
import os

rule dereplicate:
    input:
        lambda w: os.path.join(config["trimmed_fastq_dir"], f"barcode{w.barcode}.fastq.gz")
    output:
        "dada2/derep/barcode{barcode}.rds"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/dereplicate.R"


rule learn_errors:
    input:
        expand("dada2/derep/barcode{barcode}.rds", barcode=BARCODES_ALL)
    output:
        "dada2/errors.rds"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/errors.R"


rule dada:
    input:
        derep="dada2/derep/barcode{barcode}.rds",
        err="dada2/errors.rds"
    output:
        "dada2/dd/barcode{barcode}.rds"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/dada.R"


rule make_seqtab:
    input:
        expand("dada2/dd/barcode{barcode}.rds", barcode=BARCODES_ALL)
    output:
        "dada2/seqtab.rds"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/seqtab.R"


rule remove_chimeras:
    input:
        "dada2/seqtab.rds"
    output:
        "dada2/seqtab.nochim.rds"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/chimeras.R"


rule summary_plots:
    input:
        seqtab="dada2/seqtab.nochim.rds"
    output:
        track="dada2/track.tsv",
        chao="dada2/chao_curves.pdf"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/summary_plots.R"
