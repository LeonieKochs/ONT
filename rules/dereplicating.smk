# dereplication dada2
import os

rule dereplicate:
    input:
        os.path.join(TRIMMED_DIR, "barcode{barcode}.fastq.gz")
    output:
        os.path.join(OUTPUT_DIR, "dada2", "derep", "barcode{barcode}.rds")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/dereplicate.R"


rule learn_errors:
    input:
        expand(os.path.join(OUTPUT_DIR, "dada2", "derep", "barcode{barcode}.rds"), barcode=BARCODES_ALL)
    output:
        os.path.join(OUTPUT_DIR, "dada2", "errors.rds")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/errors.R"


rule dada:
    input:
        derep = os.path.join(OUTPUT_DIR, "dada2", "derep", "barcode{barcode}.rds"),
        err   = os.path.join(OUTPUT_DIR, "dada2", "errors.rds")
    output:
        os.path.join(OUTPUT_DIR, "dada2", "dd", "barcode{barcode}.rds")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/dada.R"


rule make_seqtab:
    input:
        expand(os.path.join(OUTPUT_DIR, "dada2", "dd", "barcode{barcode}.rds"), barcode=BARCODES_ALL)
    output:
        os.path.join(OUTPUT_DIR, "dada2", "seqtab.rds")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/seqtab.R"


rule remove_chimeras:
    input:
        os.path.join(OUTPUT_DIR, "dada2", "seqtab.rds")
    output:
        os.path.join(OUTPUT_DIR, "dada2", "seqtab.nochim.rds")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/chimeras.R"


rule summary_plots:
    input:
        seqtab = os.path.join(OUTPUT_DIR, "dada2", "seqtab.nochim.rds")
    output:
        track = os.path.join(OUTPUT_DIR, "dada2", "track.tsv"),
        chao  = os.path.join(OUTPUT_DIR, "dada2", "chao_curves.pdf")
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/summary_plots.R"
