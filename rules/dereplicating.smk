# dereplication dada2


rule dereplicate:
    input:
        lambda wildcards: os.path.join(config["trimmed_fastq_dir"], f"barcode{wildcards.barcode}.fastq.gz")
    output: 
        "dereplicated/barcode{barcode}.rds"
        #"dereplicated/barcode{barcode}.fastq.gz"
    conda:
        "../envs/dada2.yaml"
    script:
        "scripts/dereplicate.R"