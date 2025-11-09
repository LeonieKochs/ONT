# do I want to trimm for minlen and maxlen? What would be good min and maxlen?
# do I want to discard untrimmed?
# Decide when to gzip !
import os

SET1 = config["primer_sets"]["set1"]
SET2 = config["primer_sets"]["set2"]
BARCODES_SET1 = SET1["barcodes"]
BARCODES_SET2 = SET2["barcodes"]


# Primer trimming
# Rule for barcodes in set1
rule trim_primers_set1:
    conda: "../envs/cutadapt.yaml"
    input:
        os.path.join(config["demux_fastq_dir"], "barcode{barcode}","*.fastq")
    output:
        "trimmed/barcode{barcode}.fastq.gz"
    params:
        fwd = SET1["fwd"],
        rev = SET1["rev"],
        minlen = SET1["minlen"],
        maxlen = SET1["maxlen"]
    threads: 2
    # Limit this rule to the configured barcodes
    wildcard_constraint:
        # joint to one string 01|02...
        barcode="|".join(BARCODES_SET1)
    shell:
        # Regular 3' adapter -a
        # Regular 5' adapter -g 
        # -m discard processed reads that are shorter than LENGTH
        # -M discard processed reads that are longer than LENGTH
        # --discard-untrimmed discard reads in which no adapter was found 
        r"""
        mkdir -p trimmed logs
        cutadapt -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -o {output} {input} > logs/barcode{wildcards.barcode}_cutadapt.log
        """

# Rule for barcodes in set2
rule trim_primers_set2:
    conda: "../envs/cutadapt.yaml"
    input:
        # change like above ?!
        os.path.join(config["demux_fastq_dir"], "sample_barcode{barcode}.fastq.gz")

    output:
        "trimmed/barcode{barcode}.fastq"
    params:
        fwd = SET2["fwd"],
        rev = SET2["rev"],
        minlen = SET2["minlen"],
        maxlen = SET2["maxlen"]
    threads: 2
    wildcard_constraint:
        barcode="|".join(BARCODES_SET2)
    shell:
        r"""
        mkdir -p trimmed logs
        cutadapt -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -o {output} {input} > logs/barcode{wildcards.barcode}_cutadapt.log
        """



# add histogram
rule length_stats:
    input:
        "trimmed/{s}.fastq.gz"
    output:
        "stats/{s}_lengths.txt"
    threads: 1
    shell:
        """
        mkdir -p stats
        seqkit stats {input} > {output}
        """
