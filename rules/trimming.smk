
import os
#from glob import glob

SET1 = config["primer_sets"]["set1"]
SET2 = config["primer_sets"]["set2"]
BARCODES_SET1 = SET1["barcodes"]
BARCODES_SET2 = SET2["barcodes"]


# Primer trimming
# Rule for barcodes in set1
rule trim_primers_set1:
    conda: "../envs/cutadapt.yaml"
    input: 
        os.path.join(DEMUX_DIR, barcode{barcode}", FBD92602_pass_barcode{barcode}_d57d61d8_00000000_0.fastq.gz") # change in config file or add sample id etc
    output:
        os.path.join(TRIMMED_DIR, barcode{wc.barcode}.fastq.gz")
    params:
        fwd = SET1["fwd"],
        rev = SET1["rev"],
        minlen = SET1["minlen"],
        maxlen = SET1["maxlen"]
    threads: 2
    # Limit this rule to the configured barcodes
    wildcard_constraints:
        # joint to one string 01|02...
        barcode="|".join(BARCODES_SET1)
    shell:
        # Regular 3' adapter -a
        # Regular 5' adapter -g 
        # -m discard processed reads that are shorter than LENGTH
        # -M discard processed reads that are longer than LENGTH
        # --discard-untrimmed discard reads in which no adapter was found 
        r"""
        mkdir -p {TRIMMED_DIR} logs
        cutadapt -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -q 10,10
          -o {output} {input} > logs/barcode{wildcards.barcode}_cutadapt.log
        """

# Rule for barcodes in set2
rule trim_primers_set2:
    conda: "../envs/cutadapt.yaml"
    input:
        os.path.join(DEMUX_DIR, sample_barcode{barcode}.fastq.gz")
    output:
        os.path.join(TRIMMED_DIR, barcode{barcode}.fastq.gz")
    params:
        fwd = SET2["fwd"],
        rev = SET2["rev"],
        minlen = SET2["minlen"],
        maxlen = SET2["maxlen"]
    threads: 2
    wildcard_constraints:
        barcode="|".join(BARCODES_SET2)
    shell:
        r"""
        mkdir -p trimmed logs
        cutadapt -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -q 10,10 \
          -o {output} {input} > logs/barcode{wildcards.barcode}_cutadapt.log
        """



# add histogram
rule length_stats:
    conda: "../envs/cutadapt.yaml"
    input:
        os.path.join(TRIMMED_DIR, {.s}.fastq.gz") # what's the input?
    output:
        "stats/{s}_lengths.txt"
    threads: 1
    shell:
        r"""
        mkdir -p stats
        seqkit stats {input} > {output}
        """
