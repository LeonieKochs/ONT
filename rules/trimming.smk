
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
        demux_dir = DEMUX_DIR
        #demux_dir = rules.dorado_demultiplex.output.demux_dir # change in config file or add sample id etc
    output:
        os.path.join(TRIMMED_DIR, "barcode{barcode}.fastq.gz")
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
        set -euo pipefail

        LOGDIR="{TRIMMED_DIR}/logs"
        mkdir -p {TRIMMED_DIR} "$LOGDIR"

        demux_root="{input.demux_dir}"

        # demux_root can be:
        #  - the directory that CONTAINS fastq_pass (amplicons2-style),
        #  - OR fastq_pass itself (amplicons1-style)
        if [ "$(basename "$demux_root")" = "fastq_pass" ]; then
            run_fastq_pass="$demux_root"
        elif [ -d "$demux_root/fastq_pass" ]; then
            run_fastq_pass="$demux_root/fastq_pass"
        else
            run_fastq_pass=$(find "$demux_root" -maxdepth 8 -type d -name fastq_pass | head -n 1)
        fi
        
        if [ -z "${run_fastq_pass:-}" ] || [ ! -d "$run_fastq_pass" ]; then
            echo "ERROR: No fastq_pass directory found (or provided) at/under $demux_root" >&2
            exit 1
        fi


        fq_dir="$run_fastq_pass/barcode{wildcards.barcode}"

        if [ ! -d "$fq_dir" ]; then
            echo "No directory for barcode {wildcards.barcode} in $run_fastq_pass - creating empty trimmed file."
            gzip -c </dev/null > {output}
            exit 0
        fi

        # if there are no FASTQ files (fastq or fastq.gz), also produce an empty trimmed file
        if ! ls "$fq_dir"/*.fastq* >/dev/null 2>&1; then
            echo "No FASTQ files for barcode {wildcards.barcode} in $fq_dir - creating empty trimmed file."
            gzip -c </dev/null > {output}
            exit 0
        fi

        cutadapt -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -q 10,10 \
          -o {output} "$fq_dir"/*.fastq* \
          > "$LOGDIR/barcode{wildcards.barcode}_cutadapt.log"
        """

# Rule for barcodes in set2
rule trim_primers_set2:
    conda: "../envs/cutadapt.yaml"
    input:
        demux_dir = DEMUX_DIR
        # demux_dir = rules.dorado_demultiplex.output.demux_dir
    output:
        os.path.join(TRIMMED_DIR, "barcode{barcode}.fastq.gz")
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
        set -euo pipefail

        LOGDIR="{TRIMMED_DIR}/logs"
        mkdir -p {TRIMMED_DIR} "$LOGDIR"

        demux_root="{input.demux_dir}"

        # demux_root can be:
        #  - the directory that CONTAINS fastq_pass (amplicons2-style),
        #  - OR fastq_pass itself (amplicons1-style)
        if [ "$(basename "$demux_root")" = "fastq_pass" ]; then
            run_fastq_pass="$demux_root"
        elif [ -d "$demux_root/fastq_pass" ]; then
            run_fastq_pass="$demux_root/fastq_pass"
        else
            run_fastq_pass=$(find "$demux_root" -maxdepth 8 -type d -name fastq_pass | head -n 1)
        fi
        
        if [ -z "${run_fastq_pass:-}" ] || [ ! -d "$run_fastq_pass" ]; then
            echo "ERROR: No fastq_pass directory found (or provided) at/under $demux_root" >&2
            exit 1
        fi


        fq_dir="$run_fastq_pass/barcode{wildcards.barcode}"

        if [ ! -d "$fq_dir" ]; then
            echo "No directory for barcode {wildcards.barcode} in $run_fastq_pass - creating empty trimmed file."
            gzip -c </dev/null > {output}
            exit 0
        fi

        if ! ls "$fq_dir"/*.fastq* >/dev/null 2>&1; then
            echo "No FASTQ files for barcode {wildcards.barcode} in $fq_dir - creating empty trimmed file."
            gzip -c </dev/null > {output}
            exit 0
        fi

        cutadapt \
          -g {params.fwd} -a {params.rev} \
          --discard-untrimmed \
          -m {params.minlen} -M {params.maxlen} \
          -q 10,10 \
          -o {output} \
          "$fq_dir"/*.fastq* \
          > "$LOGDIR/barcode{wildcards.barcode}_cutadapt.log"
        """



# add histogram
rule length_stats:
    conda: "../envs/cutadapt.yaml"
    input:
        os.path.join(TRIMMED_DIR, "{s}.fastq.gz")
    output:
        "stats/{s}_lengths.txt"
    threads: 1
    shell:
        r"""
        mkdir -p stats
        seqkit stats {input} > {output}
        """
