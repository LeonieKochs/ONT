
ruleorder:  # just to emphasize model before basecall
    dorado_model > dorado_basecall

# Download/ensure model
rule dorado_model:
    output:
        # a marker file
        touch("models/DORADO_MODEL.ready")
    params:
        model = DORADO_MODEL,
        dorado = DORADO_BIN
    shell:
        # module load dorado/0.9.1-foss-2023a-CUDA-12.1.1 # in shell script?
        r"""
        mkdir -p models
        "{params.dorado}" download --model {params.model} --directory models
        touch {output}
        """

# Basecalling from POD5 directory to BAM
rule dorado_basecall:
    input:
        reads_dir = READS_DIR,
        model_ok = "models/DORADO_MODEL.ready"
    output:
        basecalls=os.path.join(OUTPUT_DIR, "dorado", "basecalls.bam")
    threads: 4
    resources:
        gpus = 1
    #envmodules: # run without conda envs
	    #"CUDA/12.1.1"
    params:
        model = DORADO_MODEL,
        dorado = DORADO_BIN,
	    kit=DORADO_KIT
    shell:
        r"""
        mkdir -p $(dirname {output.basecalls})
        "{params.dorado}" basecaller {params.model} --kit-name {params.kit} --no-trim {input.reads_dir} > {output.basecalls}
        """
	# deleted --emit-bam, should be default
	# --models-directory {params.model_dir}


# summary for basecalling        
rule dorado_basecall_summary:
    input:
        bam = rules.dorado_basecall.output.basecalls
    output:
        summary = os.path.join(OUTPUT_DIR, "dorado", "basecall_summary.tsv")
    params:
        dorado = DORADO_BIN
    shell:
        r"""
        "{params.dorado}" summary {input.bam} > {output.summary}
        """

# demultiplex BAM to per-barcode Fastq
rule dorado_demultiplex:
    input:
        bam = rules.dorado_basecall.output.basecalls
    output:
        demux_dir=directory(DEMUX_DIR)
    params:
        kit=DORADO_KIT,
        dorado = DORADO_BIN
    threads: 2
    shell:
        r"""
        "{params.dorado}" demux --kit-name {params.kit} --emit-summary --emit-fastq --output-dir {output.demux_dir} {input.bam}
        """

    # add cd {output.dir}
    #     rename 's/barcode([0-9])_/barcode0$1_/' *.fastq.gz
    # for renaming the files with single digits into two digits (9 -> 09)
