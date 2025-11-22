
ruleorder:  # just to emphasize model before basecall
    dorado_model > dorado_basecall

# Download/ensure model
rule dorado_model:
    output:
        # a marker file
        touch(f"models/{config['model']}.ready")
    params:
        model = config["model"],
        dorado = config["dorado_software"]
    shell:
        # module load dorado/0.9.1-foss-2023a-CUDA-12.1.1
        r"""
        mkdir -p models
        "{params.dorado}" download --model {params.model} --directory models
        touch {output}
        """

# Basecalling from POD5 directory to BAM
rule dorado_basecall:
    input:
        reads_dir = config["reads_dir"],
        model_ok = f"models/{config['model']}.ready"
    output:
	# test run, delete test/ later
        basecalls = f"{config['output_dir']}/dorado/basecalls.bam"
    threads: 4
    #resources: # test without slurm/gpu
        #gpus = 1
    #envmodules:
	#"CUDA/12.1.1" # ran module load CUDA before hand
    params:
        model = config["model"],
        dorado = config["dorado_software"],
	device = config.get("dorado_device", "cpu"),
	models_dir = "models"
    shell:
        r"""
        mkdir -p $(dirname {output.basecalls})
        "{params.dorado}" basecaller {params.model} --device {params.device} --models-directory {params.models_dir} --no-trim {input.reads_dir} > {output.basecalls}
        """
	# deleted --emit-bam, should be default
	# --models-directory {params.model_dir}


# summary for basecalling        
rule dorado_basecall_summary:
    input:
        bam = rules.dorado_basecall.output.basecalls
    output:
	# test run, delete test/ later
        summary = "test/results/dorado/basecall_summary.tsv"
    params:
        dorado = config["dorado_software"]
    shell:
        r"""
        "{params.dorado}" summary {input.bam} > {output.summary}
        """

# demultiplex BAM to per-barcode Fastq
rule dorado_demultiplex:
    input:
        bam = rules.dorado_basecall.output.basecalls
    output:
	# test run, delete test/ later
        demux_dir = directory(f"{config['output_dir']}/dorado/demux_fastq")
    params:
        kit = config["kit"],
        dorado = config["dorado_software"]
    threads: 2
    shell:
        r"""
        "{params.dorado}" demux --kit-name {params.kit} --emit-summary --emit-fastq --output-dir {output.demux_dir} {input.bam}
        """

    # add cd {output.dir}
    #     rename 's/barcode([0-9])_/barcode0$1_/' *.fastq.gz
    # for renaming the files with single digits into two digits (9 -> 09)
