__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepmosaic_input:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        sex="qc/peddy/peddy.sex_check.csv",
        vcf="snv_indels/deepsomatic/{sample}_{type}.vcf.gz",
    output:
        txt=temp("snv_indels/deepmosaic/{sample}_{type}.input.txt"),
    params:
        extra=config.get("deepmosaic_input", {}).get("extra", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
    log:
        temp("snv_indels/deepmosaic/{sample}_{type}.input.log"),
    benchmark:
        repeat(
            "snv_indels/deepmosaic/{sample}_{type}.input.benchmark.tsv",
            config.get("deepmosaic_input", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_input", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_input", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_input", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_input", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_input", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_input", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_input", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Creates an input file {output.txt} for DeepMosaic with info: {input.bam}, {input.vcf}, depth and sex
        """
    shell:
        "python workflow/scripts/deepmosaic_input.py {input.bam} {input.vcf} {output.txt} {params.name} {input.sex}"


rule deepmosaic_draw:
    input:
        annovar=config.get("reference", {}).get("annovar", ""),
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        txt="snv_indels/deepmosaic/{sample}_{type}.input.txt",
        vcf="snv_indels/deepsomatic/{sample}_{type}.vcf.gz",
    output:
        outdir=temp(directory("snv_indels/deepmosaic/{sample}_{type}/")),
        txt=temp("snv_indels/deepmosaic/{sample}_{type}/features.txt"),
    params:
        extra=config.get("deepmosaic_draw", {}).get("extra", ""),
    log:
        "snv_indels/deepmosaic/{sample}_{type}.draw.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic/{sample}_{type}.draw.benchmark.tsv",
            config.get("deepmosaic_draw", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_draw", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_draw", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_draw", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_draw", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_draw", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_draw", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_draw", {}).get("container", config["default_container"])
    message:
        """
        {rule}: DeepMosaic draw takes input from {input.vcf} and {input.bam} to find mosaic variants
        """
    shell:
        """deepmosaic-draw \
        -i {input.txt} \
        -o {output.outdir} \
        -a {input.annovar} \
        -b hg38 \
        -db gnomad41_genome \
        &> {log}"""


rule deepmosaic_predict:
    input:
        txt="snv_indels/deepmosaic/{sample}_{type}/features.txt",
    output:
        txt=temp("snv_indels/deepmosaic/{sample}_{type}/final_predictions.txt"),
    params:
        extra=config.get("deepmosaic_predict", {}).get("extra", ""),
    log:
        "snv_indels/deepmosaic/{sample}_{type}.predict.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic_predict/{sample}_{type}.predict.benchmark.tsv",
            config.get("deepmosaic_predict", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_predict", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_predict", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_predict", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_predict", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_predict", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_predict", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_predict", {}).get("container", config["default_container"])
    message:
        """
        {rule}:  DeepMosaic predict takes output from draw {input.txt} to predict mosaic variants {output.txt}.
        """
    shell:
        """deepmosaic-predict \
        -i {input.txt} \
        -o {output.txt} \
        -gb hg38 \
        -b 10 \
        &> {log}"""

