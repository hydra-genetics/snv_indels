__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepvariant:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="snv_indels/deepvariant/{sample}_{type}_{chr}.vcf",
    params:
        model=config.get("deepvariant", {}).get("model", "WGS"),
        extra=config.get("deepvariant", {}).get("extra", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("deepvariant", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Call variants with deepvariant on {wildcards.sample}_{wildcards.type}"
    shell:
        "run_deepvariant "
        "--model_type={params.model} "
        "--ref={input.ref} "
        "--reads={input.bam} "
        "{params.extra} "
        "--output_vcf={output.vcf} "
        "--num_shards={threads} "


rule deepvariant_gvcf:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf="snv_indels/deepvariant/{sample}_{type}_{chr}.vcf",
        gvcf="snv_indels/deepvariant/{sample}_{type}_{chr}.g.vcf",
    params:
        model=config.get("deepvariant", {}).get("model", "WGS"),
        extra=config.get("deepvariant", {}).get("extra", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("deepvariant", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Call variants with deepvariant on {wildcards.sample}_{wildcards.type}"
    shell:
        "run_deepvariant "
        "--model_type={params.model} "
        "--ref={input.ref} "
        "--reads={input.bam} "
        "{params.extra} "
        "--output_vcf={output.vcf} "
        "--output_gvcf={output.gvcf} "
        "--num_shards={threads} "
