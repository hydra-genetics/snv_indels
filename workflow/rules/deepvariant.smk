__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepvariant:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz"),
        gvcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.g.vcf.gz")
        if config.get("deepvariant", {}).get("output_gvcf", False)
        else [],
    params:
        model_type=config.get("deepvariant", {}).get("model_type", ""),
        output_gvcf=lambda wildcards: get_gvcf_output(wildcards, "deepvariant"),
        int_res=lambda wildcards: f"snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}_{wildcards.chr}",
        regions=lambda wildcards: f" --regions {wildcards.chr} ",
        extra=config.get("deepvariant", {}).get("extra", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz.benchmark.tsv",
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
    message:
        "{rule}: Run deepvariant on {input.bam}"
    shell:
        "(run_deepvariant "
        "--model_type {params.model_type} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "{params.regions} "
        "{params.extra} "
        "--output_vcf {output.vcf} "
        "{params.output_gvcf} "
        "--intermediate_results_dir {params.int_res} "
        "--num_shards {threads} ) &> {log}"


rule deepvariant_pacbio:
    input:
        bam="alignment/pbmm2_align/{sample}_{type}.bam",
        bai="alignment/pbmm2_align/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz"),
        gvcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.g.vcf.gz")
        if config.get("deepvariant", {}).get("output_gvcf", False)
        else [],
    params:
        model_type=config.get("deepvariant", {}).get("model_type", "PACBIO"),
        output_gvcf=lambda wildcards: get_gvcf_output(wildcards, "deepvariant"),
        int_res=lambda wildcards: f"snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}_{wildcards.chr}",
        regions=lambda wildcards: f" --regions {wildcards.chr} ",
        extra=config.get("deepvariant", {}).get("extra", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}.vcf.gz.benchmark.tsv",
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
    message:
        "{rule}: Run deepvariant on {input.bam}"
    shell:
        "(run_deepvariant "
        "--model_type {params.model_type} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "{params.regions} "
        "{params.extra} "
        "--output_vcf {output.vcf} "
        "{params.output_gvcf} "
        "--intermediate_results_dir {params.int_res} "
        "--num_shards {threads} ) &> {log}"
