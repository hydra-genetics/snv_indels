__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutect2:
    input:
        map="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config["reference"]["fasta"],
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        bam=temp("snv_indels/mutect2/{sample}_{type}_{chr}.bam"),
        bai=temp("snv_indels/mutect2/{sample}_{type}_{chr}.bai"),
        stats=temp("snv_indels/mutect2/{sample}_{type}_{chr}.vcf.gz.stats"),
        vcf=temp("snv_indels/mutect2/{sample}_{type}_{chr}.vcf.gz"),
        tbi=temp("snv_indels/mutect2/{sample}_{type}_{chr}.vcf.gz.tbi"),
        f1f2=temp("snv_indels/mutect2/{sample}_{type}_{chr}.f1r2.tar.gz"),
    params:
        extra=lambda wildcards: get_mutect2_extra(wildcards, "mutect2"),
    log:
        "snv_indels/mutect2/{sample}_{type}_{chr}.vcf.gz.log",
    benchmark:
        repeat("snv_indels/mutect2/{sample}_{type}_{chr}.vcf.gz.benchmark.tsv", config.get("mutect2", {}).get("benchmark_repeats", 1))
    threads: config.get("mutect2", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mutect2", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutect2", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutect2", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mutect2", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutect2", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mutect2", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2.yaml"
    message:
        "{rule}: call variants in {input.map}"
    wrapper:
        "v1.3.1/bio/gatk/mutect"


rule mutect2_gvcf:
    input:
        map="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config["reference"]["fasta"],
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        stats=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.vcf.gz.stats"),
        vcf=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.vcf.gz"),
        tbi=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.vcf.gz.tbi"),
    params:
        extra=lambda wildcards: get_mutect2_extra(wildcards, "mutect2_gvcf"),
    log:
        "snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.vcf.gz.benchmark.tsv",
            config.get("mutect2_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mutect2_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutect2_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutect2_gvcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutect2_gvcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mutect2_gvcf", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2.yaml"
    message:
        "{rule}: generate gvcf from {input.map}"
    wrapper:
        "0.78.0/bio/gatk/mutect"
