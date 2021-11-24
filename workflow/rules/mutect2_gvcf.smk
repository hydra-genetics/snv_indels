# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutect2_gvcf:
    input:
        map="alignment/mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config["reference"]["fasta"],
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        stats=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.gvcf.gz.stats"),
        vcf=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.gvcf.gz"),
        vcf_tbi=temp("snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.gvcf.gz.tbi"),
    params:
        extra=config.get("mutect2_gvcf", {}).get("extra", "")
        + " --intervals snv_indels/bed_split/design_bedfile_{chr}.bed -ERC BP_RESOLUTION ",
    log:
        "snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.log",
    benchmark:
        repeat(
            "snv_indels/mutect2_gvcf/{sample}_{type}_{chr}.benchmark.tsv",
            config.get("mutect2_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutect2_gvcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("mutect2_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutect2_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutect2_gvcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("mutect2_gvcf", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2_gvcf.yaml"
    message:
        "{rule}: Use mutect2 to generate a gvcf, snv_indels/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}"
    wrapper:
        "0.78.0/bio/gatk/mutect"
