# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule ensemble_vcf:
    input:
        vcfs=expand(
            "snv_indels/{caller}/{{sample}}_{{type}}.normalized.sorted.vcf.gz",
            caller=config.get("ensemble_vcf", {}).get("callers", []),
        ),
        tabix=expand(
            "snv_indels/{caller}/{{sample}}_{{type}}.normalized.sorted.vcf.gz.tbi",
            caller=config.get("ensemble_vcf", {}).get("callers", []),
        ),
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vcf.gz"),
    params:
        support=config.get("ensemble_vcf", {}).get("support", "1"),
        sort_order=config.get("ensemble_vcf", {}).get("sort_order", ""),
    log:
        "snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/ensemble_vcf/{sample}_{type}.ensembled.vcf.gz.benchmark.tsv",
            config.get("ensemble_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ensemble_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("ensemble_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ensemble_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("ensemble_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ensemble_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("ensemble_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("ensemble_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/ensemble_vcf.yaml"
    message:
        "{rule}: Ensamble vcfs from different callers using recall resulting in {wildcards.sample}_{wildcards.type}.ensembled.vcf.gz"
    shell:
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.sort_order} {output.vcf} {input.ref} {input.vcfs}) &> {log}"
