# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule sort_vcf:
    input:
        vcf="snv_indels/{caller}/{file}.vcf.gz",
        vcf_indexed="snv_indels/{caller}/{file}.vcf.gz.tbi",
    output:
        vcf=temp("snv_indels/{caller}/{file}.sorted.vcf.gz"),
    log:
        "snv_indels/{caller}/{file}.sorted.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{file}.sorted.vcf.gz.benchmark.tsv",
            config.get("sort_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sort_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("sort_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/sort_vcf.yaml"
    message:
        "{rule}: Sort vcf snv_indels/{wildcards.caller}/{wildcards.file}.vcf.gz"
    wrapper:
        "0.79.0/bio/bcftools/sort"
