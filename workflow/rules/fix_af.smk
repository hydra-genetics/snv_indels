# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule fix_af:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.merged.vcf",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.okAF.vcf"),
    log:
        "snv_indels/{caller}/{sample}_{type}.okAF.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.okAF.vcf.benchmark.tsv",
            config.get("fix_af", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fix_af", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("fix_af", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_af", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("fix_af", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_af", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_af", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("fix_af", {}).get("container", config["default_container"])
    conda:
        "../envs/fix_af.yaml"
    message:
        "{rule}: Fix missing AF in format in snv_indels/{wildcards.caller}/{wildcards.sample}_{wildcards.type}.merged.vcf"
    script:
        "../scripts/fix_af.py"
