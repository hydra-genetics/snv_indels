# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule bed_split:
    input:
        bed=config["reference"]["design_bed"],
    output:
        "snv_indels/bed_split/design_bedfile_{chr}.bed",
    log:
        "snv_indels/bed_split/bed_split.log",
    benchmark:
        repeat("snv_indels/bed_split/bed_split.benchmark.tsv", config.get("bed_split", {}).get("benchmark_repeats", 1))
    threads: config.get("bed_split", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("bed_split", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bed_split", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("bed_split", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bed_split", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bed_split", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("bed_split", {}).get("container", config["default_container"])
    conda:
        "../envs/bed_split.yaml"
    message:
        "{rule}: Split design bed file into chromosomes"
    shell:
        "(awk '{{if(/^{wildcards.chr}\t/) print($0)}}' {input.bed}  > {output}) &> {log}"
