# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_vcf:
    input:
        calls=expand(
            "snv_indels/{{caller}}/{{sample}}_{{type}}_{chr}.vcf.gz",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
        calls_index=expand(
            "snv_indels/{{caller}}/{{sample}}_{{type}}_{chr}.vcf.gz.tbi",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
    output:
        temp("snv_indels/{caller}/{sample}_{type}.merged.vcf"),
    params:
        extra=config.get("merge_vcf", {}).get("extra", ""),
    log:
        "snv_indels/{caller}/{sample}_{type}.merged.vcf.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.merged.vcf.benchmark.tsv",
            config.get("merge_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("merge_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("merge_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge_vcf", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("merge_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge_vcf", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("merge_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/merge_vcf.yaml"
    message:
        "{rule}: Merge vcf on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}"
    wrapper:
        "0.79.0/bio/bcftools/concat"
