# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_gvcf:
    input:
        calls=expand(
            "snv_indels/{{caller}}/{{sample}}_{{type}}_{chr}.gvcf.gz",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
        calls_index=expand(
            "snv_indels/{{caller}}/{{sample}}_{{type}}_{chr}.gvcf.gz.tbi",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
    output:
        temp("snv_indels/{caller}/{sample}_{type}.merged.gvcf.gz"),
    params:
        extra=config.get("merge_gvcf", {}).get("extra", "") + " -O z",
    log:
        "snv_indels/{caller}/{sample}_{type}.merged.gvcf.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.merged.gvcf.benchmark.tsv",
            config.get("merge_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("merge_gvcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("merge_gvcf", {}).get("container", config["default_container"])
    conda:
        "../envs/merge_gvcf.yaml"
    message:
        "{rule}: Merge gvcf on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}"
    wrapper:
        "0.70.0/bio/bcftools/concat"
