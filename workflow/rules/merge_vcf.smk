# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_vcf:
    input:
        vcf=expand("snv_indels/{caller}/{sample}_{type}_{chr}.unfilt.vcf.gz", chr=wildcards.chr),
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.unfilt.merged.vcf.gz"),
    params:
        extra="-O z " + config.get("merge_vcf", {}).get("extra", "")
    log:
        "snv_indels/{caller}/{sample}_{type}_{chr}.log",
    benchmark:
        repeat("snv_indels/{caller}/{sample}_{type}_{chr}.benchmark.tsv", config.get("vardict", {}).get("benchmark_repeats", 1))
    threads: config.get("merge_vcf", config["default_resources"]).get("threads", config["default_resources"]['threads'])
    container:
        config.get("merge_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/merge_vcf.yaml"
    message:
        "{rule}: Merge vcf on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}"
    wrapper:
        "0.79.0/bio/bcftools/merge"
