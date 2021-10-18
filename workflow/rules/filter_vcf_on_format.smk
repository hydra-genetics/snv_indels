# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule filter_vcf_on_format:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.merged.vcf.gz",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.format_filt.vcf.gz"),
    params:
        filter=config.get("filter_vcf_on_format", {}).get("filter", "-i 'FORMAT/DP > 100 & FORMAT/AD > 20 & FORMAT/AF > 0.05'"),
        extra=config.get("filter_vcf_on_format", {}).get("extra", "--mode '+' --soft-filter 'DP_AD_AF'"),
    log:
        "snv_indels/{caller}/{sample}_{type}.format_filt.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.format_filt.vcf.gz.benchmark.tsv",
            config.get("filter_vcf_on_format", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("filter_vcf_on_format", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("filter_vcf_on_format", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf_on_format.yaml"
    message:
        "{rule}: Filter vcf snv_indels/{wildcards.caller}/{wildcards.sample}_{wildcards.type}.unfilt.merged.vcf.gz based on format"
    wrapper:
        "0.72.0/bio/bcftools/filter"
