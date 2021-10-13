# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule filter_vcf_on_format:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.unfilt.merged.vcf.gz",
    output:
        vcf="snv_indels/{caller}/{sample}_{type}.merged.format_filt.vcf.gz",
    params:
        filter=config.get("filter_vcf_on_format", {}).get("filter", "-i 'FORMAT/DP > 100 & FORMAT/AD > 20 & FORMAT/AF > 0.05'"),
        extra=config.get("filter_vcf_on_format", {}).get("extra", "-m + -s +"),
    log:
        "snv_indels/filter_vcf_on_format/{sample}_{type}.log",
    benchmark:
        repeat(
            "snv_indels/filter_vcf_on_format/{sample}_{type}.benchmark.tsv",
            config.get("filter_vcf_on_format", {}).get("benchmark_repeats", 1)
            )
    threads: config.get(
                    "filter_vcf_on_format", config["default_resources"]).get("threads", config["default_resources"]['threads']
                    )
    container:
        config.get("filter_vcf_on_format", {}).get("container", config["default_container"])
    conda:
        "../envs/filter_vcf_on_format.yaml"
    message:
        "{rule}: Filter vcf snv_indels/{rule}/{wildcards.sample}_{wildcards.type}.merged.format_filt.vcf.gz based on format"
    wrapper:
        "0.78.0/bio/bcftools/filter"
