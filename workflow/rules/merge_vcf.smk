# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_vcf:
    input:
        calls=expand(
            "snv_indels/{{caller}}/{{sample}}_{{type}}_{chr}.unfilt.vcf.gz",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("merge_vcf", {}).get("skip_chrs", [])
            ),
        ),
    output:
        temp("snv_indels/{caller}/{sample}_{type}.unfilt.merged.vcf.gz"),
    params:
        extra=config.get("merge_vcf", {}).get("extra", ""),
    log:
        "snv_indels/{caller}/{sample}_{type}.log",
    benchmark:
        repeat("snv_indels/{caller}/{sample}_{type}.benchmark.tsv", config.get("vardict", {}).get("benchmark_repeats", 1))
    threads: config.get("merge_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("merge_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/merge_vcf.yaml"
    message:
        "{rule}: Merge vcf on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}"
    wrapper:
        "0.79.0/bio/bcftools/concat"
