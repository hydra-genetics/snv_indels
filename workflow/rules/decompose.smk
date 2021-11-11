# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule decompose:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.okAF.vcf.gz",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz"),
    log:
        "snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz.benchmark.tsv",
            config.get("decompose", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("decompose", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("decompose", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("decompose", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("decompose", {}).get("container", config["default_container"])
    conda:
        "../envs/decompose.yaml"
    message:
        "{rule}: Decompose vcf snv_indels/{wildcards.caller}/{wildcards.sample}_{wildcards.type}.format_filt.vcf.gz"
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output.vcf} -) &> {log}"
