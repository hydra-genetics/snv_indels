# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule normalize:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz"),
    log:
        "snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz.benchmark.tsv",
            config.get("normalize", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("normalize", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("normalize", {}).get("container", config["default_container"])
    conda:
        "../envs/normalize.yaml"
    message:
        "{rule}: Normalize vcf snv_indels/{wildcards.caller}/{wildcards.sample}_{wildcards.type}.decomposed.vcf.gz"
    shell:
        "(vt normalize -n -r {input.ref} -o {output.vcf} {input.vcf} ) &> {log}"
