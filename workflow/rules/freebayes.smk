# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule freebayes:
    input:
        ref=config["reference"]["fasta"],
        samples="alignment/mark_duplicates/{sample}_{type}_{chr}.bam",
        indexes="alignment/mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        regions="misc/bed_split/{sample}_{type}_{chr}.bed",
    output:
        vcf=temp("snv_indels/freebayes/{sample}_{type}_{chr}.vcf"),
    params:
        extra="--target misc/bed_split/{sample}_{type}_{chr}.bed %s"
        % config.get("freebayes", {}).get("extra", "--min-alternate-fraction 0.01 --genotype-qualities --strict-vcf"),
    log:
        "snv_indels/freebayes/{sample}_{type}_{chr}.unfilt.vcf.log",
    benchmark:
        repeat(
            "snv_indels/freebayes/{sample}_{type}_{chr}.benchmark.tsv", config.get("freebayes", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("freebayes", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("freebayes", {}).get("container", config["default_container"])
    conda:
        "../envs/freebayes.yaml"
    message:
        "{rule}: Use freebayes to call variants, snv_indels/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}"
    wrapper:
        "0.78.0/bio/freebayes"
