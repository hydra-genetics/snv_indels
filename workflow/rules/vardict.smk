# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule vardict:
    input:
        bam="alignment/extract/{sample}_{type}_{chr}.bam",
        reference=config["reference"]["fasta"],
        bed="misc/bed_split/{sample}_{type}_{chr}.bed",
    output:
        vcf=temp("snv_indels/vardict/{sample}_{type}_{chr}.unfilt.vcf.gz"),
    params:
        extra=config.get("vardict", {}).get("extra", "")
        bed_columns=config.get("vardict", {}).get("bed_columns", ""),
        af_th=config.get("vardict", {}).get("af_th", ""),
    log:
        "snv_indels/vardict/{sample}_{type}_{chr}.log",
    benchmark:
        repeat(
            "snv_indels/vardict/{sample}_{type}_{chr}.benchmark.tsv", config.get("vardict", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("vardict", config["default_resources"])["threads"]
    container:
        config.get("vardict", {}).get("container", config["default_container"])
    conda:
        "../envs/vardict.yaml"
    message:
        "{rule}: Use vardict to call variants, snv_indels/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}"
    wrapper:
        "0.78.0/bio/vardict"
