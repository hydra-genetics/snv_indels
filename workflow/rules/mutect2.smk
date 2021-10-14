# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutect2:
    input:
        map="alignment/extract/{sample}_{type}_{chr}.bam",
        bai="alignment/extract/{sample}_{type}_{chr}.bam.bai",
        fasta=config["reference"]["fasta"],
        bed="misc/bed_split/{sample}_{type}_{chr}.bed",
    output:
        bam=temp("snv_indels/mutect2/{sample}_{type}_{chr}.bam"),
        bai=temp("snv_indels/mutect2/{sample}_{type}_{chr}.bai"),
        stats=temp("snv_indels/mutect2/{sample}_{type}_{chr}.unfilt.vcf.gz.stats"),
        vcf=temp("snv_indels/mutect2/{sample}_{type}_{chr}.unfilt.vcf.gz"),
        vcf_tbi=temp("snv_indels/mutect2/{sample}_{type}_{chr}.unfilt.vcf.gz.tbi"),
        f1f2=temp("snv_indels/mutect2/{sample}_{type}_{chr}.f1f2.tar.gz"),
    params:
        extra=config.get("mutect2", {}).get("extra", "")
        +" --intervals misc/bed_split/{sample}_{type}_{chr}.bed "
        +"--f1r2-tar-gz snv_indels/mutect2/{sample}_{type}_{chr}.f1f2.tar.gz ",
    log:
        "snv_indels/mutect2/{sample}_{type}_{chr}.log",
    benchmark:
        repeat("snv_indels/mutect2/{sample}_{type}_{chr}.benchmark.tsv", config.get("mutect2", {}).get("benchmark_repeats", 1))
    threads: config.get("mutect2", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("mutect2", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2.yaml"
    message:
        "{rule}: Use mutect2 to call variants, snv_indels/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.input"
    wrapper:
        "0.78.0/bio/gatk/mutect"
