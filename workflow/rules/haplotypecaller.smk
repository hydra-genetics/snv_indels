__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule haplotypecaller:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config["reference"]["fasta"],
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        vcf=temp("snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf"),
    params:
        extra=config.get("haplotypecaller", {}).get("extra", ""),
        java_opts=config.get("haplotypecaller", {}).get("java_opts", ""),
    log:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.log",
    benchmark:
        repeat(
            "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.benchmark.tsv",
            config.get("haplotypecaller", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("haplotypecaller", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("haplotypecaller", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("haplotypecaller", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("haplotypecaller", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("haplotypecaller", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("haplotypecaller", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("haplotypecaller", {}).get("container", config["default_container"])
    conda:
        "../envs/haplotypecaller.yaml"
    message:
        "{rule}: call variants in {wildcards.chr} in {input.bam}"
    shell:
        "(gatk --java-options '{params.java_opts}' HaplotypeCaller "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-O {output.vcf} "
        "-L {input.bed} "
        "-A AlleleFraction "
        "{params.extra}) &> {log}"
