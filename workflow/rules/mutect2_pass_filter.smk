__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutect2_pass_filter:
    input:
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.merged.softfitered.vcf.gz",
    output:
        vcf=temp("snv_indels/gatk_mutect2/{sample}_{type}.merged.vcf"),
    log:
        "snv_indels/gatk_mutect2/{sample}_{type}.merged.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2/{sample}_{type}.merged.vcf.gz.benchmark.tsv",
            config.get("mutect2_pass_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mutect2_pass_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mutect2_pass_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutect2_pass_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutect2_pass_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mutect2_pass_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutect2_pass_filter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mutect2_pass_filter", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2_pass_filter.yaml"
    message:
        "{rule}: hardfilter mutect2 variants in {input.vcf} to {output.vcf}"
    script:
        "../scripts/mutect2_pass_filter.py"
