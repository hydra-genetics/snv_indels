__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule bed_split:
    input:
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        bed=temp("snv_indels/bed_split/design_bedfile_{chr}.bed"),
    log:
        "snv_indels/bed_split/bed_split_{chr}.bed.log",
    benchmark:
        repeat("snv_indels/bed_split/bed_split_{chr}.bed.benchmark.tsv", config.get("bed_split", {}).get("benchmark_repeats", 1))
    threads: config.get("bed_split", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bed_split", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bed_split", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bed_split", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bed_split", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bed_split", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bed_split", {}).get("container", config["default_container"])
    conda:
        "../envs/bed_split.yaml"
    message:
        "{rule}: generate bed file containing only {wildcards.chr} from {input.bed}"
    shell:
        "(awk '{{if(/^{wildcards.chr}\t/) print($0)}}' {input.bed}  > {output}) &> {log}"
