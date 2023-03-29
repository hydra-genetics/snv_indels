__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_concat:
    input:
        calls=expand(
            "{{file}}_{chr}.{{vcf}}.gz",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
        index=expand(
            "{{file}}_{chr}.{{vcf}}.gz.tbi",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
    output:
        vcf=temp("{file}.merged.{vcf}"),
    params:
        extra=config.get("bcftools_concat", {}).get("extra", ""),
    log:
        "{file}.merged.{vcf}.log",
    benchmark:
        repeat(
            "{file}.merged.{vcf}.benchmark.tsv",
            config.get("bcftools_concat", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_concat", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_concat", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_concat", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_concat", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_concat", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_concat", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_concat", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: concatenate {input.calls}"
    wrapper:
        "v1.25.0/bio/bcftools/concat"


rule bcftools_sort:
    input:
        vcf="{file}.vcf.gz",
        tbi="{file}.vcf.gz.tbi",
    output:
        vcf=temp("{file}.sorted.vcf.gz"),
    log:
        "{file}.sorted.vcf.gz.log",
    benchmark:
        repeat(
            "{file}.sorted.vcf.gz.benchmark.tsv",
            config.get("bcftools_sort", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_sort", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_sort", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_sort", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_sort", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_sort", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_sort", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_sort", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: sort {input.vcf}"
    wrapper:
        "v1.25.0/bio/bcftools/sort"


rule bcftools_view:
    input:
        bcf="{file}.bcf",
    output:
        vcf=temp("{file}.vcf.gz"),
    log:
        "{file}.vcf.gz.log",
    params:
        extra=config.get("bcftools_concat", {}).get("extra", ""),
    benchmark:
        repeat(
            "{file}.vcf.gz.benchmark.tsv",
            config.get("bcftools_view", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_view", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_view", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_view", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: convert {input.bcf} to {output.vcf}"
    wrapper:
        "v1.25.0/bio/bcftools/view"
