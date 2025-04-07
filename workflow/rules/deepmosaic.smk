__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepmosaic_deepmosaic_input:
    input:
        "...",
    output:
        "snv_indels/deepmosaic_deepmosaic_input/{sample}_{type}.output.txt",
    params:
        extra=config.get("deepmosaic_deepmosaic_input", {}).get("extra", ""),
    log:
        "snv_indels/deepmosaic_deepmosaic_input/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic_deepmosaic_input/{sample}_{type}.output.benchmark.tsv",
            config.get("deepmosaic_deepmosaic_input", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_deepmosaic_input", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_deepmosaic_input", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_deepmosaic_input", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_deepmosaic_input", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_deepmosaic_input", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_deepmosaic_input", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_deepmosaic_input", {}).get("container", config["default_container"])
    conda:
        "../envs/deepmosaic_deepmosaic_input.yaml"
    message:
        "{rule}: Do stuff on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}.input"
    wrapper:
        "..."


rule deepmosaic_deepmosaic_draw:
    input:
        "...",
    output:
        "snv_indels/deepmosaic_deepmosaic_draw/{sample}_{type}.output.txt",
    params:
        extra=config.get("deepmosaic_deepmosaic_draw", {}).get("extra", ""),
    log:
        "snv_indels/deepmosaic_deepmosaic_draw/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic_deepmosaic_draw/{sample}_{type}.output.benchmark.tsv",
            config.get("deepmosaic_deepmosaic_draw", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_deepmosaic_draw", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_deepmosaic_draw", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_deepmosaic_draw", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_deepmosaic_draw", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_deepmosaic_draw", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_deepmosaic_draw", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_deepmosaic_draw", {}).get("container", config["default_container"])
    conda:
        "../envs/deepmosaic_deepmosaic_draw.yaml"
    message:
        "{rule}: Do stuff on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}.input"
    wrapper:
        "..."


rule deepmosaic_deepmosaic_predict:
    input:
        "...",
    output:
        "snv_indels/deepmosaic_deepmosaic_predict/{sample}_{type}.output.txt",
    params:
        extra=config.get("deepmosaic_deepmosaic_predict", {}).get("extra", ""),
    log:
        "snv_indels/deepmosaic_deepmosaic_predict/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepmosaic_deepmosaic_predict/{sample}_{type}.output.benchmark.tsv",
            config.get("deepmosaic_deepmosaic_predict", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepmosaic_deepmosaic_predict", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepmosaic_deepmosaic_predict", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepmosaic_deepmosaic_predict", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepmosaic_deepmosaic_predict", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepmosaic_deepmosaic_predict", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepmosaic_deepmosaic_predict", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepmosaic_deepmosaic_predict", {}).get("container", config["default_container"])
    conda:
        "../envs/deepmosaic_deepmosaic_predict.yaml"
    message:
        "{rule}: Do stuff on snv_indels/{rule}/{wildcards.sample}_{wildcards.type}.input"
    wrapper:
        "..."
