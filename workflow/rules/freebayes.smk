__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule freebayes:
    input:
        ref=config["reference"]["fasta"],
        samples="alignment/mark_duplicates/{sample}_{type}_{chr}.bam",
        indexes="alignment/mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        regions="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        vcf=temp("snv_indels/freebayes/{sample}_{type}_{chr}.vcf"),
    params:
        extra="--target snv_indels/bed_split/design_bedfile_{chr}.bed %s"
        % config.get("freebayes", {}).get("extra", "--min-alternate-fraction 0.01 --genotype-qualities --strict-vcf"),
    log:
        "snv_indels/freebayes/{sample}_{type}_{chr}.unfilt.vcf.log",
    benchmark:
        repeat(
            "snv_indels/freebayes/{sample}_{type}_{chr}.benchmark.tsv", config.get("freebayes", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("freebayes", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("freebayes", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("freebayes", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("freebayes", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("freebayes", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("freebayes", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("freebayes", {}).get("container", config["default_container"])
    conda:
        "../envs/freebayes.yaml"
    message:
        "{rule}: call variants in {input.samples}"
    wrapper:
        "v1.3.1/bio/freebayes"
