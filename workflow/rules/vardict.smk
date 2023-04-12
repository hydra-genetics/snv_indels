__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule vardict:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        reference=config["reference"]["fasta"],
        regions="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        vcf=temp("snv_indels/vardict/{sample}_{type}_{chr}.vcf"),
    params:
        extra=config.get("vardict", {}).get("extra", "-Q 1"),
        bed_columns=config.get("vardict", {}).get("bed_columns", "-c 1 -S 2 -E 3 -g 4"),
        allele_frequency_threshold=config.get("vardict", {}).get("allele_frequency_threshold", "0.01"),
        sample_name="{sample}_{type}",
    log:
        "snv_indels/vardict/{sample}_{type}_{chr}.vcf.log",
    benchmark:
        repeat("snv_indels/vardict/{sample}_{type}_{chr}.benchmark.tsv", config.get("vardict", {}).get("benchmark_repeats", 1))
    threads: config.get("vardict", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vardict", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vardict", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vardict", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vardict", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vardict", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vardict", {}).get("container", config["default_container"])
    message:
        "{rule}: call variants in {input.bam}"
    wrapper:
        "v1.3.1/bio/vardict"
