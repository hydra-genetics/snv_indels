__author__ = "Monika Brandt"
__copyright__ = "Copyright 2024, Monika Brandt"
__email__ = "monika.brandt@scilifelab.uu.se"
__license__ = "GPL-3"


rule merge_af_complex_variants:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.normalized.sorted.vcf.gz",
        tabix="snv_indels/{caller}/{sample}_{type}.normalized.sorted.vcf.gz.tbi",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.normalized.merged_af.vcf.gz"),
    params:
        merge_method=config.get("merge_af_complex_variants", {}).get("merge_method", "skip"),
    log:
        "snv_indels/{caller}/{sample}_{type}.normalized.merged_af.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.normalized.merged_af.vcf.gz.benchmark.tsv",
            config.get("merge_af_complex_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("merge_af_complex_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("merge_af_complex_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("merge_af_complex_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("merge_af_complex_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("merge_af_complex_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("merge_af_complex_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("merge_af_complex_variants", {}).get("container", config["default_container"])
    message:
        "{rule}: For decomposed complex variants in {input.vcf} create one record with updated allele frequency"
    script:
        "../scripts/merge_af.py"
