__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcbio_variation_recall_ensemble:
    input:
        vcfs=expand(
            "snv_indels/{caller}/{{sample}}_{{type}}.normalized.sorted.vcf.gz",
            caller=config.get("bcbio_variation_recall_ensemble", {}).get("callers", []),
        ),
        tabix=expand(
            "snv_indels/{caller}/{{sample}}_{{type}}.normalized.sorted.vcf.gz.tbi",
            caller=config.get("bcbio_variation_recall_ensemble", {}).get("callers", []),
        ),
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf.gz"),
    params:
        support=config.get("bcbio_variation_recall_ensemble", {}).get("support", "1"),
        sort_order=get_bvre_params_sort_order,
    log:
        "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf.gz.benchmark.tsv",
            config.get("bcbio_variation_recall_ensemble", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcbio_variation_recall_ensemble", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("bcbio_variation_recall_ensemble", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcbio_variation_recall_ensemble", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("bcbio_variation_recall_ensemble", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcbio_variation_recall_ensemble", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcbio_variation_recall_ensemble", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("bcbio_variation_recall_ensemble", {}).get("container", config["default_container"])
    conda:
        "../envs/bcbio_variation_recall.yaml"
    message:
        "{rule}: combine vcfs from different callers into {output.vcf} {params.sort_order}"
    shell:
        "(bcbio-variation-recall ensemble -n {params.support} --names {params.sort_order} {output.vcf} {input.ref} {input.vcfs}) &> {log}"
