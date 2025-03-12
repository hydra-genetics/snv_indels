__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepsomatic_t_only:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
        pon=config.get("reference", {}).get("pon", ""),
    output:
        tmpdir=directory("snv_indels/deepsomatic/{sample}_{type}_tmp"),
        vcf="snv_indels/deepsomatic/{sample}_{type}.vcf.gz"
    params:
        extra=config.get("deepsomatic", {}).get("extra", ""),
        model=config.get("deepsomatic", {}).get("model", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
    log:
        "snv_indels/deepsomatic/{sample}_{type}.deepsomatic.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic/{sample}_{type}.output.benchmark.tsv",
            config.get("deepsomatic", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic", {}).get("container", config["default_container"])
    message:
        "{rule}: Calling small variants from short read data in tumour only sample with DeepSomatic from {input.bam}"
    shell:
        """
        run_deepsomatic \
        --model_type={params.model} \
        --ref={input.ref} \
        --reads_tumor={input.bam} \
        --output_vcf={output.vcf} \
        --sample_name_tumor={params.name} \
        --num_shards={resources.threads} \
        --logging_dir={log} \
        --vcf_stats_report=true \
        --intermediate_results_dir {output.tmpdir} \
        --pon_filtering={input.pon} --process_somatic=true \
        --regions={input.bed} \
        {params.extra} \
        """

rule deepsomatic_tn:
    input:
        normal="alignment/samtools_merge_bam/{sample}_N.bam",
        tumor="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
        pon=config.get("reference", {}).get("pon", ""),
    output:
        tmpdir=directory("snv_indels/deepsomatic/{sample}_{type}_tmp"),
        vcf="snv_indels/deepsomatic/{sample}_{type}.vcf.gz"
    params:
        extra=config.get("deepsomatic", {}).get("extra", ""),
        model=config.get("deepsomatic", {}).get("model", ""),
        name_n=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
        name_t=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
    log:
        "snv_indels/deepsomatic/{sample}_{type}.deepsomatic.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic/{sample}_{type}.output.benchmark.tsv",
            config.get("deepsomatic", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic", {}).get("container", config["default_container"])
    message:
        "{rule}: Calling small variants from short read data in tumour only sample with DeepSomatic from {input.bam}"
    shell:
        """
        run_deepsomatic \
        --model_type={params.model} \
        --ref={input.ref} \
        --reads_normal={input.normal} \
        --reads_tumor={input.tumor} \
        --output_vcf={output.vcf} \
        --sample_name_normal={params.name_n} \
        --sample_name_tumor={params.name_t} \
        --num_shards={resources.threads} \
        --logging_dir={log} \
        --vcf_stats_report=true \
        --intermediate_results_dir {output.tmpdir} \
        --regions={input.bed} \
        {params.extra} \
        """