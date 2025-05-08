__author__ = "Jessika Nordin, Camille Clouard"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepsomatic_t_only:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        tmpdir=temp(directory("snv_indels/deepsomatic_t_only/{sample}_{type}.tmp")),
        vcf=temp("snv_indels/deepsomatic_t_only/{sample}_{type}.vcf.gz"),
    params:
        extra=config.get("deepsomatic_t_only", {}).get("extra", ""),
        model=config.get("deepsomatic_t_only", {}).get("model", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
        pon=config.get("deepsomatic_t_only", {}).get("pon", ""),
    log:
        "snv_indels/deepsomatic_t_only/{sample}_{type}.deepsomatic.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic_t_only/{sample}_{type}.output.benchmark.tsv",
            config.get("deepsomatic_t_only", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepsomatic_t_only", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic_t_only", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic_t_only", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic_t_only", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic_t_only", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic_t_only", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic_t_only", {}).get("container", config["default_container"])
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
        --intermediate_results_dir \
        {output.tmpdir} \
        --regions={input.bed} \
        {params.pon} \
        {params.extra} 
        """


rule deepsomatic_tn:
    input:
        normal="alignment/samtools_merge_bam/{sample}_N.bam",
        tumor="alignment/samtools_merge_bam/{sample}_T.bam",
        bai_n="alignment/samtools_merge_bam/{sample}_N.bam.bai",
        bai_t="alignment/samtools_merge_bam/{sample}_T.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        tmpdir=temp(directory("snv_indels/deepsomatic_tn/{sample}.tmp")),
        vcf=temp("snv_indels/deepsomatic_tn/{sample}.vcf.gz"),
    params:
        extra=config.get("deepsomatic_tn", {}).get("extra", ""),
        model=config.get("deepsomatic_tn", {}).get("model", ""),
        name_n=lambda wildcards: f"{wildcards.sample}_N",
        name_t=lambda wildcards: f"{wildcards.sample}_T",
    log:
        "snv_indels/deepsomatic_tn/{sample}.deepsomatic.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic_tn/{sample}.output.benchmark.tsv",
            config.get("deepsomatic_tn", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepsomatic_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic_tn", {}).get("container", config["default_container"])
    message:
        "{rule}: Calling small variants from short read data in tumour/normal samples with DeepSomatic from {input.tumor}"
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
        {params.extra}
        """
