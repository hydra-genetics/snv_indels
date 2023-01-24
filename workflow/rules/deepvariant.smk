__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepvariant_make_examples:
    input:
        bam="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        examples=temp("snv_indels/deepvariant/{sample}_{type}_{chr}/make_examples.tfrecord@{threads}.gz"),
    params:
        extra=lambda wildcards, output: get_make_example_args(wildcards, output, "deepvariant_make_examples"),
        seq_num=config.get("deepvariant_make_examples", {}).get("threads", config["default_resources"]["threads"]) - 1,
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}/{threads}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}/{threads}.output.benchmark.tsv",
            config.get("deepvariant_make_examples", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_make_examples", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_make_examples", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_make_examples", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_make_examples", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_make_examples", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_make_examples", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_make_examples", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Deepvariant make_examples on {wildcards.sample}_{wildcards.type}"
    shell:
        """
        seq 0 {params.seq_num} | \
        parallel -q --halt 2 --line-buffer \
        make_examples \
        --mode 'calling' \
        --ref {input.ref} \
        --reads {input.bam} \
        --examples {output.examples} \
        {params.extra} --task {{}} &> {log}
        """


use rule deepvariant_make_examples as deepvariant_make_examples_gvcf with:
    output:
        examples=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/make_examples.tfrecord@{threads}.gz"),
        gvcf=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/gvcf.tfrecord@{threads}.gz")


rule deepvariant_call_variants:
    input:
        examples=lambda wildcards: get_examples_infile(wildcards, "deepvariant_make_examples", "vcf"),
        model_cpkt=config.get("deepvariant_call_variants", {}).get("model_file", ""),
    output:
        outfile="snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
    params:
        extra=config.get("deepvariant_call_variants", {}).get("extra", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant_call_variants/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("deepvariant_call_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_call_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_call_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_call_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_call_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_call_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_call_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_call_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Deepvariant call_variants on {wildcards.sample}_{wildcards.type}"
    shell:
        "(call_variants "
        "--checkpoint {input.model_cpkt} "
        "--outfile {output.outfile} "
        "--examples {input.examples} "
        "{params.extra}) &> {log}"


rule deepvariant_call_variants_gvcf:
    input:
        examples=lambda wildcards: get_examples_infile(wildcards, "deepvariant_make_examples", "gvcf"),
        model_cpkt=config.get("deepvariant_call_variants", {}).get("model_file", ""),
    output:
        outfile="snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",


rule deepvariant_postprocess_variants:
    input:
        infile="snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.vcf"),
    params:
        extra=lambda wildcards, input, output: get_postprocess_variants_args(wildcards, input, output, "deepvariant_make_examples"),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant_postprocess_variants/{sample}_{type}_{chr}.output.benchmark.tsv",
            config.get("deepvariant_postprocess_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_postprocess_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_postprocess_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_postprocess_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_postprocess_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_postprocess_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_postprocess_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_postprocess_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Deepvariant postprocess_variants on {wildcards.sample}_{wildcards.type}"
    shell:
        "(posprocess_variants "
        "--infile {input.infile} "
        "--outfile {output.vcf} "
        "{params.extra}) &> {log}"


use rule deepvariant_postprocess_variants as deepvariant_postprocess_variants_gvcf with:
    input:
        infile="snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
        gvcf_tfrecord=lambda wildcards: get_gvcf_tfrecord(wildcards, "deepvariant_make_examples"),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}.vcf"),
        gvcf=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}.g.vcf")