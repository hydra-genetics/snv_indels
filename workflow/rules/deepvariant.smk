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
        examples_dir=temp(directory("snv_indels/deepvariant/{sample}_{type}_{chr}/")),
    params:
        extra=lambda wildcards, output: get_make_example_args(wildcards, output, "deepvariant_make_examples", "vcf"),
        seq_num=config.get("deepvariant_make_examples", {}).get("threads", config["default_resources"]["threads"]) - 1,
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}/make_examples.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}/make_examples.output.benchmark.tsv",
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
        --examples {output.examples_dir}/make_examples.tfrecord@{threads}.gz \
        {params.extra} --task {{}} &> {log}
        """


use rule deepvariant_make_examples as deepvariant_make_examples_gvcf with:
    output:
        examples_dir=temp(directory("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/")),
    params:
        extra=lambda wildcards, output: get_make_example_args(wildcards, output, "deepvariant_make_examples_gvcf", "gvcf"),
        seq_num=config.get("deepvariant_make_examples", {}).get("threads", config["default_resources"]["threads"]) - 1,
    log:
        "snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/make_examples.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/make_examples.output.benchmark.tsv",
            config.get("deepvariant_make_examples_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_make_examples_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_make_examples_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_make_examples_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_make_examples_gvcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_make_examples_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_make_examples_gvcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_make_examples_gvcf", {}).get("container", config["default_container"])


rule deepvariant_call_variants:
    input:
        examples_dir="snv_indels/deepvariant/{sample}_{type}_{chr}/",
    output:
        outfile=temp("snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz"),
    params:
        examples=lambda wildcards, input: get_examples_infile(wildcards, input, "deepvariant_make_examples"),
        extra=config.get("deepvariant_call_variants", {}).get("extra", ""),
        model_dir=config.get("deepvariant_call_variants", {}).get("model_dir", ""),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants.output.benchmark.tsv",
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
        "--checkpoint {params.model_dir}/model.ckpt "
        "--outfile {output.outfile} "
        "--examples {params.examples} "
        "{params.extra}) &> {log}"


use rule deepvariant_call_variants as deepvariant_call_variants_gvcf with:
    input:
        examples_dir="snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/",
    output:
        outfile=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz"),
    params:
        examples=lambda wildcards, input: get_examples_infile(wildcards, input, "deepvariant_make_examples_gvcf"),
        extra=config.get("deepvariant_call_variants_gvcf", {}).get("extra", ""),
        model_dir=config.get("deepvariant_call_variants_gvcf", {}).get("model_dir", ""),
    log:
        "snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants.output.benchmark.tsv",
            config.get("deepvariant_call_variants_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_call_variants_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_call_variants_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_call_variants_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_call_variants_gvcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_call_variants_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_call_variants_gvcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_call_variants_gvcf", {}).get("container", config["default_container"])


rule deepvariant_postprocess_variants:
    input:
        infile="snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.vcf"),
    params:
        extra=lambda wildcards, input, output: get_postprocess_variants_args(wildcards, input, output, "deepvariant_make_examples"),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}/postprocess_variants.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}/postprocess_variants.output.benchmark.tsv",
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
        "(postprocess_variants "
        "--infile {input.infile} "
        "--ref {input.ref} "
        "--outfile {output.vcf} "
        "{params.extra}) &> {log}"


use rule deepvariant_postprocess_variants as deepvariant_postprocess_variants_gvcf with:
    input:
        examples_dir="snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/",
        infile="snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
        ref=config.get("reference", {}).get("fasta", ""),
    params:
        extra=lambda wildcards, input, output: get_postprocess_variants_args(wildcards, input, output, "deepvariant_make_examples_gvcf"),
    output:
        vcf=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}.vcf"),
        gvcf=temp("snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}.g.vcf")
    log:
        "snv_indels/deepvariant_gvcf/{sample}_{type}_{chr}/postprocess_variants.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant_gcvf/{sample}_{type}_{chr}/postprocess_variants.output.benchmark.tsv",
            config.get("deepvariant_postprocess_variants_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepvariant_postprocess_variants_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepvariant_postprocess_variants_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_postprocess_variants_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_postprocess_variants_gvcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_postprocess_variants_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_postprocess_variants_gvcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_postprocess_variants_gvcf", {}).get("container", config["default_container"])
