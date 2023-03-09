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
        examples=temp(f"snv_indels/deepvariant/{{sample}}_{{type}}_{{chr}}/make_examples.tfrecord-{{shard}}-of-{config.get('deepvariant_make_examples').get('n_shards', 10):05}.gz"),
    params:
        examples=lambda wildcards, output: get_make_examples_tfrecord(wildcards, output, 
        config.get("deepvariant_make_examples").get('n_shards', 10)),
        extra=lambda wildcards, output: deepvariant_make_example_args(wildcards, output),
        shard = lambda wildcards: int(wildcards.shard),
    log:
        "snv_indels/deepvariant/{sample}_{type}_{chr}/make_examples_{shard}.output.log",
    benchmark:
        repeat(
            "snv_indels/deepvariant/{sample}_{type}_{chr}/make_examples_{shard}.output.benchmark.tsv",
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
        "{rule}: Run deepvariant make_examples on {input.bam}"
    shell:
        "(time make_examples "
        "--mode 'calling' "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--examples {params.examples} "
        "{params.extra} --task {params.shard}) &> {log}"


rule deepvariant_call_variants:
    input:
        examples=lambda wildcards: get_example_records(wildcards, config.get("deepvariant_make_examples").get('n_shards', 10)),
    output:
        outfile=temp("snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz"),
    params:
        examples=lambda wildcards, input: get_make_examples_tfrecord(
            wildcards, input, config.get("deepvariant_make_examples").get('n_shards', 10)),
        extra=config.get("deepvariant_call_variants", {}).get("extra", ""),
        model=config.get("deepvariant_call_variants", {}).get("model", ""),
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
        "{rule}: Run deepvariant call_variants on {params.examples}"
    shell:
        "(time call_variants "
        "--checkpoint {params.model} "
        "--outfile {output.outfile} "
        "--examples {params.examples} "
        "{params.extra}) &> {log}"


rule deepvariant_postprocess_variants:
    input:
        call_variants_record="snv_indels/deepvariant/{sample}_{type}_{chr}/call_variants_output.tfrecord.gz",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.vcf"),
        gvcf=temp("snv_indels/deepvariant/{sample}_{type}_{chr}.g.vcf")
        if config.get("deepvariant_postprocess_variants", {}).get("vcf_type", "vcf") == "gvcf"
        else [],
    params:
        extra=lambda wildcards, input, output: get_postprocess_variants_args(
            wildcards, input, output, "deepvariant_make_examples",
            config.get("deepvariant_postprocess_variants", {}).get("extra", ""),
        ),
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
        mem_per_cpu=config.get("deepvariant_postprocess_variants", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("deepvariant_postprocess_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_postprocess_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_postprocess_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_postprocess_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant.yaml"
    message:
        "{rule}: Run deepvariant postprocess_variants on {input.call_variants_record}"
    shell:
        "(time postprocess_variants "
        "--infile {input.call_variants_record} "
        "--ref {input.ref} "
        "--outfile {output.vcf} "
        "{params.extra}) &> {log}"
