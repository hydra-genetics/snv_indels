__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule deeptrio_make_examples:
    input:
        bams=lambda wildcards: get_trio_bams(wildcards),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        examples=temp(
            expand(
                "snv_indels/deeptrio/{{trioid}}/make_examples_{trio_member}.tfrecord-{{shard}}-of-{nshards:05}.gz",
                nshards=config.get("deeptrio_make_examples").get("n_shards", 2),
                trio_member=["child", "parent1", "parent2"],
            )
        ),
        gvcf_tfrecords=temp(
            expand(
                "snv_indels/deeptrio/{{trioid}}/gvcf_{trio_member}.tfrecord-{{shard}}-of-{nshards:05}.gz",
                nshards=config.get("deeptrio_make_examples").get("n_shards", 2),
                trio_member=["child", "parent1", "parent2"],
            )
        ),
    params:
        examples=lambda wildcards, output: get_make_examples_tfrecord(
            wildcards, output, config.get("deeptrio_make_examples").get("n_shards", 2)
        ),
        extra=config.get("deeptrio_make_examples", {}).get("extra", ""),
        shard=lambda wildcards: int(wildcards.shard),
        nshards=config.get("deeptrio_make_examples").get("n_shards", 2),
    log:
        "snv_indels/deeptrio/{trioid}/make_examples_{shard}.output.log",
    benchmark:
        repeat(
            "snv_indels/deeptrio/{trioid}/make_examples_{shard}.output.benchmark.tsv",
            config.get("deeptrio_make_examples", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deeptrio_make_examples", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deeptrio_make_examples", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deeptrio_make_examples", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deeptrio_make_examples", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deeptrio_make_examples", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deeptrio_make_examples", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deeptrio_make_examples", {}).get("container", config["default_container"])
    conda:
        "../envs/deeptrio.yaml"
    message:
        "{rule}: Run deeptrio make_examples on {input.bams} "
    shell:
        "(make_examples "
        "--mode 'calling' "
        "--ref {input.ref} "
        "--reads {input.bams[0]} "
        "--reads_parent1 {input.bams[1]}  "
        "--reads_parent2 {input.bams[2]} "
        "--examples {params.examples} "
        "--gvcf snv_indels/deeptrio/{wildcards.trioid}/gvcf.tfrecord@{params.nshards}.gz "
        "{params.extra} --task {params.shard}) &> {log}"


rule deeptrio_call_variants:
    input:
        examples=expand(
            "snv_indels/deeptrio/{{trioid}}/make_examples_{{trio_member}}.tfrecord-{shard}-of-{nshards:05}.gz",
            shard=[f"{x:05}" for x in range(config.get("deeptrio_make_examples").get("n_shards", 2))],
            nshards=config.get("deeptrio_make_examples").get("n_shards", 2),
        ),
    output:
        outfile=temp("snv_indels/deeptrio/{trioid}/call_variants_output_{trio_member}.tfrecord.gz"),
    params:
        cuda="CUDA_VISIBLE_DEVICES={}".format(os.getenv("CUDA_VISIBLE_DEVICES"))
        if os.getenv("CUDA_VISIBLE_DEVICES") is not None
        else "",
        examples=lambda wildcards, output: get_make_examples_tfrecord(
            wildcards, output, config.get("deeptrio_make_examples").get("n_shards", 2), program="deeptrio"
        ),
        extra=config.get("deeptrio_call_variants", {}).get("extra", ""),
        model=lambda wildcards: get_deeptrio_model(wildcards),
    log:
        "snv_indels/deeptrio/{trioid}/call_variants_{trio_member}.output.log",
    benchmark:
        repeat(
            "snv_indels/deeptrio/{trioid}/call_variants_{trio_member}.output.benchmark.tsv",
            config.get("deeptrio_call_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deeptrio_call_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        gres=config.get("deeptrio_call_variants", {}).get("gres", ""),
        mem_mb=config.get("deeptrio_call_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deeptrio_call_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deeptrio_call_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deeptrio_call_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deeptrio_call_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deeptrio_call_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/deeptrio.yaml"
    message:
        "{rule}: Run deeptrio call_variants on {params.examples}"
    shell:
        "({params.cuda} call_variants "
        "--checkpoint {params.model} "
        "--outfile {output.outfile} "
        "--examples {params.examples} "
        "{params.extra}) &> {log}"


rule deeptrio_postprocess_variants:
    input:
        call_variants_record="snv_indels/deeptrio/{trioid}/call_variants_output_{trio_member}.tfrecord.gz",
        gvcf_records=expand(
            "snv_indels/deeptrio/{{trioid}}/gvcf_{{trio_member}}.tfrecord-{shard}-of-{nshards:05}.gz",
            shard=[f"{x:05}" for x in range(config.get("deeptrio_make_examples").get("n_shards", 2))],
            nshards=config.get("deeptrio_make_examples").get("n_shards", 2),
        ),
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/deeptrio/{trioid}_{trio_member}.vcf"),
        gvcf=temp("snv_indels/deeptrio/{trioid}_{trio_member}.g.vcf"),
    params:
        extra=lambda wildcards, input, output: deeptrio_postprocess_variants_args(
            wildcards,
            input,
            output,
            "deeptrio_make_examples",
            config.get("deeptrio_postprocess_variants", {}).get("extra", ""),
        ),
    log:
        "snv_indels/deeptrio/{trioid}/postprocess_variants_{trio_member}.output.log",
    benchmark:
        repeat(
            "snv_indels/deeptrio/{trioid}/postprocess_variants_{trio_member}.output.benchmark.tsv",
            config.get("deeptrio_postprocess_variants", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deeptrio_postprocess_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deeptrio_postprocess_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deeptrio_postprocess_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deeptrio_postprocess_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deeptrio_postprocess_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deeptrio_postprocess_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deeptrio_postprocess_variants", {}).get("container", config["default_container"])
    conda:
        "../envs/deeptrio.yaml"
    message:
        "{rule}: Run deeptrio postprocess_variants on {input.call_variants_record}"
    shell:
        "(postprocess_variants "
        "--infile {input.call_variants_record} "
        "--ref {input.ref} "
        "--outfile {output.vcf} "
        "--gvcf_outfile {output.gvcf} "
        "--novcf_stats_report "
        "{params.extra} ) &> {log}"
