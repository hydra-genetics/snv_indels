__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule glnexus:
    input:
        gvcfs=expand("snv_indels/deeptrio/{{sample}}_{{type}}/{trio_member}.g.vcf", trio_member=["child", "parent1", "parent2"]),
    output:
        bcf=temp("snv_indels/glnexus/{sample}_{type}.bcf"),
        dir=temp(directory("snv_indels/glnexus/{sample}_{type}/GLnexus.DB")),
    params:
        extra=config.get("glnexus", {}).get("extra", ""),
        glnexus_config=config.get("glnexus", {}).get("configfile", ""),
        in_gvcf=lambda wildcards, input: get_glnexus_input(wildcards, input),
    log:
        "snv_indels/glnexus/{sample}_{type}.bcf.log",
    benchmark:
        repeat(
            "snv_indels/glnexus/{sample}_{type}.bcf.benchmark.tsv",
            config.get("glnexus", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("glnexus", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("glnexus", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("glnexus", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("glnexus", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("glnexus", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("glnexus", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("glnexus", {}).get("container", config["default_container"])
    message:
        "{rule}: Run GLNexus for joint genotyping of gVCFs"
    shell:
        "(glnexus_cli "
        "--dir {output.dir} {params.extra} "
        "--config {params.glnexus_config} "
        "{params.in_gvcf} > {output.bcf}) &> {log}"
