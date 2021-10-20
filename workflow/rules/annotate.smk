# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule annotate:
    input:
        vcf="snv_indels/{caller}/{file}.vcf.gz",
        tabix="snv_indels/{caller}/{file}.vcf.gz.tbi",
        cache=config["annotate"]["input"]["vep_cache"],
        fasta=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/{caller}/{file}.annotated.vcf"),
    params:
        extra=config["annotate"]["params"]["extra"],
        mode=config.get("annotate", {}).get("params", {}).get("mode", "--offline --cache"),
    log:
        "snv_indels/{caller}/{file}.annotated.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{file}.annotated.vcf.gz.benchmark.tsv",
            config.get("annotate", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("annotate", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("annotate", {}).get("container", config["default_container"])
    conda:
        "../envs/annotate.yaml"
    message:
        "{rule}: Sort vcf snv_indels/{wildcards.caller}/{wildcards.file}.vcf.gz"
    shell:
        "(vep --vcf --no_stats -o {output.vcf} -i {input.vcf} --dir_cache {input.cache} --fork {threads} --refseq {params.mode} --fasta {input.fasta} {params.extra} ) &> {log}"
