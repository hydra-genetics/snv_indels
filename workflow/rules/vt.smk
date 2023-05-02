__author__ = "Jonas Almlöf, Patrik Smeds"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule vt_decompose:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.fix_af.vcf.gz",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz"),
    log:
        "snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz.benchmark.tsv",
            config.get("vt_decompose", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vt_decompose", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vt_decompose", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vt_decompose", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vt_decompose", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vt_decompose", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vt_decompose", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vt_decompose", {}).get("container", config["default_container"])
    message:
        "{rule}: decompose {input.vcf}"
    shell:
        "(vt decompose -s {input.vcf} | vt decompose_blocksub -o {output.vcf} -) &> {log}"


rule vt_normalize:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}.decomposed.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz"),
    log:
        "snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{sample}_{type}.normalized.vcf.gz.benchmark.tsv",
            config.get("vt_normalize", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vt_normalize", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vt_normalize", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vt_normalize", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vt_normalize", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vt_normalize", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vt_normalize", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vt_normalize", {}).get("container", config["default_container"])
    message:
        "{rule}: normalize {input.vcf}"
    shell:
        "(vt normalize -n -r {input.ref} -o {output.vcf} {input.vcf} ) &> {log}"
