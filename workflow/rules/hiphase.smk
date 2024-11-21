__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule hiphase:
    input:
        bam="alignment/pbmm2_align/{sample}_{type}.bam",
        ref=config.get("reference", {}).get("fasta", ""),
        snv_vcf="snv_indels/deepvariant/{sample}_{type}.merged.vcf.gz",
        sv_vcf="cnv_sv/pbsv/{sample}_{type}.vcf.gz" if config.get("hiphase", {}).get("sv_caller", False) else [],
        str_vcf="cnv_sv/trgt/{sample}_{type}.vcf.gz" if config.get("hiphase", {}).get("str_caller", False) else [],
    output:
        bam=temp("snv_indels/hiphase/{sample}_{type}.haplotagged.bam"),
        snv_vcf=temp("snv_indels/hiphase/{sample}_{type}.deepvariant.phased.vcf.gz"),
        sv_vcf=temp("snv_indels/hiphase/{sample}_{type}.pbsv.phased.vcf.gz")
        if config.get("hiphase", {}).get("sv_caller", False)
        else [],
        str_vcf=temp("snv_indels/hiphase/{sample}_{type}.trgt.phased.vcf.gz")
        if config.get("hiphase", {}).get("str_caller", False)
        else [],
    params:
        extra=config.get("hiphase", {}).get("extra", ""),
        in_sv_vcf=lambda wildcards: f"--vcf cnv_sv/pbsv/{wildcards.sample}_{wildcards.type}.vcf.gz"
        if config.get("hiphase", {}).get("sv_caller", False)
        else "",
        out_sv_vcf=lambda wildcards: f"--output-vcf snv_indels/hiphase/{wildcards.sample}_{wildcards.type}.pbsv.phased.vcf.gz"
        if config.get("hiphase", {}).get("sv_caller", False)
        else "",
        in_str_vcf=lambda wildcards: f"--vcf cnv_sv/trgt/{wildcards.sample}_{wildcards.type}.vcf.gz"
        if config.get("hiphase", {}).get("str_caller", False)
        else "",
        out_str_vcf=lambda wildcards: f"--output-vcf snv_indels/hiphase/{wildcards.sample}_{wildcards.type}.trgt.phased.vcf.gz"
        if config.get("hiphase", {}).get("str_caller", False)
        else "",
    log:
        "snv_indels/hiphase/{sample}_{type}.haplotagged.bam.log",
    benchmark:
        repeat("snv_indels/hiphase/{sample}_{type}.output.benchmark.tsv", config.get("hiphase", {}).get("benchmark_repeats", 1))
    threads: config.get("hiphase", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("hiphase", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("hiphase", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("hiphase", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("hiphase", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("hiphase", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("hiphase", {}).get("container", config["default_container"])
    message:
        "{rule}: phase variants and haplotag bam files using hiphase"
    shell:
        "hiphase "
        "--reference {input.ref} "
        "--vcf {input.snv_vcf} "
        "--output-vcf {output.snv_vcf} "
        "{params.in_sv_vcf} "
        "{params.out_sv_vcf} "
        "{params.in_str_vcf} "
        "{params.out_str_vcf} "
        "--bam {input.bam} "
        "--output-bam {output.bam} 2> {log}"
