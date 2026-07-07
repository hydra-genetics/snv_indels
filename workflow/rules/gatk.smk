__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"


rule gatk_mutect2:
    input:
        map="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        bam=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.bam"),
        bai=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.bai"),
        stats=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.stats"),
        vcf=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz"),
        tbi=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.tbi"),
        f1f2=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.f1r2.tar.gz"),
    params:
        extra=lambda wildcards: get_gatk_mutect2_extra(wildcards, "gatk_mutect2"),
    log:
        "snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.benchmark.tsv",
            config.get("gatk_mutect2", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_mutect2", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mutect2", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mutect2", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mutect2", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mutect2", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mutect2", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mutect2", {}).get("container", config["default_container"])
    message:
        "{rule}: call variants in {input.map}"
    wrapper:
        "v1.5.0/bio/gatk/mutect"


rule gatk_mutect2_gvcf:
    input:
        map="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_{chr}.bam.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
        bed="snv_indels/bed_split/design_bedfile_{chr}.bed",
    output:
        bam=temp("snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.bam"),
        bai=temp("snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.bai"),
        stats=temp("snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.stats"),
        vcf=temp("snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz"),
        tbi=temp("snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.tbi"),
    params:
        extra=lambda wildcards: get_gatk_mutect2_extra(wildcards, "gatk_mutect2_gvcf"),
    log:
        "snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.benchmark.tsv",
            config.get("gatk_mutect2_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mutect2_gvcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mutect2_gvcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mutect2_gvcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mutect2_gvcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mutect2_gvcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mutect2_gvcf", {}).get("container", config["default_container"])
    message:
        "{rule}: generate gvcf from {input.map}"
    wrapper:
        "v1.5.0/bio/gatk/mutect"


rule gatk_selectvariants_gvcf_to_vcf:
    """Extract variant-only sites from a Mutect2 gVCF, replacing a separate
    VCF-mode Mutect2 call. The gVCF stats and the locally reassembled bam are
    passed through unchanged so that downstream gatk_mutect2_merge_stats /
    gatk_mutect2_filter / alignment_samtools_merge_bam_mutect2-style rules
    continue to work without modification.

    Enable in a pipeline with:
        ruleorder: gatk_selectvariants_gvcf_to_vcf > gatk_mutect2
    """
    input:
        gvcf="snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz",
        tbi="snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.tbi",
        stats="snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.gz.stats",
        bam="snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.bam",
        bai="snv_indels/gatk_mutect2_gvcf/{sample}_{type}_{chr}.g.vcf.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz"),
        tbi=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.tbi"),
        stats=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.stats"),
        bam=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.bam"),
        bai=temp("snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.bai"),
    params:
        extra=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("extra", "--exclude-non-variants --remove-unused-alternates"),
    log:
        "snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.selectvariants.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2/{sample}_{type}_{chr}.unfiltered.vcf.gz.selectvariants.benchmark.tsv",
            config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_selectvariants_gvcf_to_vcf", {}).get("container", config["default_container"])
    message:
        "{rule}: extract variants from gVCF {input.gvcf}"
    shell:
        """
        gatk SelectVariants -R {input.fasta} -V {input.gvcf} -O {output.vcf} {params.extra} 2> {log}
        cp {input.stats} {output.stats}
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        """


rule gatk_mutect2_filter:
    input:
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.merged.unfiltered.vcf.gz",
        tbi="snv_indels/gatk_mutect2/{sample}_{type}.merged.unfiltered.vcf.gz.tbi",
        stats="snv_indels/gatk_mutect2/{sample}_{type}.unfiltered.vcf.gz.stats",
        ref=config.get("reference", {}).get("fasta", ""),
    output:
        vcf=temp("snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz"),
    params:
        extra=lambda wildcards, input: "%s --stats %s" % (config.get("gatk_mutect2_filter", {}).get("extra", ""), input.stats),
    log:
        "snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz.benchmark.tsv",
            config.get("gatk_mutect2_filter", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_mutect2_filter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mutect2_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mutect2_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mutect2_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mutect2_filter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mutect2_filter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mutect2_filter", {}).get("container", config["default_container"])
    message:
        "{rule}: softfilter mutect2 variants in {input.vcf} to {output.vcf}"
    wrapper:
        "v1.5.0/bio/gatk/filtermutectcalls"


rule gatk_mutect2_merge_stats:
    input:
        stats=expand(
            "snv_indels/gatk_mutect2/{{sample}}_{{type}}_{chr}.unfiltered.vcf.gz.stats",
            chr=extract_chr(
                "%s.fai" % (config["reference"]["fasta"]), filter_out=config.get("reference", {}).get("skip_chrs", [])
            ),
        ),
    output:
        stats=temp("snv_indels/gatk_mutect2/{sample}_{type}.unfiltered.vcf.gz.stats"),
    params:
        stats=lambda wildcards, input: " -stats ".join(input.stats),
    log:
        "snv_indels/gatk_mutect2/{sample}_{type}.unfiltered.vcf.gz.stats.log",
    benchmark:
        repeat(
            "snv_indels/gatk_mutect2/{sample}_{type}.unfiltered.vcf.gz.stats.benchmark.tsv",
            config.get("gatk_mutect2_merge_stats", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_mutect2_merge_stats", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("gatk_mutect2_merge_stats", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_mutect2_merge_stats", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_mutect2_merge_stats", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_mutect2_merge_stats", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("gatk_mutect2_merge_stats", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_mutect2_merge_stats", {}).get("container", config["default_container"])
    message:
        "{rule}: merge mutect2 stats files into {output.stats}"
    shell:
        "(gatk MergeMutectStats "
        "-O {output.stats} "
        "-stats {params.stats}) &> {log}"
