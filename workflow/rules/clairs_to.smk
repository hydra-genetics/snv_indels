__author__ = "Camille Clouard"
__copyright__ = "Copyright 2025, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


rule clairs_to_call:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref=config.get("reference", {}).get("fasta", ""),
        bed=config.get("reference", {}).get("design_bed", ""),
    output:
        snv=temp("snv_indels/clairs_to/{sample}_{type}_snv.vcf.gz"),
        indel=temp("snv_indels/clairs_to/{sample}_{type}_indel.vcf.gz"),
    params:
        extra=config.get("clairs_to_call", {}).get("extra", ""),
        platform=config.get("clairs_to_call", {}).get("platform", ""),
        snv_min_af=config.get("clairs_to_call", {}).get("snv_min_af", 0.05),
        indel_min_af=config.get("clairs_to_call", {}).get("indel_min_af", 0.1),
        outdir=directory(lambda w, output: os.path.dirname(output[0])),
    log:
        "snv_indels/clairs_to/{sample}_{type}.varcall.log",
    benchmark:
        repeat(
            "snv_indels/clairs_to/{sample}_{type}.varcall.benchmark.tsv",
            config.get("clairs_to_call", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("clairs_to_call", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("clairs_to_call", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("clairs_to_call", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("clairs_to_call", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("clairs_to_call", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("clairs_to_call", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("clairs_to_call", {}).get("container", config["default_container"])
    message:
        "{rule}: Long-read somatic small variant calling in only tumor samples with ClairS-TO."
    shell:
        "run_clairs_to "
        "--tumor_bam_fn {input.bam} "
        "--ref_fn {input.ref} "
        "--threads {resources.threads} " 
        "--platform {params.platform} "
        "--output_dir {params.outdir} "
        "-s {wildcards.sample} "
        "--bed_fn {input.bed} "
        "--snv_min_af {params.snv_min_af} "
        "--indel_min_af {params.indel_min_af} "
        "--disable_verdict "
        "--snv_output_prefix {wildcards.sample}_{wildcards.type}_snv "
        "--indel_output_prefix {wildcards.sample}_{wildcards.type}_indel "
        "{params.extra} "
        "> {log}"


rule clairs_to_concat:
    input:
        snv="snv_indels/clairs_to/{sample}_{type}_snv.vcf.gz",
        indel="snv_indels/clairs_to/{sample}_{type}_indel.vcf.gz",
    output:
        vcf=temp("snv_indels/clairs_to/{sample}_{type}.snv-indels.vcf.gz"),
        tmp=temp("snv_indels/clairs_to/{sample}_{type}.snv-indels.unsorted.vcf.gz"),
    params:
        extra=config.get("clairs_to_concat", {}).get("extra", ""),
    log:
        "snv_indels/clairs_to/{sample}_{type}.concat.log",
    benchmark:
        repeat(
            "snv_indels/clairs_to/{sample}_{type}.concat.benchmark.tsv",
            config.get("clairs_to_concat", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("clairs_to_concat", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("clairs_to_concat", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("clairs_to_concat", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("clairs_to_concat", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("clairs_to_concat", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("clairs_to_concat", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("clairs_to_concat", {}).get("container", config["default_container"])
    message:
        "{rule}: Concatenate the output of ClairS-TO into a single VCF."
    shell:
        """
        bcftools concat {params.extra} -a -Oz -o {output.tmp} {input.snv} {input.indel} 2> {log}
        bcftools sort -Oz -o {output.vcf} {output.tmp}
        """
