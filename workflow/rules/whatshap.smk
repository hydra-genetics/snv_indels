__author__ = "Magdalena Z"
__copyright__ = "Copyright 2024, Uppsala Universitet"
__email__ = "magdalena.z@scilifelab.uu.se"
__license__ = "GPL-3"


rule whatshap_haplotag:
    input:
        vcf="snv_indels/whatshap/{sample}_{type}.whatshap_phased.vcf.gz",
        tbi="snv_indels/whatshap/{sample}_{type}.whatshap_phased.vcf.gz.tbi",
        aln="alignment/picard_mark_duplicates/{sample}_{type}_chr1.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_chr1.bam.bai",
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fai"],
    output:
        output="snv_indels/whatshap/{sample}_{type}.whatshap_haplotagged.bam",
    params:
        extra=config.get("whatshap_haplotag", {}).get("extra", ""),  # optionally use --ignore-linked-read, --tag-supplementary, etc.
    log:
        "snv_indels/whatshap/{sample}_{type}.whatshap_haplotagged.log",
    benchmark:
        "snv_indels/whatshap/{sample}_{type}.whatshap_haplotagged.tsv"
    threads: config.get("whatshap_haplotag", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("whatshap_haplotag", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("whatshap_haplotag", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("whatshap_haplotag", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("whatshap_haplotag", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("whatshap_haplotag", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("whatshap_phase", {}).get("container", config["default_container"])
    message:
        "{rule}: Phases reads from {input.aln} using whatshap phase results {input.vcf}"
    wrapper:
        "v3.5.2/bio/whatshap/haplotag"


rule whatshap_phase:
    input:
        reference=config["reference"]["fasta"],
        vcf="snv_indels/deepvariant/{sample}_{type}_chr1.vcf",
        tbi="snv_indels/deepvariant/{sample}_{type}_chr1.vcf.idx",
        aln="alignment/picard_mark_duplicates/{sample}_{type}_chr1.bam",
        bai="alignment/picard_mark_duplicates/{sample}_{type}_chr1.bam.bai",
    output:
        out="snv_indels/whatshap/{sample}_{type}.whatshap_phased.vcf.gz",
        #outindex="snv_indels/whatshap/{sample}_{type}.whatshap_phased.vcf.gz.tbi",
    log:
        "snv_indels/whatshap/{sample}_{type}.whatshap_phased.log",
    benchmark:
        "snv_indels/whatshap/{sample}_{type}.whatshap_phased.tsv"
    params:
        extra=config.get("whatshap_phase", {}).get("extra", ""),
    threads: config.get("whatshap_phase", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("whatshap_phase", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("whatshap_phase", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("whatshap_phase", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("whatshap_phase", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("whatshap_phase", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("whatshap_phase", {}).get("container", config["default_container"])
    message:
        "{rule}: Phases SNVs from {input.vcf} on {input.aln} with whatshap phase to resolve SNVs into haplotypes"
    shell:
        """
        (
        whatshap phase {params.extra} \
            --output {output.out} \
            --reference {input.reference} \
            {input.vcf} \
            {input.aln} \
        ) > {log} 2>&1 
        """


rule whatshap_index:
    input:
        "snv_indels/whatshap/{sample}_{type}.whatshap_haplotagged.bam",
    output:
        "snv_indels/whatshap/{sample}_{type}.whatshap_haplotagged.bam.bai",
    params:
        extra=config.get("whatshap_haplotag", {}).get("extra", ""),  # optionally use --ignore-linked-read, --tag-supplementary, etc.
    log:
        "snv_indels/whatshap/{sample}_{type}.whatshap_index.log",
    threads: config.get("whatshap_haplotag", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("whatshap_haplotag", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("whatshap_haplotag", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("whatshap_haplotag", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("whatshap_haplotag", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("whatshap_haplotag", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("whatshap_index", {}).get("container", config["default_container"])
    wrapper:
        "v3.7.0/bio/samtools/index"

