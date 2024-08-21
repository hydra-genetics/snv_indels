__copyright__ = "Copyright 2023, Uppsala Universitet"
__email__ = "magdalena.z@scilifelab.uu.se"
__license__ = "GPL-3"

# compile_paraphase_file_list(wildcards)
GENE = ["smn1", "CR1", "AMY1A", "CTAG1A", "BOLA2"]
# SAMPLE = get_samples(samples) # Example names
# TYPES = get_unit_types(units, samples) # Example types


rule paraphase:
    input:
        bam="alignment/minimap2/{sample}_{type}.bam",
        bai="alignment/minimap2/{sample}_{type}.bam.bai",
        fasta=config.get("paraphase", {}).get("fasta", ""),
        faidx=config.get("paraphase", {}).get("fai", ""),
    output:
        merged_vcf="snv_indels/paraphase/{sample}_{type}.paraphase.vcf.gz",
        bam="snv_indels/paraphase/{sample}_{type}_realigned.paraphase.bam",
        #bai="snv_indels/paraphase/{sample}_{type}_realigned.paraphase.bam.bai",
        vcf_header="snv_indels/paraphase/{sample}_{type}.vcf_chromosome_header.vcf",
    params:
        genome=config.get("paraphase", {}).get("genome", ""),
        extra=config.get("paraphase", {}).get("extra", ""),
    log:
        "snv_indels/paraphase/{sample}_{type}.paraphase.log",
    benchmark:
        repeat(
            "snv_indels/paraphase/{sample}_{type}.out.benchmark.tsv",
            config.get("paraphase", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("paraphase", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("paraphase", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("paraphase", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("paraphase", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("paraphase", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("paraphase", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("paraphase", {}).get("container", config["default_container"])
    message:
        "{rule}: Calls SNVs on {input.bam} with paraphase to resolve SNVs in gene families"
    wrapper:
        "my-new-snakemake-wrapper-paraphrase/bio/paraphase"
        
        
# "file:///beegfs-storage/projects/wp3/nobackup/Workspace/magz_testing/snakemake-wrappers/bio/paraphase"
# "ppp/bio/paraphase"
# "my-new-snakemake-wrapper-paraphrase/bio/paraphase"
