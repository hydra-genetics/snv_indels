__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosaicforecast_input:
    input:
        vcf="snv_indels/deepsomatic_t_only/{sample}_{type}.vcf.gz",
    output:
        variants=temp("snv_indels/mosaicforecast/{sample}_{type}.input"),
    params:
        extra=config.get("mosaicforecast_input", {}).get("extra", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
    log:
        "snv_indels/mosaicforecast/{sample}_{type}.input.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast/{sample}_{type}.input_benchmark.tsv",
            config.get("mosaicforecast_input", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosaicforecast_input", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosaicforecast_input", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosaicforecast_input", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosaicforecast_input", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosaicforecast_input", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosaicforecast_input", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosaicforecast_input", {}).get("container", config["default_container"])
    message:
        "{rule}: make input file for mosaic forecast by making a file with candidate variants based in {input.vcf}"
    shell:
        """(bcftools view \
        {input.vcf} \
        -f PASS \
        --exclude 'FORMAT/DP<8' \
        {params.extra} \
        | bcftools query \
        -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t{params.name}\n' \
        > {output.variants} ) &> {log}"""


rule mosaicforecast:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",        
        fasta=config.get("reference", {}).get("fasta", ""),
        variants="snv_indels/mosaicforecast/{sample}_{type}.input",
    output:
        path=directory("snv_indels/mosaicforecast/{sample}_{type}"),
        phase="snv_indels/mosaicforecast/{sample}_{type}/all.phasing",
    params:
        extra=config.get("mosaicforecast", {}).get("extra", ""),
        path="alignment/samtools_merge_bam/",
    log:
        "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.benchmark.tsv",
            config.get("mosaicforecast", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosaicforecast", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosaicforecast", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosaicforecast", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosaicforecast", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosaicforecast", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosaicforecast", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosaicforecast", {}).get("container", config["default_container"])
    message:
        "{rule}: mosaicforecast evaluate candidate variants"
    shell:
        "(python /usr/local/bin/Phase.py "
        "{params.path} "
        "{output.path} "
        "{input.fasta} "
        "{input.variants} "
        "20 /usr/local/bin/k24.umap.wg.bw "
        " 4 bam "
        "{params.extra}) &> {log}"
