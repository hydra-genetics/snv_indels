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


rule mosaicforecast_genotype_prediction:
    input:
        features="snv_indels/mosaicforecast/{sample}_{type}/features.txt",
    output:
        predict="snv_indels/mosaicforecast/{sample}_{type}/SNP.predictions",  # SNP.predictions (Refine) or DEL.predictions (Phase) or INS.predictions (Phase)
    params:
        extra=config.get("mosaicforecast_genotype_prediction", {}).get("extra", ""),
        model_trained=config.get("mosaicforecast_genotype_prediction", {}).get("model_trained", ""),
        model_type=config.get("mosaicforecast_genotype_prediction", {}).get("model_type", ""),
    log:
        "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast_genotype_prediction.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast_genotype_prediction.benchmark.tsv",
            config.get("mosaicforecast_genotype_prediction", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosaicforecast_genotype_prediction", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosaicforecast_genotype_prediction", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosaicforecast_genotype_prediction", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
        partition=config.get("mosaicforecast_genotype_prediction", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosaicforecast_genotype_prediction", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosaicforecast_genotype_prediction", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosaicforecast_genotype_prediction", {}).get("container", config["default_container"])
    message:
        "{rule}: mosaicforecast predicts all input sites"
    shell:
        "(Prediction.R "
        "{input.features} "
        "{params.model_trained} "
        "{params.model_type} "
        "{output.predict} "
        "{params.extra}) &> {log}"


rule mosaicforecast_phasing:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
        variants="snv_indels/mosaicforecast/{sample}_{type}.input",
    output:
        path=directory("snv_indels/mosaicforecast/{sample}_{type}"),
        phase="snv_indels/mosaicforecast/{sample}_{type}/all.phasing",
    params:
        extra=config.get("mosaicforecast_phasing", {}).get("extra", ""),
        f_format=config.get("mosaicforecast_phasing", {}).get("f_format", ""),
        path=config.get("mosaicforecast_phasing", {}).get("path", ""),
        umap=config.get("mosaicforecast_phasing", {}).get("umap", ""),
    log:
        "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.benchmark.tsv",
            config.get("mosaicforecast_phasing", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosaicforecast_phasing", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosaicforecast_phasing", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosaicforecast_phasing", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosaicforecast_phasing", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosaicforecast_phasing", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosaicforecast_phasing", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosaicforecast_phasing", {}).get("container", config["default_container"])
    message:
        "{rule}: mosaicforecast phasing evaluation of candidate variants"
    shell:
        "(Phase.py "
        "{params.path} "
        "{output.path} "
        "{input.fasta} "
        "{input.variants} "
        "20 {params.umap} "
        "{resources.threads} "
        "{params.f_format} "
        "{params.extra}) &> {log}"


rule mosaicforecast_readlevel:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        fasta=config.get("reference", {}).get("fasta", ""),
        variants="snv_indels/mosaicforecast/{sample}_{type}.input",
    output:
        features="snv_indels/mosaicforecast/{sample}_{type}/features.txt",
    params:
        extra=config.get("mosaicforecast_readlevel", {}).get("extra", ""),
        f_format=config.get("mosaicforecast_readlevel", {}).get("f_format", ""),
        path=config.get("mosaicforecast_readlevel", {}).get("path", ""),
        umap=config.get("mosaicforecast_readlevel", {}).get("umap", ""),
    log:
        "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast/{sample}_{type}.mosaicforecast.vcf.benchmark.tsv",
            config.get("mosaicforecast_readlevel", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosaicforecast_readlevel", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosaicforecast_readlevel", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosaicforecast_readlevel", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosaicforecast_readlevel", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosaicforecast_readlevel", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosaicforecast_readlevel", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosaicforecast_readlevel", {}).get("container", config["default_container"])
    message:
        "{rule}: mosaicforecast extraction of read-level features"
    shell:
        "(ReadLevel_Features_extraction.py "
        "{input.variants} "
        "{output.features} "
        "{params.path} "
        "{input.fasta} "
        "{params.umap} "
        "{resources.threads} "
        "{params.f_format} "
        "{params.extra}) &> {log}"
