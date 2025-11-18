__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosaicforecast_input:
    input:
        vcf="snv_indels/deepsomatic_t_only/{sample}_{type}.vcf.gz",
    output:
        variants=temp("snv_indels/mosaicforecast_input/{sample}_{type}.input"),
    params:
        extra=config.get("mosaicforecast_input", {}).get("extra", ""),
        name=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
    log:
        "snv_indels/mosaicforecast_input/{sample}_{type}.input.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast_input/{sample}_{type}.input_benchmark.tsv",
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
        "(bcftools query "
        "{params.extra} "
        "-f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t{params.name}\n' "
        "{input.vcf} "
        "> {output.variants} ) &> {log}"


rule mosaicforecast_genotype_prediction:
    input:
        features="snv_indels/mosaicforecast_readlevel/{sample}_{type}/features.txt",
    output:
        predict=temp("snv_indels/mosaicforecast_genotype_prediction/{sample}_{type}.{variant}.predictions"),
    params:
        extra=config.get("mosaicforecast_genotype_prediction", {}).get("extra", ""),
        model_trained=lambda wildcards: config.get("mosaicforecast_genotype_prediction", {}).get(
            f"model_trained_{wildcards.variant}", ""
        ),
        model_type=lambda wildcards: config.get("mosaicforecast_genotype_prediction", {}).get(
            f"model_type_{wildcards.variant}", ""
        ),
    log:
        "snv_indels/mosaicforecast_genotype_prediction/{sample}_{type}.{variant}.predictions.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast_genotype_prediction/{sample}_{type}.{variant}.predictions.benchmark.tsv",
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
    wildcard_constraints:
        variant="SNP|INS|DEL",
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
        variants="snv_indels/mosaicforecast_input/{sample}_{type}.input",
    output:
        all_infor_snps=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/all_candidates"),
        candidates=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/all.merged.inforSNPs.pos"),
        all_phasing=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/all.phasing"),
        m_infor_snps=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/multiple_inforSNPs.log"),
        phase_table=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/all.phasing_2by2"),
        table=temp("snv_indels/mosaicforecast_phasing/{sample}_{type}/all_2x2table"),
        tmpdir=temp(directory("snv_indels/mosaicforecast_readlevel/{sample}_{type}/tmp")),
    params:
        extra=config.get("mosaicforecast_phasing", {}).get("extra", ""),
        f_format=config.get("mosaicforecast_phasing", {}).get("f_format", ""),
        path=lambda w, input: os.path.dirname(input[0]),
        umap=config.get("mosaicforecast_phasing", {}).get("umap", ""),
    log:
        "snv_indels/mosaicforecast_phasing/{sample}_{type}.mosaicforecast_phasing.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast_phasing/{sample}_{type}.mosaicforecast_phasing.benchmark.tsv",
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
        variants="snv_indels/mosaicforecast_input/{sample}_{type}.input",
    output:
        features=temp("snv_indels/mosaicforecast_readlevel/{sample}_{type}/features.txt"),
        features_tmp=temp("snv_indels/mosaicforecast_readlevel/{sample}_{type}/features.txt.tmp"),
    params:
        extra=config.get("mosaicforecast_readlevel", {}).get("extra", ""),
        f_format=config.get("mosaicforecast_readlevel", {}).get("f_format", ""),
        path=lambda w, input: os.path.dirname(input[0]),
        umap=config.get("mosaicforecast_readlevel", {}).get("umap", ""),
    log:
        "snv_indels/mosaicforecast_readlevel/{sample}_{type}.mosaicforecast_readlevel.log",
    benchmark:
        repeat(
            "snv_indels/mosaicforecast_readlevel/{sample}_{type}.mosaicforecast_readlevel.benchmark.tsv",
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
