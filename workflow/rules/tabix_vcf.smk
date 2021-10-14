
rule tabix_vcf:
    input:
        "snv_indels/{caller}/{file}.vcf.gz",
    output:
        temp("snv_indels/{caller}/{file}.vcf.gz.tbi"),
    params:
        config.get("tabix_vcf", {}).get("extra", ""),
    log:
        "snv_indels/{caller}/{file}.unfilt.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{caller}/{file}.unfilt.vcf.gz.benchmark.tsv",
            config.get("tabix_vcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("tabix_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("tabix_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/tabix_vcf.yaml"
    message:
        "{rule}: Tabix index snv_indels/{rule}/{wildcards.file}.vcf.gz"
    wrapper:
        "0.79.0/bio/tabix"
