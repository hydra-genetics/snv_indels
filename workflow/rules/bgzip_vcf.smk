
rule bgzip_vcf:
    input:
        "snv_indels/{file}.vcf",
    output:
        temp("snv_indels/{file}.vcf.gz"),
    log:
        "snv_indels/{file}.vcf.gz.log",
    benchmark:
        repeat("snv_indels/{file}.vcf.gz.benchmark.tsv", config.get("bgzip_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("bgzip_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("bgzip_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/bgzip_vcf.yaml"
    message:
        "{rule}: bgzip snv_indels/{wildcards.file}.vcf"
    shell:
        "(bgzip -c {input} > {output}) &> {log}"
