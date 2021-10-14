
rule bgzip_vcf:
    input:
        "{file}.vcf",
    output:
        temp("{file}.vcf.gz"),
    log:
        "{file}.vcf.gz.log",
    benchmark:
        repeat("{file}.vcf.gz.benchmark.tsv", config.get("bgzip_vcf", {}).get("benchmark_repeats", 1))
    threads: config.get("bgzip_vcf", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("bgzip_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/bgzip_vcf.yaml"
    message:
        "{rule}: bgzip {wildcards.file}.vcf"
    shell:
        "(bgzip -c {input} > {output}) &> {log}"
