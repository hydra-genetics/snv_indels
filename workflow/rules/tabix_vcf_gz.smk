
rule tabix_vcf_gz:
    input:
        "snv_indels/{caller}/{sample}_{type}_{chr}.unfilt.vcf.gz",
    output:
        temp("snv_indels/{caller}/{sample}_{type}_{chr}.unfilt.vcf.gz.tbi"),
    params:
        config.get("tabix_vcf_gz", {}).get("extra", ""),
    log:
        "snv_indels/{caller}/{sample}_{type}_{chr}.unfilt.vcf.gz.log",
    benchmark:
        repeat("tabix_vcf_gz/{caller}/{sample}_{type}_{chr}.unfilt.vcf.gz.benchmark.tsv", config.get("tabix_vcf_gz", {}).get("benchmark_repeats", 1))
    threads: config.get("tabix_vcf_gz", config["default_resources"]).get("threads", config["default_resources"]["threads"])
    container:
        config.get("tabix_vcf_gz", {}).get("container", config["default_container"])
    conda:
        "../envs/tabix_vcf_gz.yaml"
    message:
        "{rule}: Tabix index snv_indels/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.unfilt.vcf.gz"
    wrapper:
        "0.79.0/bio/tabix"
