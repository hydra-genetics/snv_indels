# [SNV_indels module](https://github.com/hydra-genetics/snv_indels)

The snv_indels module consists of a number of small variant callers for short read data of both somatic and germline samples. There is also tools for variant decomposition (vt), normalization (vt), ensemble (bcbio) and sorting (bcftools). Somatic calls can be called by freebayes, mutect2, and vardict while germline calls can be produced by haplotypecaller, deepvariant, deeptrio, glnexus.

# [Hydra-genetics](https://hydra-genetics.readthedocs.io/en/latest/)

We are an organization/community with the goal of making [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, a bit more structured and of higher quality.

We do this by providing [snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) that can be combined to create a complete analysis or included in already existing pipelines. All modules are subjected to extensive testing to make sure that new releases doesn't unexpectedly break existing pipeline or deviate from guidelines and best practices on how to write code.

<p align="center" width="100%">
    <img width="10%" src="images/hydragenetics.png">
</p>