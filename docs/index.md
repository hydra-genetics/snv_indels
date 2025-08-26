# [SNV_indels module](https://github.com/hydra-genetics/snv_indels)

The **snv_indels** module consists of a number of variant callers for:
* Short-read (SR) data, both somatic and germline variants.
* Long-read (LR) data, both somatic and germline variants. Most of the variant callers for LR data are based on 
  deep-learning models.
Some tools are able to handle data from tumor-only (TO) samples, while others require matched normal (MN) samples.

Somatic calls can be made by:
- freebayes, mutect2, vardict,
- deepsomatic, clairs-to.

Germline calls can be produced by:
- haplotypecaller, deepvariant, deeptrio, glnexus.


The module also provide tools for:
* variant decomposition and normalization with _VT_ (applicable to both SR and LR data),
* aggregation of the results from different callers into an ensemble with *bc_bio* (only for SR data),
* sorting and indexing of the variants with _bcftools_ (applicable to both SR and LR data),
* predict which somatic variants are likely to be mosaicisms with _deepmosaic_ or _mosaicforecast_ (applicable to both SR and LR data).

The module is designed to be used in a snakemake pipeline and can be easily integrated into existing pipelines that
follow the logic of [**hydra-genetics**](https://hydra-genetics.readthedocs.io/en/latest/).

# [Hydra-genetics](https://hydra-genetics.readthedocs.io/en/latest/)

We are an organization/community with the goal of making [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, a bit more structured and of higher quality.

We do this by providing [snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) that can be combined to create a complete analysis or included in already existing pipelines. All modules are subjected to extensive testing to make sure that new releases doesn't unexpectedly break existing pipeline or deviate from guidelines and best practices on how to write code.

<p align="center" width="100%">
    <img width="10%" src="images/hydragenetics.png">
</p>