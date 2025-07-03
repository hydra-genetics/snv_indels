<p align="center">
<a href="https://hydra-genetics-snv-indels.readthedocs.io">https://hydra-genetics-snv-indels.readthedocs.io</a>
</p>

# :snake: hydra-genetics/snv_indels

#### Snakemake module containing steps to call snv and small indels

![Lint](https://github.com/hydra-genetics/snv_indels/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/snv_indels/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![pycodestyle](https://github.com/hydra-genetics/snv_indels/actions/workflows/pycodestyle.yaml/badge.svg?branch=develop)
![pytest](https://github.com/hydra-genetics/snv_indels/actions/workflows/pytest.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/hydra-genetics/prealignment/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)
![integration test](https://github.com/hydra-genetics/prealignment/actions/workflows/integration.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

The module contains rules to call variants from `.bam`-files per chromosome, merging
the resulting `.vcf`-files, fixing the allele frequency field followed by decomposing
and normalizing steps to finally combine the results from different callers using
an ensemble approach. Available callers in the standard setup are [Mutect2](https://gatk.broadinstitute.
org/hc/en-us/articles/360037593851-Mutect2),
[Freebayes](https://github.com/freebayes/freebayes), [VarDict](https://github.com/AstraZeneca-NGS/VarDict) and 
[Haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).
Mutect2 is also used to generate a genomic `.vcf`-file.

Other variant callers can be added by providing caller-specific configurations: 
[DeepVariant](https://github.com/google/deepvariant), 
[DeepSomatic](https://github.com/google/deepsomatic), and 
[ClairS-TO](https://github.com/HKU-BAL/ClairS-TO).
The output of those variant callers can also be processed with decomposition and normalization steps if the user 
wishes so.
Moreover, the output of DeepSomatic in tumor-only settings can be piped into [DeepMosaic](https://github.com/XiaoxuYangLab/DeepMosaic) and 
[MosaicForecast](https://github.com/parklab/MosaicForecast) to identify mosaic variants.


## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-v0.9.2-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-1.3.1-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/)
[![snakemake](https://img.shields.io/badge/snakemake-6.10.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)

*Note! Releases of snv_indels <= v0.2.0 needs tabulate<0.9.0 added in requirements.txt**

## :school_satchel: Preparations

### Sample and unit data

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/prealignment/blob/develop/config/units.tsv).
The following information need to be added to these files:

#### Short-read data

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| tumor_content | ratio of tumor cells to total cells |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

#### Long-read data

| Column Id           | Description                                                                                                  |
|---------------------|--------------------------------------------------------------------------------------------------------------|
| **`samples.tsv`**   |
| sample              | unique sample/patient id, one per row                                                                        |
| (tumor_content)     | ratio of tumor cells to total cells (not relevant if data for constitutional disease analysis)               |
| (sex)               | sex prediction from somalier (relevant if data for constitutional disease analysis)                          |
| (trioid)            | ID of the trio (relevant if data for constitutional disease analysis)                                        |
| (trio_member)       | code of the sample in the trio (relevant if data for constitutional disease analysis)                        |
| **`units.tsv`**     |
| sample              | same sample/patient id as in `samples.tsv`                                                                   |
| type                | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA                              |
| platform            | type of sequencing platform, e.g. `PACBIO` or `ONT`                                                          |
| machine             | specific machine id, e.g. PacBio instrument `REVIO`                                                          |
| processing_unit     | processing unit id, e.g. `@Axxxxx`                                                                           |
| (run_id)            | run id of the sequencing run, e.g. `run1` (not relevant if PacBio data)                                      |
| barcode             | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC`                 |
| methylation         | whether methylation is called, e.g. `true` or `false`                                                        |
| (basecalling_model) | basecalling model used, e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` for ONT data (not relevant if PacBio data) |
| bam                 | absolute path to the directory with the `.bam` file, e.g. `/path/to/bam`                                     |

### Reference data

A reference `.fasta`-file should be specified in `config.yaml` in the section `reference` and `fasta`.
In addition, the file should be indexed using `samtools faidx` and the path of the resulting
file added to the stanza `fai`. A bed file containing the covered regions shall be added
to `design_bed`.

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
$ cd .tests/integration
$ snakemake -s ../../workflow/Snakefile -j1 --configfile <config_caller>.yaml --use-singularity --singularity-args " --cleanenv"
```

To test the standard pipeline for short-read data, replace `<config_caller>` with `config.yaml`.
`<config_caller>` should be replaced with the name of the caller-specific configuration file 
if you want to use long-read and/or deep-learning-based variant callers, e.g. `config_deepvariant.yaml`.
Note that using a caller-specific configuration file will turn off the standard pipeline that uses short-read 
variant callers like Mutect2, VarDict and Freebayes.

## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module snv_indels:
    snakefile:
        github(
            "hydra-genetics/snv_indels",
            path="workflow/Snakefile",
            tag="v0.1.0",
        )
    config:
        config


use rule * from snv_indels as snv_indels_*
```

### Compatibility

Latest:
 - alignment:v0.5.1

 See [COMPATIBLITY.md](../master/COMPATIBLITY.md) file for a complete list of module compatibility.

### Output files

The following output files should be targeted via another rule:

| File | Description |
|---|---|
| `snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vcf.gz` | combined `.vcf` generated by ensemble |
| `snv_indels/{caller}/{sample}_{type}.normalized.sorted.vcf.gz` | sorted `.vcf.gz` for each caller |
| `snv_indels/gatk_mutect2_gvcf/{sample}_{type}.merged.g.vcf.gz` | genomic `.vcf` |
| `snv_indels/deepvariant/{sample}_{type}.merged.vcf.gz` | deepvariant `.vcf.gz` |
| `snv_indels/deepvariant/{sample}_{type}.merged.g.vcf.gz` | genomic `.g.vcf.gz` for deepvariant |
| `snv_indels/glnexus/{sample}_{type}.vcf.gz` | trio `.vcf.gz` with proband sample id generated by glnexus from deeptrio genomic `.g.vcf.gz` |

The following output files require to add a caller-specific configurations for the pipeline:

| File                                                          | Description                                                                  |
|---------------------------------------------------------------|------------------------------------------------------------------------------|
| `snv_indels/deepsomatic_tn/{sample}_{type}.vcf.gz`            | `.vcf.gz` with called variants by deepsomatic in tumor-normal samples        |
| `snv_indels/deepsomatic_t_only/{sample}_{type}.vcf.gz`        | `.vcf.gz` with called variants by deepsomatic in tumor-only samples          |
| `snv_indels/deepmosaic/{sample}_{type}/final_predictions.txt` | list of called variants that are identified as mosaicisms by deepmosaic      |
| `snv_indels/mosaicforecast/{sample}_{type}/all.phasing`       | list of called variants that are identified as mosaicisms by mosaicforecast  |
| `snv_indels/clairs_to/{sample}_{type}.snv-indels.vcf.gz`      | `.vcf.gz` with called variants by deepsomatic in tumor-only samples          |
| `snv_indels/deepvariant/{sample}_{type}.merged.vcf.gz`        | deepvariant `.vcf.gz` for PacBio data if `deepvariant: model_type: "PACBIO"` |
| `snv_indels/deepvariant/{sample}_{type}.merged.g.vcf.gz`      | genomic `.g.vcf.gz` for PacBio data if `deepvariant: model_type: "PACBIO"`   |


## :judge: Rule Graph

### Tools suitable for short-read data
![rule_graph_sr](docs/images/snv_indels_SR.png)

### Tools suitable for long-read data
![rule_graph_lr](docs/images/snv_indels_LR.png)