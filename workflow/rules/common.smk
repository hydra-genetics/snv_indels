__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd

from hydra_genetics.utils.misc import extract_chr
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from snakemake.utils import min_version
from snakemake.utils import validate

min_version("6.10.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def get_bvre_params_sort_order(wildcards: snakemake.io.Wildcards):
    return ",".join(config.get("bcbio_variation_recall_ensemble", {}).get("callers", ""))


def get_mutect2_extra(wildcards: snakemake.io.Wildcards, name: str):
    extra = "{} {}".format(
        config.get(name, {}).get("extra", ""),
        "--intervals snv_indels/bed_split/design_bedfile_{}.bed".format(
            wildcards.chr,
        ),
    )
    if name == "mutect2":
        extra = "{} {}".format(
            extra,
            "--f1r2-tar-gz snv_indels/mutect2/{}_{}_{}.f1r2.tar.gz".format(
                wildcards.sample,
                wildcards.type,
                wildcards.chr,
            ),
        )
    if name == "mutect2_gvcf":
        extra = "{} {}".format(extra, "-ERC BP_RESOLUTION")
    return extra


def compile_output_list(wildcards: snakemake.io.Wildcards):
    files = {
        "bcbio_variation_recall_ensemble": [
            "ensembled.vcf.gz",
        ],
        "mutect2_gvcf": [
            "merged.vcf.gz",
        ],
    }
    output_files = [
        "snv_indels/%s/%s_%s.%s" % (prefix, sample, t, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    return output_files
