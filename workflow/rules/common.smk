__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2021, Patrik Smeds"
__email__ = "patrik.smeds@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd

from hydra_genetics.utils.misc import extract_chr
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate

min_version("7.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    barcode="[A-Z+]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",
    vcf="vcf|g.vcf|unfiltered.vcf",


def get_bvre_params_sort_order(wildcards: snakemake.io.Wildcards):
    return ",".join(config.get("bcbio_variation_recall_ensemble", {}).get("callers", ""))


def get_java_opts(wildcards: snakemake.io.Wildcards):
    java_opts = config.get("haplotypecaller", {}).get("java_opts", "")
    if "-Xmx" in java_opts:
        raise WorkflowError("You are not allowed to use -Xmx in java_opts. Set mem_mb in resources instead.")
    java_opts += "-Xmx{}m".format(config.get("haplotypecaller", {}).get("mem_mb", config["default_resources"]["mem_mb"]))
    return java_opts


def get_gatk_mutect2_extra(wildcards: snakemake.io.Wildcards, name: str):
    extra = "{} {}".format(
        config.get(name, {}).get("extra", ""),
        "--intervals snv_indels/bed_split/design_bedfile_{}.bed".format(
            wildcards.chr,
        ),
    )
    if name == "gatk_mutect2":
        extra = "{} {}".format(
            extra,
            "--f1r2-tar-gz snv_indels/gatk_mutect2/{}_{}_{}.unfiltered.f1r2.tar.gz".format(
                wildcards.sample,
                wildcards.type,
                wildcards.chr,
            ),
        )
    if name == "gatk_mutect2_gvcf":
        extra = "{} {}".format(extra, "-ERC BP_RESOLUTION")
    return extra


def combine_extra_args(extra_args: dict):
    args = []
    for key in sorted(extra_args):
        value = extra_args[key]
        if value is None:
            continue
        if isinstance(value, bool):
            added_arg = "" if value else "no"
            added_arg += key
            args.extend(["--{}".format(added_arg)])
        else:
            args.extend(["--{} {}".format(key, value)])

    args_str = " ".join(args)

    return args_str


def get_make_example_args(wildcards: snakemake.io.Wildcards, output: list, name: str, vcf: str):

    model_type = config.get(name, {}).get("model", "WGS")
    special_args = {}
    if model_type == "WGS" or model_type == "WES":
        special_args["channels"] = "insert_size"
    elif model_type == "PACBIO":
        special_args = {}
        special_args["add_hp_channel"] = True
        special_args["alt_aligned_pileup"] = "diff_channels"
        special_args["max_reads_per_partition"] = 600
        special_args["min_mapping_quality"] = 1
        special_args["parse_sam_aux_fields"] = True
        special_args["partition_size"] = 25000
        special_args["phase_reads"] = True
        special_args["pileup_image_width"] = 199
        special_args["realign_reads"] = False
        special_args["sort_by_haplotypes"] = True
        special_args["track_ref_reads"] = True
        special_args["vsc_min_fraction_indels"] = 0.12

    special_args_str = combine_extra_args(special_args)
    extra = "{} {}".format(config.get(name, {}).get("extra", ""), special_args_str)

    if vcf == "gvcf":
        threads = config.get(name, {}).get("threads", config["default_resources"]["threads"])
        gvcf_path = " --gvcf {}/gvcf.tfrecord@{}.gz".format(output[0], threads)
        extra = "{} {}".format(extra, gvcf_path)

    return extra


def get_examples_infile(wildcards: snakemake.io.Wildcards, input: snakemake.io.Namedlist, name: str):
    threads = config.get(name, {}).get("threads", config["default_resources"]["threads"])
    examples_infile = "{}/make_examples.tfrecord@{}.gz".format(input.examples_dir, threads)

    return examples_infile


def get_postprocess_variants_args(
    wildcards: snakemake.io.Wildcards, input: snakemake.io.Namedlist, output: snakemake.io.Namedlist, me_config: str, extra: str
):

    if len(output) == 2:
        threads = me_config.get(name, {}).get("threads", config["default_resources"]["threads"])
        gvcf_tfrecord = "{}/gvcf.tfrecord@{}.gz".format(input[0], threads)
        gvcf_in = "--nonvariant_site_tfrecord_path {}".format(gvcf_tfrecord)
        gvcf_out = " --gvcf_outfile {}".format(output.gvcf)
        extra = "{} {} {}".format(extra, gvcf_in, gvcf_out)

    return extra


def compile_output_list(wildcards: snakemake.io.Wildcards):
    files = {
        "bcbio_variation_recall_ensemble": [
            "ensembled.vcf.gz",
        ],
        "gatk_mutect2_gvcf": [
            "merged.g.vcf.gz",
        ],
        "haplotypecaller": [
            "normalized.sorted.vcf.gz",
        ],
        "deepvariant": [
            "merged.vcf.gz",
            "merged.g.vcf.gz",
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
