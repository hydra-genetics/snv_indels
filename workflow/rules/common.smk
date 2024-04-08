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


def get_gvcf_output(wildcards, name):
    if config.get(name, {}).get("output_gvcf", False):
        return f" --output_gvcf snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}_{wildcards.chr}.g.vcf.gz "
    else:
        return ""


def get_parent_bams(wildcards):
    bam_path = "alignment/samtools_merge_bam"
    proband_sample = samples[samples.index == wildcards.sample]
    trio_id = proband_sample.at[wildcards.sample, "trioid"]

    mother_sample = samples[(samples.trio_member == "mother") & (samples.trioid == trio_id)].index[0]
    mother_bam = "{}/{}_{}.bam".format(bam_path, mother_sample, list(get_unit_types(units, mother_sample))[0])

    father_sample = samples[(samples.trio_member == "father") & (samples.trioid == trio_id)].index[0]
    father_bam = "{}/{}_{}.bam".format(bam_path, father_sample, list(get_unit_types(units, father_sample))[0])

    bam_list = [mother_bam, father_bam]

    bam_list += [f"{bam}.bai" for bam in bam_list]

    return bam_list


def get_make_examples_tfrecord(
    wildcards: snakemake.io.Wildcards, input: snakemake.io.Namedlist, nshards: int, program="deepvariant"
):
    examples_path = os.path.split(input[0])[0]

    if program == "deepvariant":
        examples_tfrecord = "{}/make_examples.tfrecord@{}.gz".format(examples_path, nshards)
    elif program == "deeptrio":
        examples_tfrecord = "{}/make_examples_{}.tfrecord@{}.gz".format(examples_path, wildcards.trio_member, nshards)

    return examples_tfrecord


def get_deeptrio_model(wildcards):
    models_config = config.get("deeptrio_call_variants", {}).get("model", "")
    if wildcards.trio_member in ["parent1", "parent2"]:
        model_file = models_config.get("parent", "")
    else:
        model_file = models_config.get("child", "")

    return model_file


def deeptrio_postprocess_variants_args(
    wildcards: snakemake.io.Wildcards, input: snakemake.io.Namedlist, me_config: str, extra: str
):
    me_path = os.path.split(input.call_variants_record)[0]
    nshards = config.get(me_config).get("n_shards", 2)
    gvcf_tfrecord = "{}/gvcf_{}.tfrecord@{}.gz".format(me_path, wildcards.trio_member, nshards)
    gvcf_in = "--nonvariant_site_tfrecord_path {}".format(gvcf_tfrecord)
    extra = "{} {} ".format(extra, gvcf_in)

    return extra


def get_glnexus_input(wildcards, input):
    gvcf_input = "-i {}".format(" -i ".join(input.gvcfs))

    return gvcf_input


def compile_output_list(wildcards: snakemake.io.Wildcards):
    if config["run_deepvariant"]:
        files = {
            "deepvariant": [
                "merged.vcf.gz",
                "merged.g.vcf.gz",
            ],
        }
        output_files = [
            "snv_indels/%s/%s_%s.%s" % (prefix, sample, t, suffix)
            for prefix in files.keys()
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for t in get_unit_types(units, sample)
            for suffix in files[prefix]
        ]

        files = {
            "glnexus": ["vcf.gz"],
        }
        output_files += [
            "snv_indels/%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
            for prefix in files.keys()
            for sample in samples[samples.trio_member == "proband"].index
            for unit_type in get_unit_types(units, sample)
            for suffix in files[prefix]
        ]
    else:
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
        }
        output_files = [
            "snv_indels/%s/%s_%s.%s" % (prefix, sample, t, suffix)
            for prefix in files.keys()
            for sample in get_samples(samples[pd.isnull(samples["trioid"])])
            for t in get_unit_types(units, sample)
            for suffix in files[prefix]
        ]

    return output_files
