# Softwares used in the snv_indels module

## [bcbio_variation_recall_ensemble](https://github.com/bcbio/bcbio.variation.recall)
Ensemble variants in the `.vcf` file format from different callers in one `.vcf` file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcbio_variation_recall__bcbio_variation_recall_ensemble#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcbio_variation_recall__bcbio_variation_recall_ensemble#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcbio_variation_recall_ensemble#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcbio_variation_recall_ensemble#

---

## [bcftools_concat](https://samtools.github.io/bcftools/bcftools.html)
Concatenate `.vcf` files from different chromosomes using bcftools concat.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_concat#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_concat#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_concat#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_concat#

---

## [bcftools_sort](https://samtools.github.io/bcftools/bcftools.html)
Sort a `.vcf` file using bcfools sort.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_sort#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_sort#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_sort#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_sort#

---

## [bcftools_view](https://samtools.github.io/bcftools/bcftools.html)
Convert bcftools `.bcf` binary file format to a `.vcf` file using bcftools.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_view#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_view#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_view#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_view#

---

## [bed_split](https://linux.die.net/man/1/awk)
Split a bed file into chromosomes using an awk command going by the first column (chromosome).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bed_split__bed_split#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bed_split__bed_split#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bed_split#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bed_split#

---

## [bgzip](http://www.htslib.org/doc/bgzip.html)
Compress a `.vcf` file using bgzip.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bgzip__bgzip#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bgzip__bgzip#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bgzip#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bgzip#

---

## [deeptrio_call_variants](https://github.com/google/deepvariant)
Step 2 of 3 in the calling of SNVs and INDELs using deeptrio.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deeptrio__deeptrio_call_variants#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deeptrio__deeptrio_call_variants#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deeptrio_call_variants#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deeptrio_call_variants#

---

## [deeptrio_make_examples](https://github.com/google/deepvariant)
Step 1 of 3 in the calling of SNVs and INDELs using deeptrio.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deeptrio__deeptrio_make_examples#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deeptrio__deeptrio_make_examples#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deeptrio_make_examples#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deeptrio_make_examples#

---

## [deeptrio_postprocess_variants](https://github.com/google/deepvariant)
Step 3 of 3 in the calling of SNVs and INDELs using deeptrio.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deeptrio__deeptrio_postprocess_variants#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deeptrio__deeptrio_postprocess_variants#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deeptrio_postprocess_variants#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deeptrio_postprocess_variants#

---

## [deepvariant](https://github.com/google/deepvariant)
Calling of SNVs and INDELs using deepvariant.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepvariant__deepvariant#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepvariant__deepvariant#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepvariant#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepvariant#

---
## [fix_af](https://github.com/hydra-genetics/snv_indels/blob/develop/workflow/scripts/fix_af.py)
Python script that fixes missing AF in the format field of a vcf file as this is not generated by all callers.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__fix_af__fix_af#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__fix_af__fix_af#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__fix_af#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__fix_af#

---

## [freebayes](https://github.com/freebayes/freebayes)
Somatic variant caller for SNVs and small indels.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__freebayes__freebayes#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__freebayes__freebayes#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__freebayes#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__freebayes#

---

## [gatk_mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2)
Step 1 of 4 of Mutect2 variant calling generating chromosome split unfiltered calls as well as variant statistic files (used for variant filtering) and INDEL-bam files.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_mutect2#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_mutect2#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_mutect2#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_mutect2#

---

## [gatk_mutect2_gvcf](https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2)
Mutect2 is used to generate a genome vcf file containing allele information for all positions in the bed file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_mutect2_gvcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_mutect2_gvcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_mutect2_gvcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_mutect2_gvcf#

---

## [gatk_mutect2_filter](https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2)
Step 3 of 4 of Mutect2 variant calling filtering the called variants using the merged statistics file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_mutect2_filter#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_mutect2_filter#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_mutect2_filter#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_mutect2_filter#

---

## [gatk_mutect2_merge_stats](https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2)
Step 2 of 4 of Mutect2 variant calling merging the statistics file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__gatk__gatk_mutect2_merge_stats#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__gatk__gatk_mutect2_merge_stats#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__gatk_mutect2_merge_stats#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__gatk_mutect2_merge_stats#

---

## [glnexus](https://github.com/dnanexus-rnd/GLnexus)
Joint variant caller based on `.g.vcf` files from deeptrio or other sources. Can be used for family trios as well as larger cohorts.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__glnexus__glnexus#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__glnexus__glnexus#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__glnexus#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__glnexus#

---

## [haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller)
Germline variant caller for SNVs and INDELs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__haplotypecaller__haplotypecaller#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__haplotypecaller__haplotypecaller#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__haplotypecaller#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__haplotypecaller#

---

## [mutect2_pass_filter](https://github.com/hydra-genetics/snv_indels/blob/develop/workflow/scripts/mutect2_pass_filter.py)
Step 4 of 4 of Mutect2 somatic variant calling. A python script that hard filters the soft filtered vcf file from Mutect2 based on the FILTER column.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__haplotypecaller__haplotypecaller#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__haplotypecaller__haplotypecaller#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__haplotypecaller#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__haplotypecaller#

---

## [paraphase](https://github.com/PacificBiosciences/paraphase)
Many medically relevant genes fall into 'dark' regions where variant calling is limited due to high sequence homology with paralogs or pseudogenes. Paraphase is a Python tool that takes HiFi aligned BAMs as input (whole-genome or enrichment), phases haplotypes for genes of the same family, determines copy numbers and makes phased variant calls for these genes.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__paraphase__paraphase#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__paraphase__paraphase#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__paraphase#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__paraphase#

---

## [tabix](http://www.htslib.org/doc/tabix.html)
Creates an index file for faster processing of positions in a bgzipped vcf file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__tabix__tabix#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__tabix__tabix#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__tabix#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__tabix#

---

## [vardict](https://github.com/AstraZeneca-NGS/VarDict)
Somatic variant caller for SNVs and INDELs.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__vardict__vardict#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__vardict__vardict#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__vardict#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__vardict#

---

## [vt_decompose](https://github.com/atks/vt)
Decomposition of vcf files. Uses both the decompose and the decompose_blocksub command. Decompose_blocksub divide clumped variants into separate records, while decompose divide multiallelic variants into separate records.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__vt__vt_decompose#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__vt__vt_decompose#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__vt_decompose#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__vt_decompose#

---

## [vt_normalize](https://github.com/atks/vt)
Normalization of vcf files. Left aligns INDELs and adds one reference allele. Affects variant position.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__vt__vt_normalize#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__vt__vt_normalize#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__vt_normalize#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__vt_normalize#

---

## [whatshap_haplotag](https://github.com/whatshap/whatshap/)
WhatsHap is a software for phasing genomic variants using DNA sequencing reads, also called read-based phasing or haplotype assembly. It is especially suitable for long reads, but works also well with short reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__whatshap__whatshap_haplotag#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__whatshap__whatshap_haplotag#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__whatshap_haplotag#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__whatshap_haplotag#

---

## [whatshap_phase](https://github.com/whatshap/whatshap/)
WhatsHap is a software for phasing genomic variants using DNA sequencing reads, also called read-based phasing or haplotype assembly. It is especially suitable for long reads, but works also well with short reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__whatshap__whatshap_phase#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__whatshap__whatshap_phase#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__whatshap_phase#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__whatshap_phase#
