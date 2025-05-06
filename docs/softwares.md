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

## [deepmosaic_draw](https://github.com/shishenyxx/DeepMosaic)
DeepMosaic is a deep-learning-based mosaic single nucleotide classification tool without the need of matched control information. Firstly, feature extraction and visualization of the candidate mosaic variants (Visualization Module)

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepmosaic__deepmosaic_draw#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepmosaic__deepmosaic_draw#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepmosaic_draw#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepmosaic_draw#

---

## [deepmosaic_input]
Run python script to create the information txt-file needed in deepmosaic_draw. It contains sample_id, path to .bam and .vcf file, average depth and sex (predicted with Peddy).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepmosaic__deepmosaic_input#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepmosaic__deepmosaic_input#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepmosaic_input#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepmosaic_input#

---

## [deepsomatic_predict](https://github.com/shishenyxx/DeepMosaic)
DeepMosaic is a deep-learning-based mosaic single nucleotide classification tool without the need of matched control information. Secondly, prediction for mosaicism (Classification Module)

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepmosaic__deepmosaic_predict#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepmosaic__deepmosaic_predict#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepmosaic_predict#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepmosaic_predict#

---

## [deepsomatic_t_only](https://github.com/google/deepsomatic)
Using only tumor to call somatic SNV and indel variant from both short and long-read with DeepSomatic. The tool is the somatic version of DeepVariant. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepsomatic__deepsomatic_t_only#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepsomatic__deepsomatic_t_only#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepsomatic_t_only#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepsomatic_t_only#

---

## [deepsomatic_tn](https://github.com/google/deepsomatic)
Tumor/normal analysis to call somatic SNV and indel variant from both short and long-read with DeepSomatic. The tool is the somatic version of DeepVariant. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__deepsomatic__deepsomatic_tn#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__deepsomatic__deepsomatic_tn#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepsomatic_tn#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepsomatic_tn#

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

## [hiphase](https://github.com/PacificBiosciences/HiPhase)
Hiphase jointly phases small, structural, and tandem repeat variants for PacBio sequencing data

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__hiphase__hiphase#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__hiphase__hiphase#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__hiphase#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__hiphase#

---

## [merge_af_complex_variants](https://github.com/hydra-genetics/snv_indels/blob/develop/workflow/scripts/merge_af.py)
Python script for handling complex variants with several vcf record.

Some variant callers (e.g. vardict) will compose variants within close physical distance and report it as one complex variant.
vt_decompose will separate these complex variants into separate records. However, during decomposition the same allele might be
reported in several record, one originating from the single variant and one or more records reported from one or more complex variants,
even if these records are corresponding to the same allele at the same position. The allele frequencies will also be different for these
records since it might be derived from the frequency of the allele in combination with a specific allele at another position whithin the 
complex variant. 

This python script can turn several records from the same allele at the same position into one record. Depending on the method given by the user
the allele frequency and metrics derived from this will be reported differently. "skip" is the default method and in this case no alterations of 
the records will be made. All records for a complex variant will be returned in the output vcf. The method "max" will return the record with
the highest allele frequency and discard any additional records, with the same allele and position, from the output vcf. The "sum" method will
sum the allele frequencies from all records with the same allele and position. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__merge_af_complex_variants__merge_af_complex_variants#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__merge_af_complex_variants__merge_af_complex_variants#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__merge_af_complex_variants#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__merge_af_complex_variants#

---

## [mosaicforecast](https://github.com/parklab/MosaicForecast)
A machine learning method that leverages read-based phasing and read-level features to accurately detect mosaic SNVs (SNPs, small indels) from NGS data. It builds on existing algorithms to achieve a multifold increase in specificity.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__mosaicforecast__mosaicforecast#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__mosaicforecast__mosaicforecast#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__mosaicforecast#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__mosaicforecast#

---

## [mosaicforecast_input]
Getting the information needed from the .vcf-file for MosaicForecast. A list of the position that should be evaluated for their mosaic potential. 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__mosaicforecast__mosaicforecast_input#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__mosaicforecast__mosaicforecast_input#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__mosaicforecast_input#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__mosaicforecast_input#

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
