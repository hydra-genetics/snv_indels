import sys
from pysam import VariantFile

vcf_in = VariantFile(snakemake.input.vcf)
vcf_out = VariantFile(snakemake.output.vcf, 'w', header=vcf_in.header)
pass_filters = snakemake.params.pass_filters

for record in vcf_in.fetch():
    filter_list = record.filter.keys()
    found_filter = False
    for filter in filter_list:
        if filter in pass_filters:
            found_filter = True
    if found_filter:
        vcf_out.write(record)
