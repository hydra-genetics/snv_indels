#!/bin/python3.6
import sys
import re
from pysam import VariantFile

vcf_in = VariantFile(snakemake.input[0])
method = snakemake.input[0].split("/")[1]
# Add new filter descriptions to new header
new_header = vcf_in.header

if method == "freebayes" or method == "pisces":  # Byta description on AF for freebayes?
    sample = vcf_in.header.samples[0]

if method == "pisces" or method == "mutect2":
    new_header.info.add("AF", "A", "Float", "DescriptionDescription")

if method == "snver":
    new_header.info.add("AF", "A", "Float", "Allel count divided on depth, crude")

if method == "varscan":
    new_header.info.add("AF", "A", "Float", "Allel count divided on depth (Quality of bases: Phred score >= 15)")


# start new vcf with the new_header
vcf_out = VariantFile(snakemake.output[0], 'w', header=new_header)

for record in vcf_in.fetch():
    if method == "freebayes":
        ad = record.samples[sample].get("AD")
        af = []
        af = [ad[1] / (ad[0] + ad[1])]
        if len(ad) > 2:
            for item in ad[2:]:
                af.append(item / sum(ad))
    if method == "pisces":
        af = record.samples[sample].get("VF")
    if method == "mutect2":
        af = record.samples[0].get("AF")
    if method == "snver":
        dp = record.info["DP"]
        ac = record.info["AC"]
        af = ac / dp
    if method == "varscan":
        ac = record.samples[0].get("AD")  # AD = AC
        dp = record.samples[0].get("DP")
        af = ac / dp
    if method == "vardict":
        af = record.samples[0].get("AF")

    record.info["AF"] = af

    vcf_out.write(record)
