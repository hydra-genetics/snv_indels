#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import pandas as pd
import re 
import sys

if not os.path.exists("snv_indels/deepmosaic"):
    os.makedirs("snv_indels/deepmosaic")

#file=snakemake.output.txt
file=sys.argv[3]
with open(file, "w") as output:
    output.write("#sample_name\tbam\tvcf\tdepth\tsex\n")

#bam_path=snakemake.input.bam
#vcf_path=snakemake.input.vcf
#sample=snakemake.params.name

bam_path=sys.argv[1]
vcf_path=sys.argv[2]
sample=sys.argv[4]

path="{PWD}"
if re.search("TE", path):
    depth=100
elif re.search("TC", path):
    depth=200
else:
    depth=40

sex_file=sys.argv[5]
sex=pd.read_csv(sex_file, sep = ',')

output=open(file, "a")

for s in range(len(sex.get("sample_id"))):
    if re.search(sex.get("sample_id")[s], sample):
        output.write(str(sample) + "\t" + str(bam_path) + "\t" + str(vcf_path) + "\t" + str(depth) + "\t" + str(sex.get("predicted_sex")[s]) + "\n")

output.close()
