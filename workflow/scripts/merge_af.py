#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Monika Brandt"
__copyright__ = "Copyright 2024, Monika Brandt"
__email__ = "monika.brandt@scilifelab.uu.se"
__license__ = "GPL-3"

import logging
import pysam
import statistics
import shutil
from collections import Counter


def merge_records_complex_positions(vcf, method):
    '''
    Parse a vcf file with simple and decomposed complex variants.
    Dependent on the method given, select records for the same complex
    variant and report it back as one record with an updated allele frequency.

    param vcf: Input vcf file with all complex variants decomposed into
    several records.
    param method: str method to use, max or sum, for calculation of AF
    return record_dict: A dictionary following the format,
    {"<chr>-<pos>-<REF allele>-<ALT allele>": [<pysam.VariantRecord object>]},
    where key is based on the variant and the value is a list with only one element/record.
    return vcf.header: The header of the input vcf with the new info key,
    "COMPLEXAF" describing the method used to select AF for complex variants.
    '''

    # Add info field COMPLEXAF to header
    vcf.header.info.add("COMPLEXAF", "1", "String", "Method used to select AF for complex variants max or sum.")

    # Add info field TYPE to header
    if "TYPE" not in vcf.header.info:
        vcf.header.info.add("TYPE", "1", "String", "Variant Type: SNV Insertion Deletion Complex")

    chr_pos_ref_alt = []
    record_dict = {}

    for record in vcf.fetch():
        var_key = "{}-{}-{}-{}".format(record.chrom, record.pos, record.ref, record.alts[0])
        chr_pos_ref_alt.append(var_key)
        if var_key in record_dict:
            record_dict[var_key].append(record)
        else:
            record_list = [record]
            record_dict[var_key] = record_list

    # Select variants (on chromosome, position, ref.allele and alt. allele) that occurs in more than one record.
    complex_positions_list = [item for item, count in Counter(chr_pos_ref_alt).items() if count > 1]

    for var_key in complex_positions_list:
        # if method "max" is given by the user the code will return the record with the highest allele frequency
        # as is without any updates on for instance AF, VF or AD.
        # A new info key will be added: COMPLEXAF=max
        if method == "max":
            complex_AF = max(record.info["AF"][0] for record in record_dict[var_key])
            merged_record = record_dict[var_key][[record.info["AF"][0] == complex_AF for record in record_dict[var_key]].index(True)]  # noqa
            merged_record.info["COMPLEXAF"] = method
            merged_record.info["TYPE"] = "Complex"

        # if method "sum" is given by the user, return one record for the given variant where
        # AF is the sum of AF for all records for the variant.
        # A new info key will be added: COMPLEXAF=sum
        if method == "sum":
            # Pick first record in list to use as a template for updates.
            merged_record = record_dict[var_key][0]
            merged_record.info["COMPLEXAF"] = method
            merged_record.info["TYPE"] = "Complex"

            complex_AF = sum(record.info["AF"][0] for record in record_dict[var_key])
            merged_record.info["AF"] = complex_AF

            try:
                complex_VD = sum(record.info["VD"] for record in record_dict[var_key])
                merged_record.info["VD"] = complex_VD
            except KeyError:
                continue

            try:
                complex_HIAF = sum(record.info["HIAF"] for record in record_dict[var_key])
                merged_record.info["HIAF"] = complex_HIAF
            except KeyError:
                continue

            try:
                complex_HICNT = sum(record.info["HICNT"] for record in record_dict[var_key])
                merged_record.info["HICNT"] = complex_HICNT
            except KeyError:
                continue

            try:
                complex_QUAL = statistics.mean(record.info["QUAL"] for record in record_dict[var_key])
                merged_record.info["QUAL"] = complex_QUAL
            except KeyError:
                continue

            try:
                complex_MQ = statistics.mean(record.info["MQ"] for record in record_dict[var_key])
                merged_record.info["MQ"] = complex_MQ
            except KeyError:
                continue

            try:
                complex_VARBIAS_forw = sum(int(record.info["VARBIAS"].split(":")[0]) for record in record_dict[var_key])
                complex_VARBIAS_revs = sum(int(record.info["VARBIAS"].split(":")[1]) for record in record_dict[var_key])
                merged_record.info["VARBIAS"] = "{}:{}".format(complex_VARBIAS_forw, complex_VARBIAS_revs)
            except KeyError:
                continue

            # Parse samples.items to get/update format field
            for item in merged_record.samples.items():
                item[1]["AF"] = complex_AF
                if "complex_VD" in locals():
                    try:
                        item[1]["VD"] = complex_VD
                        AD_ref = item[1]["AD"][0]
                        item[1]["AD"] = (AD_ref, complex_VD)
                    except KeyError:
                        continue

                if "complex_VARBIAS_forw" in locals() and "complex_VARBIAS_revs" in locals():
                    try:
                        item[1]["ALD"] = (complex_VARBIAS_forw, complex_VARBIAS_revs)
                    except KeyError:
                        continue

        record_dict[var_key] = [merged_record]

    return record_dict, vcf.header, complex_positions_list


def writeVCFOut(vcf_out_path, vcf_as_dict, input_header, complex_pos_list, method):
    '''
    Write a dictionary with pysam.VariantRecords to a vcf-file.

    param vcf_out_path: Path to the vcf-file to be created.
    param vcf_as_dict: A dictionary with all variants following the format,
    {"<chr>-<pos>-<REF allele>-<ALT allele>": [<pysam.VariantRecord object>]},
    where key is based on the variant and the value is a list with only one element/record.
    param input_header: pysam.VariantHeader object
    param complex_pos_list: A list with all variants that has more than one record in the vcf-file,
    ["<chr>-<pos>-<REF allele>-<ALT allele>", "<chr>-<pos>-<REF allele>-<ALT allele>", ...]
    param method: str method to use, max or sum, for calculation of AF
    return: None
    '''
    output_vcf = pysam.VariantFile(vcf_out_path, 'w', header=input_header)
    for key, value in vcf_as_dict.items():
        if key in complex_pos_list and method == "sum":
            record = value[0]
            new_record = output_vcf.new_record(contig=record.contig,
                                               start=(record.pos-1),
                                               stop=record.pos,
                                               )
            new_record.id = record.id
            new_record.ref = record.ref
            new_record.alts = record.alts

            # Set info
            new_record.info["AF"] = record.info["AF"]
            new_record.info["TYPE"] = record.info["TYPE"]
            new_record.info["COMPLEXAF"] = record.info["COMPLEXAF"]
            info_keys = ["SAMPLE", "DP", "VD", "BIAS", "REFBIAS", "VARBIAS", "QUAL", "MQ", "HIAF", "HICNT"]
            for key in info_keys:
                try:
                    new_record.info[key] = record.info[key]
                except KeyError:
                    continue

            # Set format
            format_keys = ["GT", "DP", "AF", "VD", "AD", "ALD", "RD"]
            for key in format_keys:
                try:
                    new_record.samples.items()[0][1][key] = record.samples.items()[0][1][key]
                except KeyError:
                    continue

            output_vcf.write(new_record)
        else:
            output_vcf.write(value[0])
    output_vcf.close()
    return


def main(vcf_in, vcf_out, method):
    if method == "max" or method == "sum":
        logging.info(f"The method for calculating allele frequencies for complex variants was set to {method}.")
        logging.info(f"Read {vcf_in}")
        vcf = pysam.VariantFile(vcf_in)
        vcf_as_dict, header, complex_positions = merge_records_complex_positions(vcf, method)
        writeVCFOut(vcf_out, vcf_as_dict, header, complex_positions, method)
        logging.info(f"Output will be written to {vcf_out}", )
    elif method == "skip":
        logging.info(f"The method for calculating allele frequencies for complex variants was set to {method}."
                     f"This step will be skipped!")
        logging.info(f"Output will be written to {vcf_out}")
        shutil.copyfile(vcf_in, vcf_out)
    else:
        raise ValueError(f"Invalid input. The method for calculating allele frequencies must be max, sum or skip."
                         f"The method given by user was: {method}")


if __name__ == "__main__":
    main(vcf_in=snakemake.input.vcf, vcf_out=snakemake.output.vcf, method=snakemake.params.merge_method)
