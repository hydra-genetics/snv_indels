import pysam
import statistics
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

    for pos in complex_positions_list:
        # if method "max" is given by the user the code will return the record with the highest allele frequency
        # as is without any updates on for instance AF, VF or AD.
        # A new info key will be added: COMPLEXAF=max
        if method == "max":
            complex_AF = max(record.info["AF"][0] for record in record_dict[pos])
            merged_record = record_dict[pos][[record.info["AF"][0] == complex_AF for record in record_dict[pos]].index(True)]
            merged_record.info["COMPLEXAF"] = method
            merged_record.info["TYPE"] = "Complex"

        # if method "sum" is given by the user, return one record for the given variant where
        # AF is the sum of AF for all records for the variant.
        # A new info key will be added: COMPLEXAF=sum
        if method == "sum":
            complex_AF = sum(record.info["AF"][0] for record in record_dict[pos])
            complex_VD = sum(record.info["VD"] for record in record_dict[pos])
            complex_HIAF = sum(record.info["HIAF"] for record in record_dict[pos])
            complex_HICNT = sum(record.info["HICNT"] for record in record_dict[pos])
            complex_QUAL = statistics.mean(record.info["QUAL"] for record in record_dict[pos])
            complex_MQ = statistics.mean(record.info["MQ"] for record in record_dict[pos])
            complex_VARBIAS_forw = sum(int(record.info["VARBIAS"].split(":")[0]) for record in record_dict[pos])
            complex_VARBIAS_revs = sum(int(record.info["VARBIAS"].split(":")[1]) for record in record_dict[pos])
            # Pick first record in list to use as a template for updates.
            merged_record = record_dict[pos][0]
            merged_record.info["COMPLEXAF"] = method
            merged_record.info["TYPE"] = "Complex"
            merged_record.info["AF"] = complex_AF
            merged_record.info["VD"] = complex_VD
            merged_record.info["HIAF"] = complex_HIAF
            merged_record.info["HICNT"] = complex_HICNT
            merged_record.info["QUAL"] = complex_QUAL
            merged_record.info["MQ"] = complex_MQ
            merged_record.info["VARBIAS"] = "{}:{}".format(complex_VARBIAS_forw, complex_VARBIAS_revs)
            # Parse samples.items to get/update format field
            for item in merged_record.samples.items():
                item[1]["AF"] = complex_AF
                item[1]["VD"] = complex_VD
                AD_ref = item[1]["AD"][0]
                item[1]["AD"] = (AD_ref, complex_VD)
                item[1]["ALD"] = (complex_VARBIAS_forw, complex_VARBIAS_revs)

        record_dict[pos] = [merged_record]

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
            new_record.info["SAMPLE"] = record.info["SAMPLE"]
            new_record.info["TYPE"] = record.info["TYPE"]
            new_record.info["DP"] = record.info["DP"]
            new_record.info["VD"] = record.info["VD"]
            new_record.info["AF"] = record.info["AF"]
            new_record.info["BIAS"] = record.info["BIAS"]
            new_record.info["REFBIAS"] = record.info["REFBIAS"]
            new_record.info["VARBIAS"] = record.info["VARBIAS"]
            new_record.info["QUAL"] = record.info["QUAL"]
            new_record.info["MQ"] = record.info["MQ"]
            new_record.info["HIAF"] = record.info["HIAF"]
            new_record.info["HICNT"] = record.info["HICNT"]
            new_record.info["COMPLEXAF"] = record.info["COMPLEXAF"]

            # Set format
            new_record.samples.items()[0][1]["GT"] = record.samples.items()[0][1]["GT"]
            new_record.samples.items()[0][1]["DP"] = record.samples.items()[0][1]["DP"]
            new_record.samples.items()[0][1]["AF"] = record.samples.items()[0][1]["AF"]
            new_record.samples.items()[0][1]["VD"] = record.samples.items()[0][1]["VD"]
            new_record.samples.items()[0][1]["AD"] = record.samples.items()[0][1]["AD"]
            new_record.samples.items()[0][1]["ALD"] = record.samples.items()[0][1]["ALD"]
            new_record.samples.items()[0][1]["RD"] = record.samples.items()[0][1]["RD"]

            output_vcf.write(new_record)
        else:
            output_vcf.write(value[0])
    output_vcf.close()
    return


if __name__ == "__main__":
    vcfFile = snakemake.input.vcf
    merge_method = snakemake.params.merge_method

    if merge_method == "max" or merge_method == "sum":
        vcf = pysam.VariantFile(vcfFile)
        vcf_as_dict, header, complex_positions = merge_records_complex_positions(vcf, merge_method)
        vcf_out_path = snakemake.output.vcf
        writeVCFOut(vcf_out_path, vcf_as_dict, header, complex_positions, merge_method)
    else:
        merge_method = "skip"
        snakemake.output.vcf = vcfFile
