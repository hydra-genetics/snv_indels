#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Monika Brandt"
__copyright__ = "Copyright 2024, Monika Brandt"
__email__ = "monika.brandt@scilifelab.uu.se"
__license__ = "GPL-3"

import sys
import os
import unittest
import pysam

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_af import merge_records_complex_positions, writeVCFOut  # noqa


class TestMergeAf(unittest.TestCase):

    def test_merge_records_complex_positions_max(self):
        # Use the max allele frequency of all observed af:s from the same position.
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        vcf_as_dict, header, complex_list = merge_records_complex_positions(vcf, "max")
        variant_key = "chr21-36432243-A-G"
        variant_record = vcf_as_dict[variant_key][0]
        self.assertEqual(variant_record.info["COMPLEXAF"], "max")
        self.assertEqual(variant_record.info["TYPE"], "Complex")
        self.assertAlmostEqual(variant_record.info["AF"][0], 0.0355)

    def test_merge_records_complex_positions_sum(self):
        # Use the sum of all allele frequencies of all observed af:s from the same position.
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        vcf_as_dict, header, complex_list = merge_records_complex_positions(vcf, "sum")
        variant_key = "chr21-36432243-A-G"
        variant_record = vcf_as_dict[variant_key][0]

        self.assertEqual(variant_record.info["COMPLEXAF"], "sum")
        self.assertEqual(variant_record.info["TYPE"], "Complex")
        self.assertAlmostEqual(variant_record.info["AF"][0], 0.0461)

        expected_format = {'GT': (0, 1),
                           'DP': 282,
                           'VD': 13,
                           'AD': (268, 13),
                           'AF': (0.0461,),
                           'RD': (182, 86),
                           'ALD': (5, 8)}

        for item in variant_record.samples.items():
            for key, val in item[1].items():
                if key == 'AF':
                    self.assertAlmostEqual(expected_format[key][0], val[0])
                else:
                    self.assertEqual(expected_format[key], val)

    def test_writeVCFOut_header(self):
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        # Create content for output file based on test file
        vcf_as_dict, header, complex_list = merge_records_complex_positions(vcf, "sum")
        vcf_out_path = ".tests/units/.tests/merge_af.normalized.sorted.merged.af.vcf"
        # Create and write to output file
        writeVCFOut(vcf_out_path, vcf_as_dict, header, complex_list, "sum")

        # Make sure that writeVCFOut() creates an output vcf.
        self.assertTrue(os.path.exists(vcf_out_path))

        # Get header from vcf-file created and written by writeVCFOut().
        vcf_out_file = pysam.VariantFile(vcf_out_path)
        vcf_out_header = vcf_out_file.header

        # Write header parsed by writeVCFOut() to a test-file
        test_vcf_1 = ".tests/units/.tests/af_merge.header1.vcf"
        pysam.VariantFile(test_vcf_1, 'w', header=vcf_out_header)

        # Write header NOT parsed by writeVCFOut() to a test-file
        test_vcf_2 = ".tests/units/.tests/af_merge.header2.vcf"
        pysam.VariantFile(test_vcf_2, 'w', header=header)

        self.assertEqual(open(test_vcf_1).read(), open(test_vcf_2).read())

        # Cleanup
        os.remove(test_vcf_1)
        os.remove(test_vcf_2)
        os.remove(vcf_out_path)

    def test_info_and_format_keys(self):
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        # Create content for output file based on test file
        vcf_as_dict, header, complex_list = merge_records_complex_positions(vcf, "sum")
        vcf_out_path = ".tests/units/.tests/merge_af.normalized.sorted.merged.af.vcf"
        # Create and write to output file
        writeVCFOut(vcf_out_path, vcf_as_dict, header, complex_list, "sum")

        # Make sure that writeVCFOut() creates an output vcf.
        self.assertTrue(os.path.exists(vcf_out_path))

        # Variant to select for testing
        variant_key = "chr21-36432243-A-G"

        # Get record from output file.
        vcf_out_file = pysam.VariantFile(vcf_out_path)

        expected_info = ['SAMPLE',
                         'TYPE',
                         'DP',
                         'VD',
                         'AF',
                         'BIAS',
                         'REFBIAS',
                         'VARBIAS',
                         'QUAL',
                         'MQ',
                         'HIAF',
                         'HICNT',
                         'COMPLEXAF']

        expected_format = {'GT': (0, 1),
                           'DP': 282,
                           'VD': 13,
                           'AD': (268, 13),
                           'AF': (0.0461,),
                           'RD': (182, 86),
                           'ALD': (5, 8)}

        # Make sure that only keys updated by the code is given in the output.
        for record in vcf_out_file.fetch():
            var_key_record = "{}-{}-{}-{}".format(record.chrom, record.pos, record.ref, record.alts[0])
            if var_key_record == variant_key:
                record_info = record.info.keys()
                self.assertEqual(expected_info.sort(), record_info.sort())
                self.assertEqual(record.samples.items()[0][1].keys().sort(), list(expected_format.keys()).sort())

        # Cleanup
        os.remove(vcf_out_path)

    def test_missing_key(self):
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.missing_VD.vcf"
        vcf = pysam.VariantFile(vcfFile)
        # Create content for output file based on test file
        method = "sum"
        vcf_as_dict, input_header, complex_pos_list = merge_records_complex_positions(vcf, method)
        vcf_out_path = ".tests/units/.tests/merge_af.normalized.sorted.missing_VD.merged.af.vcf"
        # Create and write to output file
        writeVCFOut(vcf_out_path, vcf_as_dict, input_header, complex_pos_list, method)

        # Make sure that writeVCFOut() creates an output vcf.
        self.assertTrue(os.path.exists(vcf_out_path))

        # Variant to select for testing
        variant_key = "chr21-36432243-A-G"

        # Get record from output file.
        vcf_out_file = pysam.VariantFile(vcf_out_path)

        expected_info = ['SAMPLE',
                         'TYPE',
                         'DP',
                         'AF',
                         'BIAS',
                         'REFBIAS',
                         'VARBIAS',
                         'QUAL',
                         'MQ',
                         'HIAF',
                         'HICNT',
                         'COMPLEXAF']

        expected_format = {'GT': (0, 1),
                           'DP': 282,
                           'AD': (268, 13),
                           'AF': (0.0461,),
                           'RD': (182, 86),
                           'ALD': (5, 8)}

        # Make sure that only keys updated by the code is given in the output.
        for record in vcf_out_file.fetch():
            var_key_record = "{}-{}-{}-{}".format(record.chrom, record.pos, record.ref, record.alts[0])
            if var_key_record == variant_key:
                record_info = record.info.keys()
                self.assertEqual(expected_info.sort(), record_info.sort())
                self.assertEqual(record.samples.items()[0][1].keys().sort(), list(expected_format.keys()).sort())

        # Cleanup
        os.remove(vcf_out_path)
