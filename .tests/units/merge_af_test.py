import sys
import os
import unittest
import pysam

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from merge_af import merge_records_complex_positions, writeVCFOut


class TestMergeAf(unittest.TestCase):

    def test_merge_records_complex_positions_max(self):
        # Use the max allele frequency of all observed af:s from the same position. 
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        vcf_as_dict,header = merge_records_complex_positions(vcf, "max")
        variant_key = "chr21-36432243-A-G"
        variant_record = vcf_as_dict[variant_key][0]
        self.assertEqual(variant_record.info["COMPLEXAF"], "max")
        self.assertAlmostEqual(variant_record.info["AF"][0], 0.0355)


    def test_merge_records_complex_positions_sum(self):
        # Use the sum of all allele frequencies of all observed af:s from the same position. 
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        vcf_as_dict,header = merge_records_complex_positions(vcf, "sum")
        variant_key = "chr21-36432243-A-G"
        variant_record = vcf_as_dict[variant_key][0]

        self.assertEqual(variant_record.info["COMPLEXAF"], "sum")

        expected_format= {'GT': (0 ,1), 
                          'DP': 282, 
                          'VD': 13, 
                          'AD': (268, 13), 
                          'AF': (0.0461,), 
                          'RD': (182, 86),
                          'ALD': (5, 8)}
        
        for item in variant_record.samples.items():
            for key,val in item[1].items():
                print(key, val)
                if key == 'AF':
                    self.assertAlmostEqual(expected_format[key][0], val[0])
                else:
                    self.assertEqual(expected_format[key], val)

    
    def test_writeVCFOut(self):
        vcfFile = ".tests/units/.tests/merge_af.normalized.sorted.vcf"
        vcf = pysam.VariantFile(vcfFile)
        # Create content for output file based on test file
        vcf_as_dict,header = merge_records_complex_positions(vcf, "sum")
        vcf_out_path = ".tests/units/.tests/merge_af.normalized.sorted.merged.af.vcf"
        # Create and write to output file
        writeVCFOut(vcf_out_path, vcf_as_dict, header)

        self.assertTrue(os.path.exists(vcf_out_path))

        vcf_out_file = pysam.VariantFile(vcf_out_path)
        vcf_out_header = vcf_out_file.header
        test_vcf_1 = "/home/monika/Projects/vcf-test/af_merge.header1.vcf"
        pysam.VariantFile(test_vcf_1, 'w', header=header)
        test_vcf_2 = "/home/monika/Projects/vcf-test/af_merge.header2.vcf"
        pysam.VariantFile(test_vcf_2, 'w', header=vcf_out_header)
 

        self.assertEqual(open(test_vcf_1).read(), open(test_vcf_2).read())
        # Cleanup
        os.remove(test_vcf_1)
        os.remove(test_vcf_2)
        os.remove(vcf_out_path)
