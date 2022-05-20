#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


import io
import os
import pysam
import sys
import unittest
from dataclasses import dataclass

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.abspath(os.path.join(TEST_DIR, "../../workflow/scripts"))
sys.path.insert(0, SCRIPT_DIR)

from fix_af import getCaller, modifyHeader, fixFreebayes, writeNewVcf  # noqa


class TestGetCaller(unittest.TestCase):
    def test_getcaller(self):
        @dataclass
        class TestCase:
            name: str
            input: str
            expected: str

        testcases = [
                TestCase(
                    name="Identify caller successfully",
                    input="snv_indels/gatk_mutect2/NA12878_T.vcf",
                    expected="gatk_mutect2",
                ),
                TestCase(
                    name="Input path too short",
                    input="gatk_mutect2/NA12878_T.vcf",
                    expected="ValueError",
                ),
        ]

        for case in testcases:
            try:
                actual = getCaller(case.input)
                self.assertEqual(
                    case.expected,
                    actual,
                    "failed test '{}': expected {}, got {}".format(
                        case.name, case.expected, actual
                    ),
                )
            except ValueError:
                if case.expected == "ValueError":
                    assert True
                else:
                    assert False


class TestModifyHeader(unittest.TestCase):
    def test_modifyheader(self):
        @dataclass
        class TestCase:
            name: str
            caller: str
            header: pysam.libcbcf.VariantHeader
            expected: str

        testcases = [
                TestCase(
                    name="Caller is pisces",
                    caller="pisces",
                    header=pysam.libcbcf.VariantHeader(),
                    expected=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.modifyHeader.pisces_mutect2.expected.vcf"
                    ).header,
                ),
                TestCase(
                    name="Caller is gatk_mutect2",
                    caller="gatk_mutect2",
                    header=pysam.libcbcf.VariantHeader(),
                    expected=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.modifyHeader.pisces_mutect2.expected.vcf"
                    ).header,
                ),
                TestCase(
                    name="Caller is varscan",
                    caller="varscan",
                    header=pysam.libcbcf.VariantHeader(),
                    expected=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.modifyHeader.varscan.expected.vcf"
                    ).header,
                ),
        ]

        for case in testcases:
            actual = modifyHeader(case.caller, case.header)
            self.assertEqual(
                str(case.expected),
                str(actual),
                "failed test '{}': expected {}, got {}".format(
                    case.name, case.expected, actual
                ),
            )


class TestFixFreebayes(unittest.TestCase):
    def test_fixfreebayes(self):
        @dataclass
        class TestCase:
            name: str
            input: pysam.libcbcf.VariantFile
            expected: tuple

        testcases = [
                TestCase(
                    name="Fix ad field in freebayes vcf",
                    input=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.fixFreebayes.input.vcf"
                    ),
                    expected=(0.2,),
                ),
        ]

        for case in testcases:
            for row in case.input.fetch():
                actual = fixFreebayes(case.input.header, row)
            self.assertTupleEqual(
                case.expected,
                actual,
                "failed test '{}': expected {}, got {}".format(
                    case.name, case.expected, actual
                ),
            )


class TestWriteNewVcf(unittest.TestCase):
    def test_modifyheader(self):
        @dataclass
        class TestCase:
            name: str
            path: str
            vcf: pysam.libcbcf.VariantFile
            caller: str
            expected: str

        testcases = [
                TestCase(
                    name="Caller is freebayes",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.freebayes.vcf.vcf"
                    ),
                    caller="freebayes",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.freebayes.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller is haplotypecaller",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.haplotypecaller.vcf.vcf"
                    ),
                    caller="haplotypecaller",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.haplotypecaller.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller is gatk_mutect2",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.mutect2_vardict.vcf.vcf"
                    ),
                    caller="gatk_mutect2",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.mutect2_vardict.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller is pisces",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.pisces.vcf.vcf"
                    ),
                    caller="pisces",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.pisces.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller is vardict",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.mutect2_vardict.vcf.vcf"
                    ),
                    caller="vardict",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.mutect2_vardict.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller is varscan",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.varscan.vcf.vcf"
                    ),
                    caller="varscan",
                    expected=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.varscan.expected.vcf"
                    ),
                ),
                TestCase(
                    name="Caller not available",
                    path=(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.actual.vcf"
                    ),
                    vcf=pysam.VariantFile(
                        ".tests/units/.tests/"
                        "fix_af.writeNewVcf.freebayes.vcf.vcf"
                    ),
                    caller="snver",
                    expected="ValueError",
                ),
        ]

        for case in testcases:
            try:
                writeNewVcf(case.path, case.vcf.header, case.vcf, case.caller)
                self.assertListEqual(
                    list(io.open(case.expected)),
                    list(io.open(case.path)),
                    "failed test '{}': files are different".format(
                        case.name
                    ),
                )
                os.remove(case.path)
            except ValueError:
                if case.expected == "ValueError":
                    assert True
                else:
                    assert False
