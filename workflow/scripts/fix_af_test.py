#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


import io
import os
import pysam
import unittest
from dataclasses import dataclass
from fix_af import getCaller, modifyHeader, fixFreebayes, writeNewVcf


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
                    input="snv_indels/mutect2/NA12878_T.vcf",
                    expected="mutect2",
                ),
                TestCase(
                    name="Input path too short",
                    input="mutect2/NA12878_T.vcf",
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
                    expected=pysam.VariantFile("workflow/scripts/.tests/fix_af.modifyHeader.pisces_mutect2.expected.vcf").header,
                ),
                TestCase(
                    name="Caller is mutect2",
                    caller="mutect2",
                    header=pysam.libcbcf.VariantHeader(),
                    expected=pysam.VariantFile("workflow/scripts/.tests/fix_af.modifyHeader.pisces_mutect2.expected.vcf").header,
                ),
                TestCase(
                    name="Caller is varscan",
                    caller="varscan",
                    header=pysam.libcbcf.VariantHeader(),
                    expected=pysam.VariantFile("workflow/scripts/.tests/fix_af.modifyHeader.varscan.expected.vcf").header,
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
                    input=pysam.VariantFile("workflow/scripts/.tests/fix_af.fixFreebayes.input.vcf"),
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
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.freebayes.vcf.vcf"),
                    caller="freebayes",
                    expected="workflow/scripts/.tests/fix_af.writeNewVcf.freebayes.expected.vcf",
                ),
                TestCase(
                    name="Caller is mutect2",
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.mutect2_vardict.vcf.vcf"),
                    caller="mutect2",
                    expected="workflow/scripts/.tests/fix_af.writeNewVcf.mutect2_vardict.expected.vcf",
                ),
                TestCase(
                    name="Caller is pisces",
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.pisces.vcf.vcf"),
                    caller="pisces",
                    expected="workflow/scripts/.tests/fix_af.writeNewVcf.pisces.expected.vcf",
                ),
                TestCase(
                    name="Caller is vardict",
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.mutect2_vardict.vcf.vcf"),
                    caller="vardict",
                    expected="workflow/scripts/.tests/fix_af.writeNewVcf.mutect2_vardict.expected.vcf",
                ),
                TestCase(
                    name="Caller is varscan",
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.varscan.vcf.vcf"),
                    caller="varscan",
                    expected="workflow/scripts/.tests/fix_af.writeNewVcf.varscan.expected.vcf",
                ),
                TestCase(
                    name="Caller not available",
                    path="workflow/scripts/.tests/fix_af.writeNewVcf.actual.vcf",
                    vcf=pysam.VariantFile("workflow/scripts/.tests/fix_af.writeNewVcf.freebayes.vcf.vcf"),
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
