# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from gffparser import parse_gff
from gffparser import gff_parser


class TestPartGenes(unittest.TestCase):
    @unittest.expectedFailure
    def test_NCBI(self):
        filepath = os.path.join("gffparser", "test", "data", "FNCA01000004.1.fragment.gff3")
        record = parse_gff(filepath)[0]
        non_part = record.parent_features["gene-SAMN04488589_1459"]
        assert str(non_part.location) == ""

        part = record.parent_features["gene-SAMN04488589_1460"]
        assert str(part.location) == ""


class TestDuplicates(unittest.TestCase):
    def test_exact_duplicates(self):
        filepath = os.path.join("gffparser", "test", "data", "duplicate_lines.gff3")
        record = parse_gff(filepath)[0]
        with open(filepath) as handle:
            lines = handle.read().splitlines()
        assert len(record.all_features) == 2
        assert len(lines) == 4


class TestAttributes(unittest.TestCase):
    def test_repeated_separator(self):
        attribute_section = "a=long thing;;B=7"
        attributes = gff_parser.build_attributes(attribute_section)
        assert attributes == {"a": "long thing", "B": "7"}

    def test_extra_separator(self):
        attribute_section = ";a=long thing;B=7"
        attributes = gff_parser.build_attributes(attribute_section)
        assert attributes == {"a": "long thing", "B": "7"}

        attribute_section = "a=long thing;B=7;"
        attributes = gff_parser.build_attributes(attribute_section)
        assert attributes == {"a": "long thing", "B": "7"}
