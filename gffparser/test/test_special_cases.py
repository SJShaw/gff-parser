# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from Bio.SeqFeature import CompoundLocation

from gffparser import gff_parser


def get_path_of_file(filename):
    return os.path.join("gffparser", "test", "data", filename)


def records_from_local_file(filename):
    return gff_parser.parse_gff(get_path_of_file(filename))


class TestPartGenes(unittest.TestCase):
    def test_NCBI(self):
        record = records_from_local_file("FNCA01000004.1.fragment.gff3")[0]
        non_part = record.parent_features["gene-SAMN04488589_1459"]
        assert str(non_part.location) == "[130580:131342](-)"

        part = record.parent_features["gene-SAMN04488589_1460"]
        assert str(part.location) == "join{[131545:131582](-), [131701:131740](-)}"


class TestDuplicates(unittest.TestCase):
    def test_exact_duplicates(self):
        filename = "duplicate_lines.gff3"
        record = records_from_local_file(filename)[0]
        with open(get_path_of_file(filename)) as handle:
            lines = handle.read().splitlines()
        assert len(record.all_features) == 2
        assert len(lines) == 4

    def test_cds_merging(self):
        record = records_from_local_file("cds_sharing_ids.gff3")[0]
        assert len(record.cds_features) == 1
        assert isinstance(record.cds_features[0].location, CompoundLocation)


class TestAttributes(unittest.TestCase):
    def build(self, text):
        return gff_parser.build_attributes(text)

    def test_repeated_separator(self):
        attributes = self.build("a=long thing;;B=7")
        assert attributes == {"a": "long thing", "B": "7"}

    def test_extra_separator(self):
        attributes = self.build(";a=long thing;B=7")
        assert attributes == {"a": "long thing", "B": "7"}

        attributes = self.build("a=long thing;B=7;")
        assert attributes == {"a": "long thing", "B": "7"}

    def test_quoted_separators(self):
        attributes = self.build('a="two; sections";b=one " section')
        assert attributes == {"a": "two; sections", "b": 'one " section'}

        attributes = self.build('a="two; sections";b=one " section')
        assert attributes == {"a": "two; sections", "b": 'one " section'}

        attributes = self.build('weird "key=normal value')
        assert attributes == {'weird "key': "normal value"}

        with self.assertRaisesRegex(ValueError, "unclosed quote contains an attribute separator"):
            self.build('bad="missing; close')
        with self.assertRaisesRegex(ValueError, "unclosed quote contains an attribute separator"):
            self.build('a="missing close;b=good value')
