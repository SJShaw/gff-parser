# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from gffparser import gff_parser

from .test_special_cases import records_from_local_file


class TestStrand(unittest.TestCase):
    TO_SPEC = ["+", ".", "-"]
    KNOWN_EXTRA = ["+1", "1", "0", "-1"]

    def test_spec(self):
        for strict in [True, False]:
            for strand in TestStrand.TO_SPEC:
                assert gff_parser.interpret_strand(strand, strict=strict) in {-1, 0, 1}

    def test_extra(self):
        for strand in TestStrand.KNOWN_EXTRA:
            assert gff_parser.interpret_strand(strand, strict=False) in {-1, 0, 1}
            with self.assertRaisesRegex(gff_parser.GFFParseError, "Unknown strand"):
                gff_parser.interpret_strand(strand, strict=True)

    def test_other(self):
        for strict in [True, False]:
            for strand in [-1, 'x', None]:
                with self.assertRaisesRegex(gff_parser.GFFParseError, "Unknown strand"):
                    gff_parser.interpret_strand(strand, strict=strict)


class TestMultiRecord(unittest.TestCase):
    def test_parsing(self):
        records = records_from_local_file("multi_record.gff3")
        records = sorted(records, key=lambda x: x.name)
        assert len(records) == 2
        assert records[0].name == "seq1"
        assert records[1].name == "seq2"

    def test_conversion(self):
        records = records_from_local_file("multi_record.gff3")
        assert len(records) == 2

        biopython_recs = gff_parser.convert_gff_to_biopython(records)
        assert len(biopython_recs) == 2
        assert [rec.name for rec in records] == [bio.id for bio in biopython_recs]
