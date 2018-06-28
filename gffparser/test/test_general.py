# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from gffparser import gff_parser


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
            with self.assertRaisesRegex(ValueError, "Unknown strand"):
                gff_parser.interpret_strand(strand, strict=True)

    def test_other(self):
        for strict in [True, False]:
            for strand in [-1, 'x', None]:
                with self.assertRaisesRegex(ValueError, "Unknown strand"):
                    gff_parser.interpret_strand(strand, strict=strict)
