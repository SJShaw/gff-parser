# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from gffparser import parse_gff

class TestPartGenes(unittest.TestCase):
    def test_NCBI(self):
        filepath = os.path.join("gffparser", "test", "data", "FNCA01000004.1.fragment.gff3")
        print(os.getcwd() + "/" + filepath)
        record = parse_gff(filepath)
        non_part = record.parent_features["gene-SAMN04488589_1459"]
        assert str(non_part.location) == ""

        part = record.parent_features["gene-SAMN04488589_1460"]
        assert str(part.location) == ""
