# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from Bio.SeqFeature import FeatureLocation, CompoundLocation

from gffparser import parse_gff


class TestParentLinking(unittest.TestCase):
    def test_simple(self):
        filepath = os.path.join("gffparser", "test", "data", "simple_parent.gff3")
        record = parse_gff(filepath)[0]
        gene = record.parent_features["gene"]
        cds = record.cds_features[0]
        assert cds.parent == gene

    def test_wrong_ordering(self):
        filepath = os.path.join("gffparser", "test", "data", "out_of_order_parent.gff3")
        record = parse_gff(filepath)[0]
        gene = record.parent_features["gene"]
        cds = record.cds_features[0]
        assert cds.parent == gene

    def test_RNA_gene(self):
        filepath = os.path.join("gffparser", "test", "data", "RNA_gene.gff3")
        record = parse_gff(filepath)[0]
        gene = record.parent_features["A"]
        assert len(gene.generics) == 1
        assert gene.gff_type == "RNA_gene"
        rna = record.parent_features["B"]
        assert rna.parent == gene

class TestCDSMerging(unittest.TestCase):
    def test_distinct(self):
        filepath = os.path.join("gffparser", "test", "data", "distinct_cds_same_parent.gff3")
        record = parse_gff(filepath)[0]
        assert len(record.cds_features) == 2
        assert record.cds_features[0].gff_id != record.cds_features[1].gff_id
        assert isinstance(record.cds_features[0].location, FeatureLocation)
        assert isinstance(record.cds_features[1].location, FeatureLocation)

    def test_indistinct(self):
        filepath = os.path.join("gffparser", "test", "data", "indistinct_cds_same_parent.gff3")
        record = parse_gff(filepath)[0]
        assert len(record.cds_features) == 1
        assert isinstance(record.cds_features[0].location, CompoundLocation)
