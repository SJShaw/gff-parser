#!/usr/bin/env python

"""
Parser module for gff3 files (dev)
"""

from collections import defaultdict
from typing import Dict, List, Optional, Union
import urllib  # TODO: python2/3 compatibility
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature, BeforePosition, AfterPosition, ExactPosition
from Bio.SeqRecord import SeqRecord


class Feature:
    """main feature object for gff class """
    def __init__(self, location: FeatureLocation, gff_type: str, attributes: Dict[str, str], parent: Optional["Feature"] = None) -> None:
        if parent is not None:
            assert isinstance(parent, Feature), type(parent)
        assert isinstance(location, FeatureLocation), type(location)
        self.name = attributes.pop("Name", None)
        self.attributes = attributes
        self.location = location
        self.gff_type = str(gff_type)
        self.gff_id = attributes.pop("ID", None)
        self.generics = []  # type: List["Feature"]
        self.exons = []  # type: List["Exon"]
        self.parent = None  # type: Feature

    def add_generic(self, feature: "Feature") -> None:
        assert self.gff_id, "parents must have IDs"
        assert isinstance(feature, Feature), type(feature)
        self.generics.append(feature)
        feature.parent = self

    def __repr__(self) -> str:
        return "%s(%s, %s)" % (self.gff_type, self.gff_id, self.location)

    def add_exon(self, exon: "Exon") -> None:
        """Adds exons to the feature"""
        assert isinstance(exon, Exon), type(exon)
        self.exons.append(exon)

    @property
    def locus_tag(self):
        if "locus_tag" in self.attributes:
            return self.attributes["locus_tag"]
        if self.parent:
            return self.parent.locus_tag
        return None


class Gene(Feature):
    """Gene gff feature; Can contain cds and mRNA."""
    def __init__(self, location: FeatureLocation, gff_type: str, attributes: Dict[str, str]) -> None:
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []  # type: List[CDS]
        self.mrnas = []  # type: List[RNA]
        self.other_rna = []  # type: List[RNA]

    def add_rna(self, rna: "RNA") -> None:
        "Adds RNA to parent object"
        assert isinstance(rna, RNA)
        if rna.gff_type == "mRNA":
            self.mrnas.append(rna)
        else:
            self.other_rna.append(rna)

    def add_cds(self, cds: "CDS") -> None:
        "Adds CDS children to the gene"
        assert isinstance(cds, CDS)
        self.cdss.append(cds)


class RNA(Feature):
    """ mRNA, rRNA, or tRNA features. All must have a Gene as a parent.
        mRNA may have CDS and Exon children.
    """
    def __init__(self, location: FeatureLocation, gff_type: str, parent: Gene, attributes: Dict[str, str]) -> None:
        super().__init__(location, gff_type, attributes, parent=parent)
        self.name = self.attributes.pop("Name", None) or (parent.name if parent else None)
        self.cdss = []  # type: List[CDS]
        self.exons = []  # type: List[Exon]

        if parent:
            assert isinstance(parent, Gene), type(parent)

    def add_cds(self, cds: "CDS") -> None:
        """Adds CDS children to the RNA object"""
        assert self.gff_type == "mRNA", self.gff_type
        assert isinstance(cds, CDS)
        self.cdss.append(cds)


class Exon(Feature):
    """ Exon feature, may be standalone or a child of a Gene or mRNA feature.
        May have generic child features.
    """
    def __init__(self, location: FeatureLocation, parent: Union[Gene, RNA], attributes: Dict[str, str]) -> None:
        super().__init__(location, "exon", attributes)
        assert isinstance(parent, Feature), "Unknown parent type: %s" % parent
        self.parent = parent


class CDS(Feature):
    """ CDS feature, may be standalone or a child of a Gene or mRNA feature.
    """
    def __init__(self, location: FeatureLocation, parent: Union[Gene, RNA], attributes: Dict[str, str]) -> None:
        super().__init__(location, "CDS", attributes)
        if parent is not None:
            assert isinstance(parent, (RNA, Gene)), "Unknown parent type: %s" % parent
        assert location.strand in [-1, 1], "Invalid CDS strand"
        self.parent = parent


class GFFRecord:
    def __init__(self, name: str) -> None:
        self.name = name
        self.parent_features = {}  # type: Dict[str, Feature]
        self.all_features = []  # type: List[Feature]
        self.cds_features = []  # type: List[CDS]
        self.cds_features_by_parent = {}  # type: Dict[Feature, CDS]
        self.total_features = 0


def decode(text: str) -> str:
    return urllib.parse.unquote(text, errors="strict")


def build_attributes(section: str) -> Dict[str, str]:
    attributes = {}
    for attribute in section.split(";"):
        if not attribute:
            continue
        parts = attribute.split("=")
        assert len(parts) == 2, "Invalid attribute: %s" % attribute
        attributes[parts[0]] = decode(parts[1])
    return attributes


def interpret_strand(strand: str, strict: bool = False) -> int:
    possibilities = {
        "+": 1,
        "-": -1,
        ".": 0,
    }
    fallbacks = {
        "+1": 1,
        "1": 1,
        "-1": -1,
        "0": 0,
    }
    result = possibilities.get(strand)
    if result is None and not strict:
        result = fallbacks.get(strand)
    if result is None:
        raise ValueError("Unknown strand representation: %s" % strand)
    return result


def construct_feature(record: GFFRecord, line: str, parent_lines: Dict[str, str]) -> Feature:
    parts = line.strip().split(sep="\t" if "\t" in line else None, maxsplit=9)
    assert len(parts) == 9, line

    try:
        seqid, _, gff_type, start, end, _, strand_dir, _, attribute_string = parts
        start = ExactPosition(int(start) - 1)  # 0-indexed as FeatureLocation expects
        end = ExactPosition(int(end))
    except ValueError:
        print(line)
        raise

    assert seqid == record.name

    strand = interpret_strand(strand_dir)

    try:
        attributes = build_attributes(attribute_string)
    except:
        print("bad attributes line:", attribute_string, file=sys.stderr)
        raise

    if attributes.get("partial") == "true" and ("start_range" in attributes or "end_range" in attributes):
        attributes.pop("partial")
        start_range = attributes.pop("start_range", "%s,%s" % (start, end))
        end_range = attributes.pop("end_range", "%s,%s" % (start, end))
        if start_range.startswith("."):
            start = BeforePosition(int(start))
        if end_range.endswith("."):
            end = AfterPosition(int(end))
    location = FeatureLocation(start, end, strand)

    part = attributes.pop("part", None)
    if not attributes.get("ID") and part:
        if not attributes.get("Name"):
            raise ValueError("Feature specified as part but has no reference")
        attributes["ID"] = "%s-%s" % (gff_type, attributes["Name"])

    parent_name = attributes.pop("Parent", None)
    if parent_name is not None and parent_name not in parent_lines:
        raise ValueError("Parent feature referenced but is not known: %s" % parent_name)

    parent = record.parent_features.get(parent_name)
    if parent_name and not parent:
        parent = construct_feature(record, parent_lines[parent_name], parent_lines)

    feature = None

    if gff_type == "pseudogene":
        gff_type = "gene"
        attributes["pseudo"] = "true"

    if gff_type == "gene":
        if attributes.get("ID") in record.parent_features:
            gene = record.parent_features[attributes.get("ID")]
            gene.location = gene.location + location  # TODO obey the is_ordered=true attribute
            return gene
        else:
            gene = Gene(location, gff_type, attributes)
            feature = gene
    elif gff_type in ["mRNA", "rRNA", "tRNA", "CDS", "exon"]:
        if gff_type in ["mRNA", "tRNA", "rRNA"]:
            if parent is not None:
                assert isinstance(parent, Gene), str(parent)
            rna = RNA(location, gff_type, parent, attributes)
            feature = rna
            if parent:
                parent.add_rna(rna)

        elif gff_type == "exon":
            exon = Exon(location, parent, attributes)
            feature = exon
            parent.add_exon(exon)

        elif gff_type == "CDS":
            if parent:  # not all GFF files will specify genes at all
                assert isinstance(parent, (Gene, RNA)), str(parent)
            # merge if no ID and shared gene parent
            if "ID" not in attributes and parent in record.cds_features_by_parent:
                cds = record.cds_features_by_parent[parent]
                cds.location = cds.location + location
                # TODO: manage attribute merging/changes
                return cds
            cds = CDS(location, parent, attributes)
            record.cds_features_by_parent[parent] = cds
            if parent:
                parent.add_cds(cds)
            feature = cds
            record.cds_features.append(cds)
    else:
        generic = Feature(location, gff_type, attributes, parent=parent)
        feature = generic
        if parent:
            parent.add_generic(generic)
    assert feature is not None
    record.all_features.append(feature)
    # features that have an id must have a unique id
    if feature.gff_id:
        assert feature.gff_id not in record.parent_features, feature
        record.parent_features[feature.gff_id] = feature
    return feature


def populate_record(record: GFFRecord, named_lines: Dict[str, str], unnamed_lines: List[str]) -> None:
    # build named features first
    for name, line in named_lines.items():
        # skip if it was generated recursively as a dependency
        if name in record.parent_features:
            continue
        feature = construct_feature(record, line, named_lines)
        assert feature is not None, line
        assert feature.gff_id == name, "%s != %s" % (feature.gff_id, name)

    # then unnamed features
    for line in unnamed_lines:
        construct_feature(record, line, named_lines)


def parse_gff(filename: str, strict=False) -> Dict[str, Feature]:
    """
    Parses gff file
    """
    records_by_name = {}  # type: Dict[str, GFFRecord]
    named_by_record = defaultdict(dict)  # type: Dict[str, Dict[str, str]]
    anon_by_record = defaultdict(list)  # type: Dict[str, List[str]]

    with open(filename) as handle:
        for line in handle:
            # TODO: keep sequences, ensure annotations don't follow
            if line.startswith(">"):
                break
            # skip comments
            if line.startswith("#"):
                continue
            # skip empty lines
            if not line.strip():
                continue

            parts = line.strip().split(sep="\t" if "\t" in line else None, maxsplit=9)
            assert len(parts) == 9, line
            seqid = parts[0]
            # construct a record skeleton
            if seqid not in records_by_name:
                #  if len(seqid) > 16:  # TODO: doesn't work on refseq ids, need to change
                #       raise ValueError("Record name too long: %s" % seqid)
                records_by_name[seqid] = GFFRecord(seqid)

            # track feature lines
            try:
                attributes = build_attributes(parts[8])
            except:
                print("bad attributes line:", parts[8], file=sys.stderr)
                raise

            # handle multi-line features that should have IDs but can't
            part = attributes.pop("part", None)
            if not attributes.get("ID") and part:
                if not attributes.get("Name"):
                    raise ValueError("Feature specified as part but has no reference")
                gff_type = parts[2]
                attributes["ID"] = "%s-%s" % (gff_type, attributes["Name"])

            if "ID" in attributes:
                existing = named_by_record[seqid].get(attributes["ID"])
                if existing and existing != line:
                    raise ValueError("Duplicate ID detected in record %s: %s" % (seqid, attributes["ID"]))
                # exact duplicates overridden
                named_by_record[seqid][attributes["ID"]] = line
            else:
                anon_by_record[seqid].append(line)

    # populate records
    for record_name, record in records_by_name.items():
        populate_record(record, named_by_record[record_name], anon_by_record[record_name])

    print("parsed %d records" % len(records_by_name))
    return list(records_by_name.values())


QUALIFIER_MAPPING = {
    "Dbxref": "db_xref",
    "Note": "note",
}

ATTRIBUTES = {
}


TYPE_MAPPING = {
    "Gene": "gene",
    "Src": "source",
}


def convert_record_to_biopython(gff_record: GFFRecord) -> SeqRecord:
    seq_record = SeqRecord(seq=Seq("", generic_dna), id=gff_record.name, name=gff_record.name)

    old_to_new = {}

    for feature in gff_record.all_features:
        new_feature = SeqFeature(feature.location)
        new_feature.type = TYPE_MAPPING.get(feature.attributes.pop("gbkey", None) or feature.gff_type, feature.gff_type)
        assert new_feature.type

        for key, val in feature.attributes.items():
            if key in ATTRIBUTES:
                setattr(new_feature, ATTRIBUTES[key], val)
                continue
            qualifier_name = QUALIFIER_MAPPING.get(key, key)
            if qualifier_name in new_feature.qualifiers:
                new_feature.qualifiers[qualifier_name].append(val)
            else:
                new_feature.qualifiers[qualifier_name] = [val]
        seq_record.features.append(new_feature)

        if "locus_tag" not in new_feature.qualifiers:
            locus = feature.locus_tag
            if locus:
                new_feature.qualifiers["locus_tag"] = [locus]
        old_to_new[feature] = new_feature

    return seq_record


def convert_gff_to_genbank(gff_records: List[GFFRecord]) -> List[SeqRecord]:
    seq_records = []
    for gff_record in gff_records:
        seq_records.append(convert_record_to_biopython(gff_record))
    return seq_records


if __name__ == "__main__":
    seq_records = convert_gff_to_genbank(parse_gff(sys.argv[1]))
#    from helperlibs.bio import seqio
#    with open("gff_to.gbk", "w") as handle:
#        seqio.write(seq_records, handle, "genbank")
