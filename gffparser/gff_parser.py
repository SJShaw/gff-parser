#!/usr/bin/env python

""" Parser for GFF3 files, also has the capability to convert a parsed file to
    Genbank via a BioPython SeqRecord
"""

from collections import defaultdict
from typing import Dict, Iterator, List, Optional, Union
import urllib  # TODO: python2/3 compatibility
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord


ATTRIBUTE_SEPARATOR = ";"


class Feature:
    """ A generic GFF feature """
    def __init__(self, location: FeatureLocation, gff_type: str,
                 attributes: Dict[str, str], parent: Optional["Feature"] = None) -> None:
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
        """ Add a generic feature as child """
        assert self.gff_id, "parents must have IDs"
        assert isinstance(feature, Feature), type(feature)
        self.generics.append(feature)
        feature.parent = self

    def __repr__(self) -> str:
        return "%s(%s, %s)" % (self.gff_type, self.gff_id, self.location)

    @property
    def locus_tag(self):
        """ The locus tag of a feature (or the parent's, if not specified and
            the parent has a locus_tag itself
        """
        if "locus_tag" in self.attributes:
            return self.attributes["locus_tag"]
        if self.parent:
            return self.parent.locus_tag
        return None


class Gene(Feature):
    """ A gene-specific variant of Feature. """
    def __init__(self, location: FeatureLocation, gff_type: str, attributes: Dict[str, str]) -> None:
        super().__init__(location, gff_type, attributes)
        self.name = self.attributes.pop("Name", None)
        self.cdss = []  # type: List[CDS]
        self.mrnas = []  # type: List[RNA]
        self.other_rna = []  # type: List[RNA]

    def add_rna(self, rna: "RNA") -> None:
        """ Adds an RNA instance to the gene """
        assert isinstance(rna, RNA)
        if rna.gff_type == "mRNA":
            self.mrnas.append(rna)
        else:
            self.other_rna.append(rna)

    def add_cds(self, cds: "CDS") -> None:
        """ Adds a CDS feature to the gene """
        assert isinstance(cds, CDS)
        self.cdss.append(cds)


class RNA(Feature):
    """ mRNA, rRNA, or tRNA features. All may have a Gene as a parent.
    """
    def __init__(self, location: FeatureLocation, gff_type: str, parent: Gene, attributes: Dict[str, str]) -> None:
        super().__init__(location, gff_type, attributes, parent=parent)
        self.name = self.attributes.pop("Name", None) or (parent.name if parent else None)
        self.cdss = []  # type: List[CDS]

        if parent:
            assert isinstance(parent, Gene), type(parent)

    def add_cds(self, cds: "CDS") -> None:
        """ Adds CDS children to the RNA object """
        assert self.gff_type == "mRNA", self.gff_type
        assert isinstance(cds, CDS)
        self.cdss.append(cds)


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
    """ A simple container for all features within a specific record. """
    def __init__(self, name: str) -> None:
        self.name = name
        self.parent_features = {}  # type: Dict[str, Feature]
        self.all_features = []  # type: List[Feature]
        self.cds_features = []  # type: List[CDS]
        self.cds_features_by_parent = {}  # type: Dict[Feature, CDS]
        self.total_features = 0


def decode(text: str) -> str:
    """ Convert %-style escaping into full text """
    return urllib.parse.unquote(text, errors="strict")


def build_attributes(section: str) -> Dict[str, str]:
    """ Parse the attributes section of a GFF feature line """
    def separate_attributes(text: str) -> Iterator[str]:
        """ Performs the actual separation of attributes given a single block
            of text
        """
        in_quoted_section = False
        separator_in_unclosed_quote = False
        current = []
        for char in text:
            if in_quoted_section:
                if char == '"':  # and not ' because 3' and 5' are common
                    in_quoted_section = False
                    separator_in_unclosed_quote = False
                    continue
                elif char == ATTRIBUTE_SEPARATOR:
                    separator_in_unclosed_quote = True
            elif char == ATTRIBUTE_SEPARATOR:
                yield "".join(current)
                current = []
                continue
            # catch a value with a leading, but unclosed, quote marker
            elif char == '"' and current and current[-1] == "=":
                in_quoted_section = True
                continue
            current.append(char)
        # no text remaining, so quoted sections will never close
        if in_quoted_section and separator_in_unclosed_quote:
            raise ValueError("Attribute with unclosed quote contains an attribute separator:"
                             " %s" % ("".join(current)))
        yield "".join(current)

    attributes = {}
    for attribute in separate_attributes(section):
        if not attribute:
            continue
        parts = attribute.split("=")
        assert len(parts) == 2, "Invalid attribute: %s" % attribute
        attributes[parts[0]] = decode(parts[1])
    return attributes


def interpret_strand(strand: str, strict: bool = False) -> int:
    """ Converts the strand into a BioPython FeatureLocation compatible value
    """
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


def construct_location(raw_start: str, raw_end: str, raw_strand: str,
                       attributes: Dict[str, str], strict: bool = False) -> FeatureLocation:
    """ Converts the raw sections of a GFF line into a FeatureLocation.

        Some attribute keys can modify the location's values.
    """

    try:
        start = ExactPosition(int(raw_start) - 1)  # 0-indexed as FeatureLocation expects
        end = ExactPosition(int(raw_end))
    except ValueError as err:
        raise ValueError("Invalid location values: %s" % str(err))

    if start < 0 or end < 0:
        raise ValueError("Location values cannot be negative")

    strand = interpret_strand(raw_strand, strict=strict)

    # handle ambiguous positions as noted in attributes
    if attributes.get("partial") == "true" and ("start_range" in attributes or "end_range" in attributes):
        attributes.pop("partial")
        start_range = attributes.pop("start_range", "%s,%s" % (start, end))
        end_range = attributes.pop("end_range", "%s,%s" % (start, end))
        if start_range.startswith("."):
            start = BeforePosition(int(start))
        if end_range.endswith("."):
            end = AfterPosition(int(end))

    return FeatureLocation(start, end, strand)


def construct_feature(record: GFFRecord, line: str, parent_lines: Dict[str, List[str]],
                      same_id: List[str], strict: bool = False) -> Feature:
    """ Build a feature from a GFF feature line.

        Args:
            record: the record the feature will belong to
            line: the line from a GFF file
            parent_lines: a dictionary mapping ID attribute to a list of lines
                          from the GFF file with that ID
            same_id: a list of other lines from the file that use the same ID
                     as the given line
            strict: if True, errors will be raised if formats and values don't
                    match the GFF specification
    """
    parts = line.strip().split(sep="\t" if "\t" in line else None, maxsplit=8)
    parts = [part.strip() for part in parts]
    assert len(parts) == 9, parts

    seqid, _, gff_type, start, end, _, strand_dir, _, attribute_string = parts

    # sanity check ensuring that this feature belongs to the provided record
    assert seqid == record.name

    # construct the simple parts
    attributes = build_attributes(attribute_string)
    location = construct_location(start, end, strand_dir, attributes, strict=strict)

    # ensure any specified parent exists
    parent_name = attributes.pop("Parent", None)
    if parent_name is not None and parent_name not in parent_lines:
        raise ValueError("Parent feature referenced but is not known: %s" % parent_name)

    # Some software inserts genes and mRNA with the same ID, but with the mRNA
    # having that same ID as a parent. Discard if discovered.
    if parent_name is not None and parent_name == attributes.get("ID"):
        return None

    # handle compound location genes as provided sometimes by NCBI which have no
    # ID, but are referenced by others (e.g. FNCA010000004.1, SAMN04488589_1460)
    part = attributes.pop("part", None)
    if not attributes.get("ID") and part:
        if not attributes.get("Name"):
            raise ValueError("Feature specified as part but has no reference: %s" % line)
        attributes["ID"] = "%s-%s" % (gff_type, attributes["Name"])

    # ensure the parent feature has been constructed
    parent = record.parent_features.get(parent_name)
    if parent_name and not parent:
        relevant = parent_lines[parent_name]
        parent = construct_feature(record, relevant[0], parent_lines, relevant[1:], strict=strict)

    # now construct the feature itself
    feature = None

    # some pseudogenes use the class instead of an attribute, so make them all
    # attributes
    # TODO: not great for plain GFF parsing, maybe only in genbank conversion
    if gff_type == "pseudogene":
        gff_type = "gene"
        attributes["pseudo"] = "true"

    # if others exist with the same id, merge them here
    if same_id:
        other = construct_feature(record, same_id[0], parent_lines, same_id[1:])
        assert other is not None

        if gff_type != other.gff_type:
            raise ValueError("Two features of different types share the same ID"
                             ": %s  %s" % (other, line))
        if other.location == location:
            raise ValueError("Two features with the same location")

        if other.location != location:
            other.location = other.location + location  # TODO obey the is_ordered=true attribute
            assert isinstance(other.location, CompoundLocation)

        # TODO: attribute merging

        return other

    # since all the same IDs have been processed if they exist, constructing new
    # features should only happen if they have no ID or the ID doesn't exist yet
    if "ID" in attributes:
        assert attributes["ID"] not in record.parent_features, feature

    # finally, create a feature for the given line
    if gff_type == "gene":
        gene = Gene(location, gff_type, attributes)
        feature = gene
    elif gff_type in ["mRNA", "rRNA", "tRNA"]:
        if parent is not None:
            assert isinstance(parent, Gene), str(parent)
        rna = RNA(location, gff_type, parent, attributes)
        feature = rna
        if parent:
            parent.add_rna(rna)
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

    # a final check that all features with an id have a unique id
    if feature.gff_id:
        assert feature.gff_id not in record.parent_features, feature
        record.parent_features[feature.gff_id] = feature

    record.all_features.append(feature)

    return feature


def populate_record(record: GFFRecord, named_lines: Dict[str, List[str]],
                    unnamed_lines: List[str], strict: bool = False) -> None:
    """ Construct Features for each provided line from a GFF file

        Args:
            record: the GFFRecord to populate with features
            named_lines: a dictionary mapping ID attribute to a list of all
                         lines sharing that ID
            unnamed_lines: a list of lines that did not specify an ID
            strict: whether to raise an error if line contents do not match the
                    GFF3 specification
    """
    # build named features first
    for name, lines in named_lines.items():
        # skip if it was generated recursively as a dependency
        if name in record.parent_features:
            continue
        feature = construct_feature(record, lines[0], named_lines, lines[1:], strict=strict)
        assert feature is not None, lines
        assert feature.gff_id == name, "%s != %s" % (feature.gff_id, name)

    # then unnamed features
    for line in unnamed_lines:
        construct_feature(record, line, named_lines, [], strict=strict)


def parse_gff(filename: str, strict: bool = False) -> List[GFFRecord]:
    """ Parses a file in GFF3 format and returns a list of GFFRecords, one for
        each mentioned in the file.
    """
    records_by_name = {}  # type: Dict[str, GFFRecord]
    named_by_record = defaultdict(lambda: defaultdict(list))  # type: Dict[str, Dict[str, List[str]]]
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

            parts = line.strip().split(sep="\t" if "\t" in line else None, maxsplit=8)
            assert len(parts) == 9, parts
            seqid = parts[0]
            # construct a record skeleton
            if seqid not in records_by_name:
                #  if len(seqid) > 16:  # TODO: doesn't work on refseq ids, need to change
                #       raise ValueError("Record name too long: %s" % seqid)
                records_by_name[seqid] = GFFRecord(seqid)

            attributes = build_attributes(parts[8])

            # handle multi-line features that should have IDs but can't
            part = attributes.pop("part", None)
            if not attributes.get("ID") and part:
                if not attributes.get("Name"):
                    raise ValueError("Feature specified as part but has no reference")
                gff_type = parts[2]
                attributes["ID"] = "%s-%s" % (gff_type, attributes["Name"])

            if "ID" in attributes:
                if strict and named_by_record[seqid][attributes["ID"]]:
                    raise ValueError("File contains features with duplicated IDs:"
                                     " %s" % attributes["ID"])
                # skip exact duplicates
                if line not in named_by_record[seqid][attributes["ID"]]:
                    named_by_record[seqid][attributes["ID"]].append(line)
            else:
                anon_by_record[seqid].append(line)

    # populate records
    for record_name, record in records_by_name.items():
        populate_record(record, named_by_record[record_name], anon_by_record[record_name], strict=strict)

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
    """ Converts a single GFFRecord to a SeqRecord """
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
    """ Converts a list of GFFRecord to a list of SeqRecord """
    seq_records = []
    for gff_record in gff_records:
        seq_records.append(convert_record_to_biopython(gff_record))
    return seq_records


if __name__ == "__main__":
    seq_records = convert_gff_to_genbank(parse_gff(sys.argv[1]))
#    from helperlibs.bio import seqio
#    with open("gff_to.gbk", "w") as handle:
#        seqio.write(seq_records, handle, "genbank")
