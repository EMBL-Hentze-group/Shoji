"""gff3_parser
parse GFF3 formatted files
"""

import heapq
import logging
import re
from collections import namedtuple
from string import Template
from typing import Callable, Dict, Generator, Optional, Set, Tuple

from networkx import DiGraph, ancestors
from xopen import xopen

from .gene import Gene
from .output import output_writer

logger = logging.getLogger(__file__)

# named tuple for major columns in gff3
Feature = namedtuple("feature", ["chrom", "ftype", "begin", "end", "strand", "attribs"])


class GFF3parser:
    """GFF3parser
    Parse GFF3 files and extract exons features for genes
    """

    def __init__(
        self,
        gff: str,
        out: str,
        gene_like_features: Optional[Set[str]],
        idx_id: str = "ID",
        parent_id: str = "Parent",
        gene_id: str = "gene_id",
        gene_name: str = "gene_name",
        gene_type: str = "gene_type",
        use_tabix: bool = True,
    ) -> None:
        """__init__
        Args:
            gff: input gff3 file. Supports .gz files
            out: output file name. Supports writing .gz, .bz2 or .xz compressed bed files. Additional support for tabix indexing
            gene_like_features: Gene like features to parse, for example tRNA added from tRNAScan
            idx_id: "ID" attribute in attribute column. Defaults to "ID".
            parent_id: "Parent" attribute in attribute column. Defaults to "Parent".
            gene_id: "gene_id" attribute in attribute column. Defaults to "gene_id".
            gene_name: "gene_name" attribute in attribute column. Defaults to "gene_name".
            gene_type: "gene_type" attribute in attribute column. Defaults to "gene_type".
            use_tabix: boolean, if True, use tabix to index gzip file
        """
        self.gff: str = gff
        self.out: str = out
        self.gene_like_features: Set[str] = set()
        if gene_like_features is not None:
            self.gene_like_features = {f.lower() for f in gene_like_features}
        self.idx_id: str = idx_id
        self.parent_id: str = parent_id
        self.gene_id: str = gene_id
        self.gene_type: str = gene_type
        self.gene_name: str = gene_name
        self.use_tabix: bool = use_tabix
        # regular expression for gff3 attributes
        self._gffre = re.compile(r"(\w+)\=([^;]+)", re.IGNORECASE)
        # unique id template incase self.idx_id attribute is missing
        # format: chromosome|begin|end|strand
        self._uid = Template("${chrom}|${begin}|${end}|${strand}")
        # directed graph with gene -> feature dependency
        self._gene_graph: DiGraph = DiGraph()
        # save per gene features in a dictionary
        self._gene_feature_map: Dict[str, Gene] = {}

    def process(self, feature_type: str = "exon") -> None:
        """process process gff3 file
        parse gff3 file for the given feature type
        Args:
            feature_type: str, gff3 third column feature to parse. Defaults to "exon".
        """
        self._gene_feature_dependancy()
        self._parse_gene_features(feature_type=feature_type)
        self._sort_n_write()

    def _gff3_reader(
        self,
    ) -> Generator[
        Feature,
        None,
        None,
    ]:
        """_gff3_reader Generic gff3 reader
        Helper function
        Read GFF3 file and yield a named tupele

        Yields:
            named tuple with following fields: chrom, ftype, begin, end, strand and attribs
        """

        with xopen(self.gff) as _gh:
            for f in _gh:
                if f[0] == "#":
                    continue
                fdat = f.strip().split("\t")
                if len(fdat) < 9:
                    continue
                yield Feature(
                    fdat[0],  # chromosome
                    fdat[2],  # type
                    fdat[3],  # begin
                    fdat[4],  # end
                    fdat[6],  # strand
                    fdat[-1],  # attribute
                )

    def _gene_feature_dependancy(self) -> None:
        """_gene_feature_dependancy
        Helper function.
        Parse gff3 file and generate gene dependancy graph using networkx
        all nodes in this graph with in_degree 0 should be all genes and
        all nodes with out_degree 0 should be all gene features
        reverse map each feature with 0 out_degree node to its corresponding gene
        Args:
            None
        """
        logger.info("Generating gene feature dependancy from %s", self.gff)
        for feat in self._gff3_reader():
            # self._feature_names.add(feat.ftype)  # add feature names
            attribs = dict(re.findall(self._gffre, feat.attribs))
            if self.parent_id not in attribs:
                continue
            if self.idx_id in attribs:
                uid: str = attribs[self.idx_id]
            else:
                # use position as unique id if idx attrib is not found
                # chromosome|begin|end|strand
                uid: str = self._uid.substitute(
                    chrom=feat.chrom, begin=feat.begin, end=feat.end, strand=feat.strand
                )
            self._gene_graph.add_edge(attribs[self.parent_id], uid)
        if len(self._gene_graph) == 0:
            raise RuntimeWarning(
                f"Cannot parse gene and features from {self.gff}. Check your input file!"
            )

    def _parse_gene_features(self, feature_type: str = "exon") -> None:
        """_parse_gene_features parse gene features
        Helper function, parse gff3 file for the given feature type
        Args:
            feature_type: str, gff3 third column feature to parse. Defaults to "exon".
        """
        # feature type checker for gene like features
        logger.info("Parsing gene and %s features from %s", feature_type, self.gff)
        feature_checker = self._get_feature_type_checker()
        for feat in self._gff3_reader():
            try:
                begin, end = self._get_coordinates(feat.begin, feat.end)
            except ValueError as vf:
                logging.warning("%s, skipping feature", str(vf))
                continue
            attribs = dict(re.findall(self._gffre, feat.attribs))
            try:
                idx = attribs[self.idx_id]
            except KeyError:
                idx = self._uid.substitute(
                    chrom=feat.chrom, begin=feat.begin, end=feat.end, strand=feat.strand
                )
            gene_like_feature = feature_checker(feat.ftype)
            # what all to keep ?
            if (feat.ftype == feature_type) and self.is_leaf_node(idx):
                # this is a feature type and is a leaf node, keep it
                genes = self._get_genes(idx)
                if len(genes) == 0:
                    logging.warning(
                        "Cannot find 'gene' for feature %s--> %s:%s-%s(%s) with attributes %s! Skipping",
                        feat.ftype,
                        feat.chrom,
                        feat.begin,
                        feat.end,
                        feat.strand,
                        feat.attribs,
                    )
                    continue
                for g in genes:
                    try:
                        self._gene_feature_map[g].add_feature(feat.ftype, begin, end)
                    except KeyError:
                        self._gene_feature_map[g] = Gene()
                        self._gene_feature_map[g].add_feature(feat.ftype, begin, end)
            elif (self._is_parent_node(idx)) or (
                gene_like_feature and (idx not in self._gene_graph)
            ):
                # either a gene or
                # a gene like feature with id not found in the gene dependanch graph
                self._add_gene_info(feat.chrom, begin, end, feat.strand, attribs)

    def _get_coordinates(self, sbegin: str, send: str) -> Tuple[int, int]:
        """_get_coordinates
        Return co-ordinates as integers, convert begin position to 0 based system
        Args:
            sbegin: begin position as a string
            send: end position as a string

        Raises:
            ValueError: raise if begin position cannot be converted to integer
            ValueError: raise if end position cannot be converted to integer

        Returns:
            Tuple[int,int] begin and end co-ordinates as a tuple
        """
        try:
            # convert from 1 based co-ordinate system to 0 based for BED
            begin = int(sbegin) - 1
        except ValueError as v:
            raise ValueError(
                f"Cannot parse feature begin position from {sbegin}, {v}"
            ) from v
        try:
            end = int(send)
        except ValueError as v:
            raise ValueError(
                f"Cannot parse feature end position from {send}, {v}"
            ) from v
        return (begin, end)

    def _get_feature_type_checker(self) -> Callable:
        """_get_feature_type_checker
        Helper function
        Reture an approprite feature type checker function
        """

        def no_gene_like_feature(ftype: str) -> bool:
            """no_gene_like_feature
            Helper function
            self.gene_like_features is empty, always return False
            Args:
                ftype: feature type, str, 3rd column in gff3 file

            Returns:
                False
            """
            return False

        def check_gene_like_feature(ftype: str) -> bool:
            """check_gene_like_feature
            Helper function
            check if ftype is present in self.gene_like_features and return
            Args:
                ftype: feature type, str, 3rd column in gff3 file

            Returns:
                boolean
            """
            if ftype.lower() in self.gene_like_features:
                return True
            return False

        if len(self.gene_like_features) == 0:
            return no_gene_like_feature
        return check_gene_like_feature

    def is_leaf_node(self, idx: str) -> bool:
        """is_leaf_node
        Helper function
        Check whether a given node is a leaf node
        Args:
            idx: str, unique id of the featre

        Returns:
            boolean
        """
        if idx not in self._gene_graph:
            return False
        in_degree = self._gene_graph.in_degree(idx)
        out_degree = self._gene_graph.out_degree(idx)
        if in_degree > 0 and out_degree == 0:
            return True
        return False

    def _is_parent_node(self, idx: str) -> bool:
        """is_parent_node
        Helper function
        Check whether the given node is a parent node
        Args:
            idx: str, unique id of the feature

        Returns:
            boolean
        """
        if idx not in self._gene_graph:
            return False
        in_degree = self._gene_graph.in_degree(idx)
        out_degree = self._gene_graph.out_degree(idx)
        if in_degree == 0 and out_degree > 0:
            return True
        return False

    def _get_genes(self, idx: str) -> Set[str]:
        """_get_genes
        Helper function
        Get all the gene ancestors of the given node.
        Gene ancestors are the parent nodes with out_degree == 0
        Args:
            idx: str, node id, unique id from GFF3

        Returns:
            Set[str], unique ids of all gene ancestors of this node
        """
        genes = set()
        for anc in ancestors(self._gene_graph, idx):
            if self._gene_graph.in_degree(anc) == 0:
                genes.add(anc)
        return genes

    def _add_gene_info(
        self, chrom: str, begin: int, end: int, strand: str, attribs: Dict[str, str]
    ) -> None:
        """_add_gene_info add gene info
        Helper function, add gene info to self._gene_feature_map dictionary
        Args:
            chrom: str, chromosome name
            begin: int, gene begin position
            end: int, gene end position
            strand: str, strand info
            attribs: Dict[str,str] attribute dictionary, formatted 9th column

        Raises:
            KeyError: if self.gene_id is not found in the attribute dictionary
        """
        idx = attribs[self.idx_id]
        if idx not in self._gene_feature_map:
            # if loop here is far more elegant than try except block
            self._gene_feature_map[idx] = Gene()
        try:
            self._gene_feature_map[idx].gene_id = attribs[self.gene_id]
        except KeyError as k:
            raise KeyError(f"Cannot find 'gene_id' for ID {idx}! {k}") from k
        self._gene_feature_map[idx].gene_name = self._attrib_getter(
            attribs, self.gene_name, attribs[self.gene_id]
        )
        self._gene_feature_map[idx].gene_type = self._attrib_getter(
            attribs, self.gene_type, "unknown"
        )
        self._gene_feature_map[idx].chrom = chrom
        self._gene_feature_map[idx].begin = begin
        self._gene_feature_map[idx].end = end
        self._gene_feature_map[idx].strand = strand

    def _attrib_getter(
        self, attrib: Dict[str, str], attrib_key: str, default: str
    ) -> str:
        """_attrib_getter
        Helper function, return value from attrib dictionary based on key.
        If key is not found, return default value
        Args:
            attrib: Attribute dictionary
            attrib_key: key in the attribute dictionary
            default: default string to return if key is not found

        Returns:
            str
        """
        try:
            return attrib[attrib_key]
        except KeyError:
            return default

    def _sort_n_write(self) -> None:
        """_sort_n_write sort gene features and write to file
        Helper function, sort the gene features based on chromsom, begin co-ordinate and write to self.out
        """
        # Template for the name column in final bed format
        name_col = Template(
            "${idx}@${name}@${gtype}@${feat}@${ix}/${alln}@${idx}:${feat}${padix}"
        )
        # output writer function
        ow = output_writer(self.out, use_tabix=self.use_tabix, preset="bed")
        # keep sorted features
        chrom_sorted = {}
        for dat in self._gene_feature_map.values():
            if dat.chrom not in chrom_sorted:
                chrom_sorted[dat.chrom] = []
            gene_id = dat.gene_id
            gene_name = dat.gene_name
            gene_type = dat.gene_type
            strand = dat.strand
            for (
                begin,
                end,
                ftype,
                ix,
                alln,
            ) in dat.merge_overlapping_exons_add_introns_add_tags():
                heapq.heappush(
                    chrom_sorted[dat.chrom],
                    (
                        begin,
                        end,
                        name_col.substitute(
                            idx=gene_id,
                            name=gene_name,
                            gtype=gene_type,
                            feat=ftype,
                            ix=ix,
                            alln=alln,
                            padix=f"{ix:04}",
                        ),
                        strand,
                    ),
                )
        with ow(self.out) as _fh:
            for chrom in sorted(chrom_sorted.keys()):
                while True:
                    try:
                        (begin, end, name, strand) = heapq.heappop(chrom_sorted[chrom])
                        _fh.write(f"{chrom}\t{begin}\t{end}\t{name}\t0\t{strand}\n")
                    except IndexError:
                        _fh.flush()
                        break


# def print_stuff(gff3: GFF3parser) -> None:
#     for idx, dat in gff3._gene_feature_map.items():
#         print(
#             idx,
#             dat.chrom,
#             dat.begin,
#             dat.end,
#             dat.strand,
#             dat.gene_id,
#             dat.gene_name,
#             dat.gene_type,
#         )
#         print(dat.features)
#         print(dat.merge_overlapping_exons_add_introns_add_tags())
#     print(len(gff3._gene_feature_map))


# import sys


# def main():
#     root = logging.getLogger()
#     root.setLevel(logging.DEBUG)

#     handler = logging.StreamHandler(sys.stdout)
#     handler.setLevel(logging.DEBUG)
#     root.addHandler(handler)

# gencode = GFF3parser(
#     gff="/workspaces/clip_savvy/test_data/gencode.v42.annotation.plus.tRNAs.sorted.gff3",
#     out="/workspaces/clip_savvy/test_data/gencode.v42.annotation.plus.tRNAs.bed",
#     parent_id="Parent",
#     idx_id="ID",
#     gene_name="gene_name",
#     gene_id="gene_id",
#     gene_type="gene_type",
#     gene_like_features=set(["tRNA"]),
#     use_tabix=False,
# )
# gencode.process()

#     # print_stuff(gencode)

# ensembl = GFF3parser(
#     gff="/workspaces/clip_savvy/test_data/Mus_musculus.GRCm39.111.gff3.gz",
#     out="/workspaces/clip_savvy/test_data/Mus_musculus.GRCm39.test_tabix.bed.gz",
#     parent_id="Parent",
#     idx_id="ID",
#     gene_name="Name",
#     gene_id="gene_id",
#     gene_type="biotype",
#     gene_like_features=None,
#     use_tabix=True,
# )
# ensembl.process()
# print_stuff(ensembl)


#     # ncbi = GFF3parser(
#     #     gff="/workspaces/clip_savvy/test_data/GCA_000001405.29_GRCh38.p14_genomic.gff.gz",
#     #     out="blah",
#     #     parent_id="Parent",
#     #     idx_id="ID",
#     #     gene_name="Name",
#     #     gene_id="Dbxref",
#     #     gene_type="gene_biotype",
#     #     gene_like_features=None,
#     # )
#     # ncbi.process()
#     # print_stuff(ncbi)
#     # print(ncbi._gene_feature_map["gene-ND1"].features)


# if __name__ == "__main__":
#     main()
