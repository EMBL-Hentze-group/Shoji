"""gff3_parser
parse GFF3 formatted files
"""

import logging
import re
from collections import namedtuple
from string import Template
from typing import Callable, Generator, Optional, Set, Tuple

from networkx import DiGraph, ancestors
from xopen import xopen

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
    ) -> None:
        """__init__ _summary_

        Args:
            gff: input gff3 file. Supports .gz files
            out: output file name. Supports writing .gz bed files
            gene_like_features: Gene like features to parse, for example tRNA added from tRNAScan
            idx_id: "ID" attribute in attribute column. Defaults to "ID".
            parent_id: "Parent" attribute in attribute column. Defaults to "Parent".
            gene_id: "gene_id" attribute in attribute column. Defaults to "gene_id".
            gene_name: "gene_name" attribute in attribute column. Defaults to "gene_name".
            gene_type: "gene_type" attribute in attribute column. Defaults to "gene_type".
        """
        self.gff = gff
        self.out = out
        self.gene_like_features = set()
        if gene_like_features is not None:
            self.gene_like_features = [f.lower() for f in gene_like_features]
        self.idx_id = idx_id
        self.parent_id = parent_id
        self.gene_id = gene_id
        self.gene_type = gene_type
        self.gene_name = gene_name
        # regular expression for gff3 attributes
        self._gffre = re.compile(r"(\w+)\=([^;]+)", re.IGNORECASE)
        # unique id template incase self.idx_id attribute is missing
        # format: chromosome|begin|end|strand
        self._uid = Template("$chrom|$begin|$end|$strand")
        # directed graph with gene -> feature dependency
        self._gene_graph = DiGraph()
        # gff3 3rd column entries
        self._feature_names = set()

    def _gff3_reader(
        self,
    ) -> Generator[
        Feature,
        None,
        None,
    ]:
        """_gff3_reader Generic gff3 reader
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
        for feat in self._gff3_reader():
            self._feature_names.add(feat.ftype)  # add feature names
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

    def process(self, feature_type: str = "exon"):
        pass

    def _parse_gene_features(self, feature_type: str = "exon") -> None:
        # feature type checker for gene like features
        feature_checker = self._get_feature_type_checker()
        for feat in self._gff3_reader():
            try:
                begin, end = self._get_coordinates(feat.begin, feat.end)
            except ValueError as vf:
                logging.warning(f"{vf}, skipping feature")
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
                        f"Cannot find 'gene' for feature {feat.ftype}--> {feat.chrom}:{feat.begin}-{feat.end}({feat.strand}) with attributes {feat.attribs}! Skipping"
                    )
                    continue
            elif self._is_parent_node(idx):
                # this is a gene, keep it
                pass
            elif gene_like_feature and (idx not in self._gene_graph):
                # this is a gene like and is not in the graph, keep it
                pass

        pass

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
