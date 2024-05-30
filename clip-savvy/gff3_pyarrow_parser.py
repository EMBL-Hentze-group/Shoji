import heapq
import logging
import re
from string import Template
from typing import Callable, Dict, List, Optional, Set

import pyarrow as pa
from networkx import DiGraph, ancestors

from . import pyarrow_reader as pr
from .gene import Gene
from .output import output_writer

logger = logging.getLogger(__file__)


class GFF3parser:
    """GFF3parser
    Parse GFF3 files and extract exon features for genes
    experimental module using pyarrow
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
            use_tabix: boolean, if True, use bgzip to compress the file and then index the file using tabix
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
        self._uid = Template("${chrom}|${start}|${end}|${strand}")
        # directed graph with gene -> feature dependency
        self._gene_graph: DiGraph = DiGraph()
        # save per gene features in a dictionary
        self._gene_feature_map: Dict[str, Gene] = {}
        # List of features from pyarrow table
        self._feats: List[Dict] = []
        # directed graph with gene -> feature dependency
        self._gene_graph: DiGraph = DiGraph()
        # save per gene features in a dictionary
        self._gene_feature_map: Dict[str, Gene] = {}
        # heap to store and sort co-ordinates
        self._heap = []
        # output write to write outputs
        self._ow = output_writer(self.out, use_tabix=self.use_tabix, preset="bed")

    def process(self, feature_type: str = "exon") -> None:
        """process process gff3 file
        parse gff3 file for the given feature type
        Args:
            feature_type: str, gff3 third column feature to parse. Defaults to "exon".
        """
        # gff3 pyarrow table
        gff3_pa = self._read_gff3()
        if gff3_pa.shape[0] == 0:
            raise RuntimeError(f"Cannot parse annotation features from {self.gff}!")
        chroms: List[str] = sorted(gff3_pa.column("seqname").unique().tolist())
        logger.info(
            "Found %s features from %s chromosomes in %s",
            f"{gff3_pa.shape[0]:,}",
            f"{len(chroms):,}",
            self.gff,
        )
        with self._ow(self.out) as _fh:
            # _fh: Union[_io.TextIOWrapper,tempfile._TemporaryFileWrapper]
            for chrom in chroms:
                logging.info("Parsing gene and feature data from %s", chrom)
                self._feats: List[Dict] = gff3_pa.filter(
                    pa.compute.field("seqname") == chrom
                ).to_pylist()
                # generate dependancy
                self._gene_feature_dependancy()
                logging.info(
                    "   %s: found %s genes and features",
                    chrom,
                    f"{len(self._gene_graph):,}",
                )
                # parse features
                self._parse_gene_features(feature_type=feature_type)
                logging.info(
                    "   %s: mapped features to %s genes",
                    chrom,
                    f"{len(self._gene_feature_map):,}",
                )
                # sort features
                self._sort_coordinates()
                # now write to file
                self._write(_fh, chrom)

        logging.info("Wrote sorted featues to %s", self.out)

    def _read_gff3(self):
        """_read_gff3 gff3 reader
        Read GFF3 file using pyarrow
        Returns:
            pyarrow Table
        """
        gff3_reader = pr.Reader(self.gff)
        return gff3_reader.gff().drop_columns(["source", "score", "frame"])

    def _gene_feature_dependancy(self) -> None:
        """_gene_feature_dependancy generate gene feature dependancy per chromosome
        Helper function.
        generate gene dependancy graph using networkx per chromosome features
        all nodes in this graph with in_degree 0 should be all genes and
        all nodes with out_degree 0 should be all gene features
        """
        self._gene_graph = DiGraph()
        for f in self._feats:
            attribs: Dict[str, str] = dict(re.findall(self._gffre, f["attributes"]))
            if self.parent_id not in attribs:
                continue
            if self.idx_id in attribs:
                uid: str = attribs[self.idx_id]
            else:
                # use position as unique id if ID attribute is not found
                # format: chromosome|begin|end|strand
                uid: str = self._uid.substitute(
                    chrom=f["seqid"],
                    start=str(f["start"]),
                    end=str(f["end"]),
                    strand=f["strand"],
                )
            self._gene_graph.add_edge(attribs[self.parent_id], uid)

    def _parse_gene_features(self, feature_type: str = "exon") -> None:
        """_parse_gene_features parse gene features per chromosome
        Helper function, parse gene features given gene_graph dependency graph and
        list of features
        Args:
            gene_graph: DiGraph object, gene - feature dependancy per chromosome
            feats: List[Dict] feature list
            feature_type: str, gff3 third column feature to parse. Defaults to "exon".

        Returns:
            Dict[str, Gene]
        """

        self._gene_feature_map = {}
        feature_checker = self._get_feature_type_checker()
        for f in self._feats:
            start = f["start"] - 1  # to 0 based format
            attribs = dict(re.findall(self._gffre, f["attributes"]))
            try:
                uid: str = attribs[self.idx_id]
            except KeyError:
                uid: str = self._uid.substitute(
                    chrom=f["seqname"],
                    start=str(f["start"]),
                    end=str(f["end"]),
                    strand=f["strand"],
                )
            gene_like_feature: bool = feature_checker(f["type"])
            # what all to keep ?
            if (f["type"] == feature_type) and self._is_leaf_node(uid):
                genes: Set[str] = self._get_genes(uid)
                if len(genes) == 0:
                    logging.warning(
                        "Cannot find 'gene' for feature %s--> %s:%i-%i(%s) with attributes %s! Skipping",
                        f["type"],
                        f["seqname"],
                        f["start"],
                        f["end"],
                        f["strand"],
                        f["attributes"],
                    )
                    continue
                for g in genes:
                    if g not in self._gene_feature_map:
                        self._gene_feature_map[g] = Gene()
                    self._gene_feature_map[g].add_feature(f["type"], start, f["end"])
            elif (self._is_parent_node(uid)) or (
                gene_like_feature and (uid not in self._gene_graph)
            ):
                # either a gene or
                # a gene like feature with id not found in the gene dependanch graph
                self._add_gene_info(f["seqname"], start, f["end"], f["strand"], attribs)

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

    def _is_leaf_node(self, idx: str) -> bool:
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

    def _sort_coordinates(self) -> None:
        """_sort_coordinates sort co-ordinates
        Helper function
        Add features to heap to sort
        """
        # empty heap
        self._heap = []
        # Template for the name column in final bed format
        name_col = Template(
            "${idx}@${name}@${gtype}@${feat}@${ix}/${alln}@${idx}:${feat}${padix}"
        )
        for dat in self._gene_feature_map.values():
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
                    self._heap,
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
                        0,
                        strand,
                    ),
                )

    def _write(self, fh, chrom: str) -> None:
        """_write write to file
        Pop items from the list and write to the file handle
        Args:
            fh: file handle, Union[_io.TextIOWrapper,tempfile._TemporaryFileWrapper]
            chrom: chromosome name
        """
        while True:
            try:
                (begin, end, name, score, strand) = heapq.heappop(self._heap)
                fh.write(f"{chrom}\t{begin}\t{end}\t{name}\t{score}\t{strand}\n")
            except IndexError:
                fh.flush()
                break


# import sys


# def main():
#     root = logging.getLogger()
#     root.setLevel(logging.DEBUG)

#     handler = logging.StreamHandler(sys.stdout)
#     handler.setLevel(logging.DEBUG)
#     root.addHandler(handler)

#     gencode = GFF3parser(
#         gff="/workspaces/clip_savvy/test_data/gencode.v42.annotation.plus.tRNAs.sorted.gff3",
#         out="/workspaces/clip_savvy/test_data/gencode.v42.annotation.plus.tRNAs.new_format.bed.gz",
#         parent_id="Parent",
#         idx_id="ID",
#         gene_name="gene_name",
#         gene_id="gene_id",
#         gene_type="gene_type",
#         gene_like_features=set(["tRNA"]),
#         use_tabix=True,
#     )
#     gencode.process()
#     # ensembl = GFF3parser(
#     #     gff="/workspaces/clip_savvy/test_data/Mus_musculus.GRCm39.111.gff3.gz",
#     #     out="/workspaces/clip_savvy/test_data/Mus_musculus.GRCm39.pa_bz2.bed.bz2",
#     #     parent_id="Parent",
#     #     idx_id="ID",
#     #     gene_name="Name",
#     #     gene_id="gene_id",
#     #     gene_type="biotype",
#     #     gene_like_features=None,
#     #     use_tabix=False,
#     # )
#     # ensembl.process()


# if __name__ == "__main__":
#     main()
