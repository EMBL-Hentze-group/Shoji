from typing import Dict, List, Tuple, Optional, NamedTuple
from interval import Interval
import logging
# from string import Template

logger = logging.getLogger(__name__)


class Feature(NamedTuple):
    """Feature
    Return an exon or intron data as bed formatted line,
    excluding chromosome column
    """

    chromStart: int
    chromEnd: int
    name: str
    score: int
    strand: str


class Gene:
    """Gene
    class for gene and gene feature data
    """

    def __init__(self) -> None:
        """__init__
        Attributes:
            gene_id (str): The unique identifier for the gene.
            gene_name (str): The name of the gene.
            gene_type (str): The type or category of the gene.
            chrom (str): The chromosome where the gene is located.
            begin (int): The starting position of the gene on the chromosome. Default is 0.
            end (int): The ending position of the gene on the chromosome. Default is 0.
            strand (str): The strand of the DNA (e.g., "+" or "-"). Default is ".".
            _feature_map (Dict[str, Interval]): A dictionary mapping feature names to Interval object.
        """
        self._gene_id: Optional[str] = None
        self._gene_name: Optional[str] = None
        self._gene_type: Optional[str] = None
        self._chrom: Optional[str] = None
        self._start: int = 0
        self._end: int = 0
        self._strand: str = "."
        self._features: Dict[str, Interval] = {}

    @property
    def gene_id(self) -> Optional[str]:
        return self._gene_id

    @gene_id.setter
    def gene_id(self, value: str) -> None:
        self._gene_id = value

    @property
    def gene_name(self) -> Optional[str]:
        return self._gene_name

    @gene_name.setter
    def gene_name(self, value: str) -> None:
        self._gene_name = value

    @property
    def gene_type(self) -> Optional[str]:
        return self._gene_type

    @gene_type.setter
    def gene_type(self, value: str) -> None:
        self._gene_type = value

    @property
    def chrom(self) -> Optional[str]:
        return self._chrom

    @chrom.setter
    def chrom(self, value: str) -> None:
        self._chrom = value

    @property
    def start(self) -> int:
        return self._start

    @start.setter
    def start(self, value: int) -> None:
        self._start = value

    @property
    def end(self) -> int:
        return self._end

    @end.setter
    def end(self, value: int) -> None:
        self._end = value

    @property
    def strand(self) -> str:
        return self._strand

    @strand.setter
    def strand(self, value: str) -> None:
        self._strand = value

    def add_feature(self, feature: str, start: int, end: int) -> None:
        """add_feature _summary_
        Add feature start and end positions
        Args:
            feature: str, feature type
            start: feature start position
            end: feature end position
        """
        try:
            self._features[feature].add(start, end)
        except KeyError:
            self._features[feature] = Interval(start, end)

    @property
    def features(self) -> Dict[str, Interval]:
        return self._features

    @property
    def exons(self) -> List[Tuple[int, int]]:
        """exons exon intervals
        Return merged non overlapping exon intervals
        Raises:
            KeyError: if "exon" key is not found on _features

        Returns:
            List[Tuple[int,int]]
        """
        if len(self._features) == 0:
            logger.warning(
                "Gene %s does not have any features. Returning start and end co-ordinates",
                self.gene_id,
            )
            return [(self.start, self.end)]
        try:
            exons = self.features["exon"].intervals
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k
        return exons

    def tagged_exons(self) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.warning(
                "Gene %s does not have any features. Returning start and end co-ordinates",
                self.gene_id,
            )
            return [
                Feature(
                    chromStart=self.start,
                    chromEnd=self.end,
                    name=self._name_formatter(index=1, n=1, feature_type="exon"),
                    score=0,
                    strand=self.strand,
                )
            ]
        index = self._get_index(len(self.features["exon"]))
        feats: List[Feature] = []
        for i, (start, end) in enumerate(self.features["exon"].intervals):
            feats.append(
                Feature(
                    chromStart=start,
                    chromEnd=end,
                    name=self._name_formatter(
                        index=index[i],
                        n=len(self.features["exon"]),
                        feature_type="exon",
                    ),
                    score=0,
                    strand=self.strand,
                )
            )
        return feats

    def _get_index(self, length: int) -> List[int]:
        """_get_index feature index
        Get feature index
        Args:
            length: int, length/size of feature

        Returns:
            List[int]
        """
        if self.strand == "-":
            return [i for i in range(length, 0, -1)]
        return [i for i in range(1, length + 1)]

    def _name_formatter(
        self,
        index: int,
        n: int,
        feature_type: str,
        split_index: Optional[int] = None,
    ) -> str:
        """_name_formatter format name column for output bed
        Example formats:
            ENSG00000290825.1@DDX11L2@lncRNA@exon@1/3@ENSG00000290825.1:exon0001
            ENSG00000290825.1@DDX11L2@lncRNA@intron@1/2@ENSG00000290825.1:intron0001
        Args:
            index: index of this feature out of total number of intervals
            n: total number of intervals in this feature
            feature_type: type of this feature
            split_index: if feature is an intron and is
                        split into multiple chunks to avoid overlapping an exon from
                        another gene, give the index of the split
                        . Defaults to None.

        Returns:
            str
        """
        index_str: str = f"{index}/{n}"
        feature_id: str = f":{feature_type}{index:04}"
        if split_index is not None:
            index_str = f"{index}-{split_index}/{n}"
            feature_id = f":{feature_type}{index:04}-{split_index}"
        base_dat = [
            self.gene_id,
            self.gene_name,
            self.gene_type,
            feature_type,
            index_str,
            self.gene_id,
            feature_id,
        ]
        return "@".join(base_dat)

    @property
    def introns(self) -> List[Tuple[int, int]]:
        if len(self._features) == 0:
            logger.info("Gene %s does not have any features! No introns", self.gene_id)
            return []
        try:
            introns = self.features["exon"].__invert__().intervals
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k
        return introns

    def tagged_introns(self) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.warning("Gene %s does not have any exons!", self.gene_id)
            return []
        if len(self.features["exon"].intervals) == 1:
            logger.info("Gene %s has only one exon!", self.gene_id)
            return []
        introns: Interval = ~self.features["exon"]
        index = self._get_index(len(introns))
        feats: List[Feature] = []
        for i, (start, end) in enumerate(introns.intervals):
            feats.append(
                Feature(
                    chromStart=start,
                    chromEnd=end,
                    name=self._name_formatter(
                        index=index[i],
                        n=len(introns),
                        feature_type="intron",
                    ),
                    score=0,
                    strand=self.strand,
                )
            )
        return feats

    def remove_intron_exon_overlaps(self, exons: Interval) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.warning("Gene %s does not have any exons!", self.gene_id)
            return []
        if len(self.features["exon"].intervals) == 1:
            logger.info("Gene %s has only one exon!", self.gene_id)
            return []
        introns: Interval = ~self.features["exon"]
        index = self._get_index(len(introns))
        split_feats: List[Feature] = []
        for i, (start, end) in enumerate(introns.intervals):
            intron_i: Interval = Interval(start, end)
            split_introns: Interval = intron_i - exons
            if len(split_introns) == 0:
                logger.info(
                    "Skipping, %s intron %i - %i, complete overlap with an exon",
                    self.gene_id,
                    start,
                    end,
                )
                continue
            split_counter = self._get_index(len(split_introns))
            for j, (split_start, split_end) in enumerate(split_introns.intervals):  # type: ignore
                if introns.__contains__(split_start, split_end):
                    split_index = None
                else:
                    split_index = split_counter[j]
                split_feats.append(
                    Feature(
                        chromStart=split_start,
                        chromEnd=split_end,
                        name=self._name_formatter(
                            index=index[i],
                            n=len(introns),
                            feature_type="intron",
                            split_index=split_index,
                        ),
                        score=0,
                        strand=self.strand,
                    )
                )
        return split_feats
