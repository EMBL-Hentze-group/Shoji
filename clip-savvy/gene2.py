from typing import Callable, Dict, Generator, List, Set, Tuple, Optional
from .interval import Interval
import logging

logger = logging.getLogger(__name__)


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
                "Gene %s does not have exons. Returning start and end co-ordinates",
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

    def introns(self) -> List[Tuple[int, int]]:
        if len(self._features) == 0:
            logger.info("Gene %s does not have exons", self.gene_id)
            return []
        try:
            introns = self.features["exon"].__invert__().intervals
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k
        return introns
