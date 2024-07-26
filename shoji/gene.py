from bisect import bisect_left, bisect_right
from typing import Dict, List, NamedTuple, Optional, Tuple

from loguru import logger

from .interval import Interval


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
        if feature not in self._features:
            self._features[feature] = Interval()
        self._features[feature].add(start, end)

    @property
    def features(self) -> Dict[str, Interval]:
        return self._features

    @property
    def exons(self) -> Interval:
        """exons exon intervals
        Return merged non overlapping exon intervals
        Raises:
            KeyError: if "exon" key is not found on _features

        Returns:
            List[Tuple[int,int]]
        """
        if len(self._features) == 0:
            logger.warning(
                f"Gene {self.gene_id} does not have any features. Returning start and end co-ordinates"
            )
            exon = Interval()
            exon.add(self.start, self.end)
            return exon
        try:
            exons = self.features["exon"]
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k
        return exons

    @property
    def exon_start(self) -> int:
        """exon_start
        Start position of the first exon
        Returns:
            int
        """
        return self.features["exon"].first

    @property
    def exon_end(self) -> int:
        """exon_end
        End position of the last exon
        Returns:
            int
        """
        return self.features["exon"].last

    def nexons(self) -> int:
        """nexons number of exons"""
        if len(self._features) == 0:
            return 1
        try:
            return len(self.features["exon"])
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k

    def tagged_exons(self) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.warning(
                f"Gene {self.gene_id} does not have any features. Returning start and end co-ordinates"
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
        for i, (start, end) in enumerate(self.features["exon"]):
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
        feature_id: str = f"{self.gene_id}:{feature_type}{index:04}"
        if split_index is not None:
            index_str = f"{index}-{split_index}/{n}"
            feature_id = f"{self.gene_id}:{feature_type}{index:04}-{split_index}"
        base_dat = [
            self.gene_id,
            self.gene_name,
            self.gene_type,
            feature_type,
            index_str,
            feature_id,
        ]
        return "@".join(base_dat)

    @property
    def introns(self) -> List[Tuple[int, int]]:
        if len(self._features) == 0:
            logger.info(f"Gene {self.gene_id} does not have any features! No introns")
            return []
        try:
            introns = self.features["exon"].__invert__()
        except KeyError as k:
            raise KeyError(
                f"Cannot find 'exon' features for gene {self.gene_id}: {k}"
            ) from k
        return introns

    def tagged_introns(self) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.debug(f"Gene {self.gene_id} does not have any exons!")
            return []
        if len(self.features["exon"]) == 1:
            logger.debug(f"Gene {self.gene_id} has only one exon!")
            return []
        introns: Interval = ~self.features["exon"]
        index = self._get_index(len(introns))
        feats: List[Feature] = []
        for i, (start, end) in enumerate(introns):
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

    def remove_exon_overlapping_intron(self, exons: Interval) -> List[Feature]:
        if (len(self._features) == 0) or ("exon" not in self.features):
            logger.debug(f"Gene {self.gene_id} does not have any exons!")
            return []
        elif len(self.features["exon"]) == 1:
            logger.debug(f"Gene {self.gene_id} has only one exon!")
            return []
        elif self.features["exon"] == exons:
            logger.debug(f"Gene {self.gene_id}: no intron exon overlap")
            return self.tagged_introns()
        introns: Interval = ~self.features["exon"]
        index = self._get_index(len(introns))
        split_introns = introns - exons
        if len(split_introns) == 0:
            logger.info(
                f"Skipping introns from gene: {self.gene_id} introns found: {len(introns)}, complete overlap with exons of other genes",
            )
            return []
        final_exons: List[Feature] = []
        for i, (start, end) in enumerate(introns):
            if (start, end) in split_introns:
                # intron does not overlap any exon
                final_exons.append(
                    Feature(
                        chromStart=start,
                        chromEnd=end,
                        name=self._name_formatter(
                            index=index[i], n=len(introns), feature_type="intron"
                        ),
                        score=0,
                        strand=self.strand,
                    )
                )
                continue
            start_pos = bisect_right(split_introns.ends, start)
            end_pos = bisect_left(split_introns.starts, end, lo=start_pos)
            if start_pos == end_pos:
                # this interval is not present in split_introns
                # comple overlap with an exon
                continue
            for j, (split_start, split_end) in enumerate(
                split_introns.islice(start_pos, end_pos)
            ):
                final_exons.append(
                    Feature(
                        chromStart=split_start,
                        chromEnd=split_end,
                        name=self._name_formatter(
                            index=index[i],
                            n=len(introns),
                            feature_type="intron",
                            split_index=j + 1,
                        ),
                        score=0,
                        strand=self.strand,
                    )
                )
        return final_exons
