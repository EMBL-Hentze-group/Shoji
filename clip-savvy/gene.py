from typing import Dict, Set, Tuple


class Gene:
    """Gene
    Class for gene and gene feature data
    """

    def __init__(self) -> None:
        self.gene_id: str = None  # type: ignore
        self.gene_name: str = None  # type: ignore
        self.gene_type: str = None  # type: ignore
        self.chrom: str = None  # type: ignore
        self.begin: int = 0
        self.end: int = 0
        self.strand: str = "."
        self._feature_map: Dict[str, Set[Tuple[int, int]]] = {}

    def add_feature(self, feature: str, begin: int, end: int) -> None:
        """add_feature _summary_
        Add feature begin and end position
        Args:
            feature: str, feature type
            begin: feature begin position
            end: feature end position
        """
        try:
            self._feature_map[feature].add((begin, end))
        except KeyError:
            self._feature_map[feature] = set([(begin, end)])
