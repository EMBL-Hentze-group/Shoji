from typing import Callable, Dict, Generator, List, Set, Tuple


class Gene:
    """Gene
    Class for gene and gene feature data
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
            _feature_map (Dict[str, Set[Tuple[int, int]]]): A dictionary mapping feature names to sets of tuples representing start and end positions of features.
        """
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

    @property
    def features(self) -> Dict[str, Set[Tuple[int, int]]]:
        """features
        Property, return feature dictionary
        Returns:
            Dict[str, Set[Tuple[int, int]]] key: feature type, Set[Tuple[int, int]]: feature positions
        """
        return self._feature_map

    def merge_overlapping_features(self) -> Dict[str, List[Tuple[int, int]]]:
        """merge_overlapping_features
        merge overlapping features and return non overlapping intervals per feature

        Returns:
            Dict[str, List[Tuple[int, int]]] key: feature name, value: list of non overlapping begin and end positions
        """
        if len(self._feature_map) == 0:
            return {"exon": [(self.begin, self.end)]}
        merged_features: Dict[str, List[Tuple[int, int]]] = {}
        for feat, pos in self._feature_map.items():
            merged_features[feat] = _merge_pos(pos)
        return merged_features

    def merge_overlapping_exons_add_introns_add_tags(
        self,
    ) -> List[Tuple[int, int, str, int, int]]:
        """merge_overlapping_exons_add_introns_add_tags merge overlapping exons
        Merge overlapping exons, add introns and add tags
        Tags:
            - Position 2:
                * This tag will be either "exon" or "intron".
            - Position 3:
                * m th feature out of N feaures
            - Position 4:
                * N total number of reatures
        Raises:
            KeyError: raise KeyError if "exon" is not found in self._feature_map

        Returns:
            List[Tuple[begin,end,ftype,ix,alln]]
                begin: begin co-ordinate of this feature
                end: end co-ordinate of this feaure
                ftype: type of this feature (exon or intron)
                idx: this feature is in the ixth position out of all features of this type
                alln: total number of features of this type
        """
        if len(self._feature_map) == 0:
            # no exons, return gene begin and end position as exons
            return [(self.begin, self.end, "exon", 1, 1)]
        try:
            exons = _merge_pos(self._feature_map["exon"])
        except KeyError as k:
            raise KeyError(
                f"Cannot find exon features for gene {self.gene_id}, {k}"
            ) from k
        if len(exons) == 1:
            # There is only one exon after merging, no introns needed
            return [(exons[0][0], exons[0][1], "exon", 1, 1)]
        tagged_pos: List[Tuple[int, int, str, int, int]] = list()
        tagger = self._get_tagger()
        exon_tagger = tagger(len(exons))
        for i, tag in enumerate(exon_tagger):
            tagged_pos.append((exons[i][0], exons[i][1], "exon", tag[0], tag[1]))
        introns = self._get_intron_pos(exons)
        intron_tagger = tagger(len(introns))
        for i, tag in enumerate(intron_tagger):
            tagged_pos.append((introns[i][0], introns[i][1], "intron", tag[0], tag[1]))
        return sorted(tagged_pos)

    @staticmethod
    def _get_intron_pos(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """_get_intron_pos generate intron co-ordinates
        Generate intron co-ordinates from a list of exon co-ordinates
        Args:
            exons:  List[Tuple[begin,end]] list of exon begin and end positions

        Returns:
            List[Tuple[begin,end]] list of intron begin and end positions
        """
        prev_end = 0
        introns: List[Tuple[int, int]] = list()
        for i, (begin, end) in enumerate(exons):
            if i == 0:
                prev_end = end
                continue
            if begin > prev_end:
                introns.append((prev_end, begin))
            prev_end = end
        return introns

    def _get_tagger(self) -> Callable:
        """_get_tagger get strand appropriate tagger
        If there are more than one feature per feature type,
        add a tagger with the format m/N where
            m: the position of the current feature 1,2,3...
            N: total number of features
        Returns:
            Callable, strand appropriate tagger
        """

        def tagger_plus(nmax: int) -> Generator[Tuple[int, int], None, None]:
            """tagger_plus plus strand tagger
            for nmax features yield values 1,nmax, 2,nmax....nmax,nmax
            Args:
                nmax: int, maximum number of features

            Yields:
                Tuple, 1,nmax, 2,nmax....nmax,nmax
            """
            i: int = 1
            while i <= nmax:
                yield (i, nmax)
                i += 1

        def tagger_minus(nmax: int) -> Generator[Tuple[int, int], None, None]:
            """tagger_minus minus strand tagger

            for nmax features yield values nmax,nmax, nmax-1,nmax...1,nmax
            Args:
                nmax: int, maximum number of features

            Yields:
                Tuple, nmax,nmax, nmax-1,nmax...1,nmax
            """
            i: int = nmax
            while i >= 1:
                yield (i, nmax)
                i -= 1

        if self.strand == "-":
            return tagger_minus
        return tagger_plus


def _merge_pos(pos: Set[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """_merge_pos
    Merge overlapping intevals and return non overlapping intervals

    Args:
        pos: Set[Tuple[int,int]] set of overlapping intervals

    Returns:
        List[Tuple[int,int]] set of non overlapping sorted intervals
    """
    spos = sorted(pos)
    if len(spos) == 1:
        return spos
    posi = iter(spos)
    begin, end = next(posi)
    mergedpos = []
    for cbegin, cend in posi:
        if cbegin >= end:
            # positions are zero based
            mergedpos.append((begin, end))
            begin, end = cbegin, cend
        else:
            end = max(end, cend)
    mergedpos.append((begin, end))
    return mergedpos
