from .interval import Interval
from typing import Optional, List, Dict
from bisect import bisect_left, bisect_right


class ChromInterval:
    def __init__(self, strands: Optional[List[str]] = None) -> None:
        self._strand_intervals: Dict[str, Interval] = {}
        if strands is not None:
            self._init_strand_intervals(strands)

    def _init_strand_intervals(self, strands: List[str]) -> None:
        for strand in strands:
            self._strand_intervals[strand] = Interval()

    def add(self, strand: str, start: int, end: int) -> None:
        self._strand_intervals[strand].add(start, end)

    def __len__(self) -> int:
        return len(self._strand_intervals)

    def strand_len(self, strand: str) -> int:
        return len(self._strand_intervals[strand])

    @property
    def intervals(self) -> Dict[str, Interval]:
        return self._strand_intervals

    def find_overlaps(self, strand, start: int, end: int) -> Interval:
        start_pos = bisect_right(self._strand_intervals[strand].ends, start)
        end_pos = bisect_left(self._strand_intervals[strand].starts, end, lo=start_pos)
        return self._strand_intervals[strand].islice(start_pos, end_pos)
