from typing import List, Optional, Tuple
from bisect import bisect_left, bisect_right
from sortedcontainers import SortedList, SortedSet
from operator import itemgetter


class Interval(SortedSet):
    def __init__(self, start: Optional[int] = None, end: Optional[int] = None):
        super().__init__()
        if start is not None and end is not None:
            self.add(start, end)

    def __repr__(self) -> str:
        if len(self) == 0:
            return "[]"
        else:
            return "[" + ", ".join([f"({i[0]}, {i[1]})" for i in self]) + "]"

    def add(self, start: int, end: int) -> None:
        if start < 0 or end < 0:
            raise ValueError(
                f"'start' and 'end' values MUST be >= 0. Found {start} and {end}!"
            )
        if end - start <= 0:
            raise ValueError(
                f"'end' MUST be larger than 'start'. Found {start} and  {end}!"
            )
        if (start, end) not in self:
            # find the position where the existing end is >= the start
            start_pos = bisect_right(self._list, start, key=itemgetter(1))
            # find the position where the existing start is <= the end
            # no need to search from 0th index, start from start_pos
            end_pos = bisect_left(self._list, end, lo=start_pos, key=itemgetter(0))
            if start_pos == end_pos:
                # new interval does not overlap any existing intervals
                # or the list is empty
                super().add((start, end))
            else:
                # end_pos will always be the "next" index where this end can be inserted,
                # but this does not mean that this end overlaps with the interval positions at end_pos
                for i in range(end_pos - 1, start_pos - 1, -1):
                    # new interval overlaps with this/these existing interval(s)
                    start = min(start, self._list[i][0])
                    end = max(end, self._list[i][1])
                    _ = self.pop(i)
                super().add((start, end))

    @property
    def first(self) -> int:
        return self[0][0]

    @property
    def last(self) -> int:
        return self[-1][1]

    @property
    def empty(self) -> bool:
        """empty _summary_
        shamelessly stolen from
        https://github.com/AlexandreDecan/portion/blob/master/portion/interval.py#L176
        credits to the author
        Returns:
            _description_
        """
        return len(self) == 0

    def __contains__(self, interval: Tuple[int, int]) -> bool:
        """__contains__ check if an entry exists
        Helper function
        source: https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
        Args:
            start: start position of the interval
            end: end position of the interval

        Returns:
            boolean, True if the interval already exists
        """
        pos = self.bisect_left(interval)
        try:
            if pos < len(self) and self[pos] == interval:
                return True
        except IndexError:
            return False
        else:
            return False

    def __iter__(self):
        return iter(self._list)

    def __len__(self):
        return len(self._list)

    def __getitem__(self, index):
        return self._list[index]

    def __reversed__(self):
        return reversed(self._list)

    def __delitem__(self, index):
        del self._list[index]

    def __hash__(self):
        return hash(tuple(self._list))

    def __eq__(self, other):
        if not isinstance(other, Interval):
            return NotImplemented
        return hash(self) == hash(other)

    def __ne__(self, other):
        if not isinstance(other, Interval):
            return NotImplemented
        return not self == other

    def __lt__(self, other: "Interval") -> bool:
        if not isinstance(other, Interval):
            return NotImplemented
        if self.empty or other.empty:
            return False
        return self.last <= other.first

    def __gt__(self, other: "Interval") -> bool:
        if not isinstance(other, Interval):
            return NotImplemented
        if self.empty or other.empty:
            return False
        return self.first >= other.last

    def __call__(self):
        return self._list

    def _merge_overlap(self):
        """_merge_overlap merge overlapping intervals"""
        if len(self) > 1:
            merged = SortedList()
            start, end = self[0]
            for i in self[1:]:
                if i[0] < end:
                    end = max(end, i[1])
                else:
                    merged.add((start, end))
                    start, end = i
            merged.add((start, end))
            self._list = merged

    @staticmethod
    def _intersects(
        start: int, end: int, istart: int, iend: int
    ) -> Optional[Tuple[int, int]]:
        """_intersects check intersection
        Return True if two given intervals overlap
        Args:
            start: int, start position of the first interval
            end: int, end position of the first interval
            istart: int, start position of the second interval
            iend: int, end position of the second interval

        Returns:
            bool
        """
        start_max: int = max(start, istart)
        end_min: int = min(end, iend)
        if start_max < end_min:
            return start_max, end_min
        return None

    def __and__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        if (self.empty or other.empty) or (self > other) or (self < other):
            return self.__class()
        intersections: Interval = Interval()
        for ostart, oend in other:
            start_pos = bisect_right(self._list, ostart, key=itemgetter(1))
            end_pos = bisect_left(self._list, oend, lo=start_pos, key=itemgetter(0))
            if start_pos == end_pos:
                # this interval fits between two intervals in self
                # or is either smaller or larger than all intervals in self
                continue
            for i in range(end_pos - 1, start_pos - 1, -1):
                ostart = max(ostart, self._list[i][0])
                oend = min(oend, self._list[i][1])
            intersections.add(ostart, oend)
        return intersections

    def __or__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        if self.empty and other.empty:
            return self.__class__()
        elif self.empty and (not other.empty):
            return other
        elif (not self.empty) and (other.empty):
            return self
        else:
            unions: Interval = Interval()
            for ostart, oend in other:
                start_pos = bisect_right(self._list, ostart, key=itemgetter(1))
                end_pos = bisect_left(self._list, oend, lo=start_pos, key=itemgetter(0))
                if start_pos == end_pos:
                    # this interval fits between two intervals in self
                    # or is either smaller or larger than all intervals in self
                    unions.add(ostart, oend)
                    continue
                for i in range(end_pos - 1, start_pos - 1, -1):
                    ostart = min(ostart, self._list[i][0])
                    oend = max(oend, self._list[i][1])
                unions.add(ostart, oend)
            return unions

    def __invert__(self) -> "Interval":
        if len(self) <= 1:
            return self.__class__()
        prev_ends: List[int] = [i[1] for i in self][:-1]
        next_starts: List[int] = [i[0] for i in self][1:]
        inverted: Interval = Interval()
        for end, start in zip(prev_ends, next_starts):
            if end < start:
                inverted.add(end, start)
        return inverted

    def __sub__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        common = self & other
        if common.empty:
            return self
        difference: Interval = Interval()
        for start, stop in self:
            intersects: bool = False
            for cstart, cstop in common:
                if overlap := self._intersects(start, stop, cstart, cstop):
                    intersects = True
                    if abs(start - overlap[0]) > 0:
                        cstart, cend = sorted((start, overlap[0]))
                        difference.add(cstart, cend)
                    if abs(overlap[1] - stop) > 0:
                        cstart, cend = sorted((overlap[1], stop))
                        difference.add(cstart, cend)
            if not intersects:
                difference.add(start, stop)
        return difference
