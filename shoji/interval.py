from typing import List, Optional, Tuple
from bisect import bisect_left, bisect_right
from sortedcontainers import SortedList
from operator import itemgetter


class Interval(SortedList):
    def __init__(self):
        super().__init__()

    def __repr__(self) -> str:
        if len(self) == 0:
            return "[]"
        else:
            return "[" + ", ".join([f"({i[0]}, {i[1]})" for i in self]) + "]"

    def __iter__(self):
        return super().__iter__()

    def __len__(self):
        return super().__len__()

    def __getitem__(self, index):
        return super().__getitem__(index)

    def __reversed__(self):
        return super().__reversed__()

    def __delitem__(self, index):
        super().__delitem__(index)

    def __hash__(self):
        return hash(tuple(self))

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
        return self

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

    @property
    def first(self) -> int:
        return self[0][0]

    @property
    def last(self) -> int:
        return self[-1][1]

    @property
    def starts(self) -> List[int]:
        return [i[0] for i in self]

    @property
    def ends(self) -> List[int]:
        return [i[1] for i in self]

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

    def add(self, start: int, end: int) -> None:
        if start < 0 or end < 0:
            raise ValueError(
                f"'start' and 'end' values MUST be >= 0. Found {start} and {end}!"
            )
        if end - start <= 0:
            raise ValueError(
                f"'end' MUST be larger than 'start'. Found {start} and  {end}!"
            )
        # find the position where the existing end is >= the start
        start_pos = bisect_right(self, start, key=itemgetter(1))
        # find the position where the existing start is <= the end
        # no need to search from 0th index, start from start_pos
        end_pos = bisect_left(self, end, lo=start_pos, key=itemgetter(0))
        if start_pos == end_pos:
            # new interval does not overlap any existing intervals
            # or the list is empty
            super().add((start, end))
        else:
            # end_pos will always be the "next" index where this end can be inserted,
            # but this does not mean that this end overlaps with the interval positions at end_pos
            for i in range(end_pos - 1, start_pos - 1, -1):
                if self[i] == (start, end):
                    # duplicate interval
                    return
                # new interval overlaps with this/these existing interval(s)
                start = min(start, self[i][0])
                end = max(end, self[i][1])
                _ = self.pop(i)
            super().add((start, end))

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
            self = merged

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

    def islice(
        self, istart: Optional[int] = None, istop: Optional[int] = None
    ) -> "Interval":
        """islice _summary_
        override the SortedList.islice method to return an Interval object
        Args:
            istart: start index. Defaults to None.
            istop: end index. Defaults to None.

        Returns:
            Interval object
        """
        sliced: Interval = Interval()
        for start, end in super().islice(start=istart, stop=istop, reverse=False):
            sliced.add(start, end)
        return sliced

    def get_overlaps(self, start: int, end: int) -> "Interval":
        """get_overlaps
        Given a start and end position, return the overlapping intervals
        Args:
            start: int, interval start
            end: int, interval end

        Returns:
            Interval object
        """
        start_pos = bisect_right(self, start, key=itemgetter(1))
        end_pos = bisect_left(self, end, lo=start_pos, key=itemgetter(0))
        return self.islice(start_pos, end_pos)

    def __and__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        if (self.empty or other.empty) or (self > other) or (self < other):
            return self.__class()
        intersections: Interval = Interval()
        for ostart, oend in other:
            start_pos = bisect_right(self, ostart, key=itemgetter(1))
            end_pos = bisect_left(self, oend, lo=start_pos, key=itemgetter(0))
            if start_pos == end_pos:
                # this interval fits between two intervals in self
                # or is either smaller or larger than all intervals in self
                continue
            for i in range(start_pos, end_pos):
                ostart = max(ostart, self[i][0])
                oend = min(oend, self[i][1])
            if oend > ostart:
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
            for start, end in self:
                unions.add(start, end)
            for ostart, oend in other:
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
        difference: Interval = Interval()
        for start, end in self:
            start_pos = bisect_right(other, start, key=itemgetter(1))
            end_pos = bisect_left(other, end, lo=start_pos, key=itemgetter(0))
            if start_pos == end_pos:
                # no overlaps
                difference.add(start, end)
                continue
            overlaps: List[Tuple[int, int]] = list(other.islice(start_pos, end_pos))
            if start < overlaps[0][0]:
                difference.add(start, overlaps[0][0])
            if end > overlaps[-1][1]:
                difference.add(overlaps[-1][1], end)
            if len(overlaps) > 1:
                prev_ends: List[int] = [i[1] for i in overlaps][:-1]
                next_starts: List[int] = [i[0] for i in overlaps][1:]
                for pend, nstart in zip(prev_ends, next_starts):
                    if pend < nstart:
                        difference.add(pend, nstart)
        return difference
