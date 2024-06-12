from bisect import bisect_left, insort
from functools import total_ordering
from typing import List, Optional, Tuple


@total_ordering
class Interval:
    """
    Class to hold interval data
    """

    def __init__(self, start: Optional[int] = None, end: Optional[int] = None) -> None:
        self._intervals: List[Tuple[int, int]] = []
        if (start is not None) and (end is not None):
            self.add(start, end)

    def __len__(self) -> int:
        return len(self._intervals)

    def __hash__(self) -> int:
        return hash(tuple(self._intervals))

    def __eq__(self, other: "Interval") -> bool:
        if not isinstance(other, Interval):
            return NotImplemented
        return hash(self) == hash(other)

    def __ne__(self, other: "Interval") -> bool:
        if not isinstance(other, Interval):
            return NotImplemented
        return not self == other

    def __lt__(self, other: "Interval") -> bool:
        """__lt__ less than
        Implementation based on portions package interval module
        source: https://github.com/AlexandreDecan/portion/blob/master/portion/interval.py#L549
        Args:
            other: _description_

        Returns:
            _description_
        """
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

    def __and__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        if (self.empty or other.empty) or (self > other) or (self < other):
            return self.__class__()
        intersections: List[Tuple[int, int]] = []
        for sstart, sstop in self.intervals:
            for ostart, ostop in other.intervals:
                if not self._intersects(sstart, sstop, ostart, ostop):
                    continue
                istart: int = max(sstart, ostart)
                istop: int = min(sstop, ostop)
                intersections.append((istart, istop))
        if len(intersections) == 0:
            return self.__class__()
        if len(intersections) == 1:
            return Interval(start=intersections[0][0], end=intersections[0][1])
        return self._merge_overlaps(intersections)

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
            return self._merge_overlaps(self.intervals + other.intervals)

    def __sub__(self, other: "Interval") -> "Interval":
        if not isinstance(other, Interval):
            return NotImplemented
        common = self & other
        if common.empty:
            return self
        difference: List[Tuple[int, int]] = []
        for sstart, sstop in self.intervals:
            overlap = False
            for cstart, cstop in common.intervals:
                if self._intersects(sstart, sstop, cstart, cstop):
                    overlap = True
                    if abs(sstart - cstart) > 0:
                        difference.append(tuple(sorted((sstart, cstart))))  # type: ignore
                    if abs(sstop - cstop) > 0:
                        difference.append(tuple(sorted((sstop, cstop))))  # type: ignore
            if not overlap:
                difference.append((sstart, sstop))
        if len(difference) == 0:
            return self.__class__()
        if len(difference) == 1:
            return Interval(start=difference[0][0], end=difference[0][1])
        return self._merge_overlaps(difference)

    def __invert__(self) -> "Interval":
        if self.empty:
            return self.__class__()
        prev_ends: List[int] = [i[1] for i in self.intervals][:-1]
        cur_starts: List[int] = [i[0] for i in self.intervals][1:]
        inverts: List[Tuple[int, int]] = []
        for end, start in zip(prev_ends, cur_starts):
            if end >= start:
                continue
            inverts.append((end, start))
        if len(inverts) == 0:
            return self.__class__()
        if len(inverts) == 1:
            return Interval(inverts[0][0], inverts[0][1])
        return self._merge_overlaps(inverts)

    @staticmethod
    def _merge_overlaps(overlaps: List[Tuple[int, int]]) -> "Interval":
        """_merge_overlaps merge overlapping intervals

        Args:
            overlaps: List[Tuple[int,int]]

        Returns:
            Interval object
        """
        intervals = sorted(set(overlaps))
        if len(intervals) == 1:
            intv = Interval(start=intervals[0][0], end=intervals[0][1])
            return intv
        iintervals = iter(intervals)
        start, end = next(iintervals)
        merged = []
        for cstart, cend in iintervals:
            if cstart >= end:
                merged.append((start, end))
                start, end = cstart, cend
            else:
                end = max(end, cend)
        merged.append((start, end))
        intv = Interval()
        for start, end in merged:
            intv.add(start=start, end=end)
        return intv

    @property
    def intervals(self) -> List[Tuple[int, int]]:
        if self.empty:
            return []
        return self._intervals

    @property
    def first(self) -> int:
        return self._intervals[0][0]

    @property
    def last(self) -> int:
        return self._intervals[-1][1]

    @property
    def empty(self) -> bool:
        """empty _summary_
        shamelessly stolen from
        https://github.com/AlexandreDecan/portion/blob/master/portion/interval.py#L176
        credits to the author
        Returns:
            _description_
        """
        return len(self._intervals) == 0

    def add(self, start: int, end: int) -> None:
        """add add an interval
        add an interval to the list of intervals
        Args:
            start: int, start position of the new interval
            end: int, end position of the new interval

        Raises:
            ValueError: if either start or end position is < 0
            ValueError: if end == start or end is smaller than start
        """
        if start < 0 or end < 0:
            raise ValueError(
                f"'start' and 'end' values MUST be >= 0. Found {start} and {end}!"
            )
        if end - start <= 0:
            raise ValueError(
                f"'end' MUST be larger than 'start'. Found {start} and  {end}!"
            )
        if len(self._intervals) == 0:
            # no entries so far
            self._intervals.append((start, end))
        elif not self.__contains__(start, end):
            # insert the interval if it is not a duplicate
            self._insert(start, end)

    def __contains__(self, start: int, end: int) -> bool:
        """__contains__ check if an entry exists
        Helper function
        source: https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
        Args:
            start: start position of the interval
            end: end position of the interval

        Returns:
            boolean, True if the interval already exists
        """
        pos = bisect_left(self._intervals, (start, end))
        try:
            if pos < len(self._intervals) and self._intervals[pos] == (start, end):
                return True
        except IndexError:
            return False
        else:
            return False

    def __overlaps__(self, start: int, end: int) -> List[int]:
        """__overlaps__ check for overlaps
        Given start and end positions, return indices of all intervals
        that overlaps with the given positions
        Args:
            start: int, start position of an interval
            end: int, end position of an interl

        Returns:
            List[int]
        """
        # see https://stackoverflow.com/questions/20908047/using-bisect-in-a-list-of-tuples
        ovpos: List[int] = []
        for i, (istart, iend) in enumerate(self._intervals):
            if end <= istart:
                # reached a non overlapping interval
                break
            if self._intersects(start, end, istart, iend):
                # all intervals with which this interval overlaps
                ovpos.append(i)
        return ovpos

    @staticmethod
    def _intersects(start: int, end: int, istart: int, iend: int) -> bool:
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
        return max(start, istart) < min(end, iend)

    def _insert(self, start: int, end: int) -> None:
        """_insert insert new interval or modify exisiting
        Given a non duplicate interval, insert a new interval or modify and extend existing ones
        Args:
            start: int, start position of the new interval
            end: int, end position of then new interval
        """
        if end <= self.first or self.last <= start:
            # simplest case: incoming interval is  either
            # * smaller than the current smallest interval OR
            # * larger than the current largest interval
            insort(self._intervals, (start, end))
        else:
            ovpos = self.__overlaps__(start, end)
            if len(ovpos) == 0:
                # there are no overlaps, insert
                insort(self._intervals, (start, end))
            else:
                new_start, new_end = start, end
                for p in sorted(ovpos, reverse=True):
                    # wow!
                    # source: https://stackoverflow.com/questions/11303225/how-to-remove-multiple-indexes-from-a-list-at-the-same-time
                    new_start = min(new_start, self._intervals[p][0])
                    new_end = max(new_end, self._intervals[p][1])
                    self._intervals.pop(p)
                insort(self._intervals, (new_start, new_end))

    def __repr__(self) -> str:
        if len(self._intervals) == 0:
            return "[]"
        else:
            return "[" + ", ".join([f"({i[0]}, {i[1]})" for i in self._intervals]) + "]"
