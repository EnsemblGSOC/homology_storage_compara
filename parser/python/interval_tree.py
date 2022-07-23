from typing import List, Dict, Tuple, Union, Any
import bisect

class IntervalTree:
    def __init__(self, left, right, mid, intervals):
        self.left = left
        self.right = right
        self.mid = mid
        self.intervals = intervals
        self.intervals_leftsorted = sorted(intervals, key=lambda x: x[0])
        self.intervals_rightsorted = sorted(intervals, key=lambda x: x[1])

    def height(self):
        def height_rec(t):
            if not t:
                return 0
            else:
                return 1 + max(height_rec(t.left), height_rec(t.right))
        return height_rec(self)

    def print_tree(self, width=64):
        height = self.height()
        nodes = [(self, 0)]
        prev_level = 0
        repr_str = ''
        while nodes:
            n, level = nodes.pop(0)
            if prev_level != level:
                prev_level = level
                repr_str += '\n'
            if not n:
                if level < height-1:
                    nodes.extend([(None, level+1), (None, level+1)])
                repr_str += '{val:^{width}}'.format(val='-', width=width//2**level)
            elif n:
                if n.left or level < height-1:
                    nodes.append((n.left, level+1))
                if n.right or level < height-1:
                    nodes.append((n.right, level+1))
                repr_str += '{val:^{width}}'.format(val=str(n.intervals) + ", " + str(n.mid), width=width//2**level)
        print(repr_str)

    def search(self, x) -> List[Tuple[int, int]]:
        if x == self.mid:
            return self.intervals
        if x < self.mid:
            # search in left subtree
            intervals = []
            if self.left:
                intervals += self.left.search(x)
            for i in self.intervals_leftsorted:
                if i[0] <= x:
                    intervals.append(i)
                else:
                    break
            return intervals
        if x > self.mid:
            intervals = []
            if self.right:
                intervals += self.right.search(x)
            for i in self.intervals_rightsorted[::-1]:
                if i[1] >= x:
                    intervals.append(i)
                else:
                    break
            return intervals


def construct_interval_tree(intervals: List[Tuple]):
    left, right = min([x[0] for x in intervals]), max([x[1] for x in intervals])
    mid = (left + right) // 2
    overlapping_intervals = [x for x in intervals if (x[0] <= mid and x[1] >= mid)]
    left_intervals = [x for x in intervals if x[1] < mid]
    right_intervals = [x for x in intervals if x[0] > mid]
    left_subtree = None
    right_subtree = None
    if left_intervals:
        left_subtree = construct_interval_tree(left_intervals)
    if right_intervals:
        right_subtree = construct_interval_tree(right_intervals)
    return IntervalTree(left_subtree, right_subtree, mid, overlapping_intervals)


def exclude_intervals(interval, excluded: List[Tuple]) -> List[Tuple[int, int]]:
    if not excluded:
        return [interval]
    # sort the excluded intervals
    excluded = sorted(excluded)
    left, right = interval
    # discard intervals that do not overlap with the input interval
    start_idx = bisect.bisect_left([b[0] for b in excluded], left)
    end_idx = bisect.bisect_left([b[0] for b in excluded], right)
    excluded = excluded[start_idx:end_idx]
    new_intervals = []
    # trim the starting and ending points of the excluded intervals
    if excluded[0][0] < left:
        excluded[0][0] = left
    if excluded[-1][1] > right:
        excluded[-1][1] = right
    # include the first gap, if any
    if left < excluded[0][0]:
        new_intervals.append((left, excluded[0][0] - 1))
    end = excluded[0][1]
    for x in excluded[1:]:
        # there is a gap between the current position and the lhs
        # of the next interval
        if end < x[0]:
            new_intervals.append((end + 1, x[0] - 1))
            end = x[1]
        elif end < x[1]:
            end = x[1]
    # include the last gap, if any
    if right > excluded[-1][1]:
        new_intervals.append((excluded[-1][1] + 1, right))
    return new_intervals
