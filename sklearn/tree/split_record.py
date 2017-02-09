from __future__ import division, print_function

from numpy import nan
from numpy import inf

from .stats_node import StatsNode
from .criterion import _impurity_mse


class SplitRecord(object):
    """Container of the parent and children statistics at a node

    Attributes
    ----------
    feature: int,
        The index of the feature to apply the split.

    pos: int,
        The sample index where to apply the split.

    threshold: float,
        The threshold of X[pos, feature] to make the split.

    impurity: float,
        The impurity for that specific split.

    nid: int,
        The id of the current node.

    c_stats: StatsNode,
        Current node statistics.

    l_stats: StatsNode,
        Left node statistics.

    r_stats: StatsNode,
        Right node statistics.
    """

    def __init__(self):
        self._feature = 0
        self._pos = 0
        self._threshold = nan
        self._impurity = inf
        self._nid = 0

        # statistics related to current and children node
        self._c_stats = StatsNode(0., 0., 0, 0.)
        self._l_stats = StatsNode(0., 0., 0, 0.)
        self._r_stats = StatsNode(0., 0., 0, 0.)

    def reset(self, feature, pos, threshold, impurity,
              nid, c_stats, l_stats, r_stats):
        """Reset the split record"""
        self._feature = int(feature)
        self._pos = (pos)
        self._threshold = float(threshold)
        self._impurity = float(impurity)
        self._nid = int(nid)
        self._c_stats = c_stats
        self._l_stats = l_stats
        self._r_stats = r_stats

    def clear(self):
        """Clear the split record"""
        self._feature = 0
        self._pos = 0
        self._threshold = nan
        self._impurity = inf
        self._nid = 0
        self._c_stats.clear()
        self._l_stats.clear()
        self._r_stats.clear()

    @property
    def feature(self):
        return self._feature

    @feature.setter
    def feature(self, val):
        self._feature = int(val)

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, val):
        self._pos = int(val)

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, val):
        self._threshold = float(val)

    @property
    def impurity(self):
        return self._impurity

    @impurity.setter
    def impurity(self, val):
        self._impurity = float(val)

    @property
    def nid(self):
        return self._nid

    @nid.setter
    def nid(self, val):
        self._nid = int(val)

    @property
    def c_stats(self):
        return self._c_stats

    @c_stats.setter
    def c_stats(self, val):
        self._c_stats = val

    @property
    def l_stats(self):
        return self._l_stats

    @l_stats.setter
    def l_stats(self, val):
        self._l_stats = val

    @property
    def r_stats(self):
        return self._r_stats

    @r_stats.setter
    def r_stats(self, val):
        self._r_stats = val

    def expand_record(self):
        """Create two new records from the left and right stats"""
        # create the left child split record
        left_sr = SplitRecord()
        left_sr.c_stats = self._l_stats
        # FIXME stuck with impurity mse for the moment
        left_sr.impurity = _impurity_mse(left_sr.c_stats)

        # create the right child split record
        right_sr = SplitRecord()
        right_sr.c_stats = self._l_stats
        # FIXME stuck with impurity mse for the moment
        right_sr.impurity = _impurity_mse(right_sr.c_stats)

        return left_sr, right_sr

    def __str__(self):
        info = ("feature: {}\n"
                "position: {}\n"
                "threshold: {}\n"
                "impurity: {}\n"
                "node id: {}\n"
                "current stats: {}\n"
                "left stats: {}\n"
                "right stats: {}\n".format(self._feature,
                                           self._pos,
                                           self._threshold,
                                           self._impurity,
                                           self._nid,
                                           self.c_stats,
                                           self.l_stats,
                                           self.r_stats))
        return info
