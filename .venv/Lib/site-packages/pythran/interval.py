""" Module with facilities to represent range values. """

from math import isinf, isnan
import itertools

import numpy


class Interval(object):

    """ Representation for a range of values. """

    def __init__(self, low, high):
        """ Set initial bound of the range object. """
        if isnan(low):
            low = -float('inf')
        if isnan(high):
            high = +float('inf')
        self._low = low
        self._high = high

    @property
    def low(self):
        return self._low

    @property
    def high(self):
        return self._high

    def __repr__(self):
        """ Return a nicely formatted representation string. """
        return "Interval(low={low}, high={high})".format(low=self.low,
                                                         high=self.high)

    def bounds(self):
        return self.low, self.high

    def __contains__(self, value):
        return self.low <= value <= self.high

    def union(self, other):
        """ Intersect current range with other."""
        return Interval(min(self.low, other.low), max(self.high, other.high))

    def intersect(self, other):
        return Interval(max(self.low, other.low), min(self.high, other.high))

    def copy(self):
        return Interval(self.low, self.high)

    def widen(self, other):
        """ Widen current range. """
        if self.low < other.low:
            low = -float("inf")
        else:
            low = self.low
        if self.high > other.high:
            high = float("inf")
        else:
            high = self.high
        return Interval(low, high)

    def __mul__(self, other):
        """
        Combiner for Multiplication operation.

        >>> Interval(1, 5) * Interval(-5, -4)
        Interval(low=-25, high=-4)
        >>> Interval(-1, 5) * Interval(-5, 3)
        Interval(low=-25, high=15)
        >>> Interval(1, 5) * Interval(3, 8)
        Interval(low=3, high=40)
        """

        def all_bounds():
            return itertools.chain(self.bounds(), other.bounds())

        if any(map(isinf, all_bounds())) and any(x == 0 for x in all_bounds()):
            return UNKNOWN_RANGE

        res = [v1 * v2 for v1, v2 in
               itertools.product(self.bounds(), other.bounds())]
        return Interval(min(res), max(res))

    __mult__ = __mul__

    def __div__(self, other):
        """
        Combiner for Divide operation.

        >>> Interval(-1, 5) / Interval(3, 8)
        Interval(low=-1, high=1)
        >>> Interval(-1, 5) / Interval(-5, -4)
        Interval(low=-2, high=0)
        >>> Interval(-1, 5) / Interval(-5, 3)
        Interval(low=-inf, high=inf)
        """

        if other.low <= 0 and other.high >= 0:
            return UNKNOWN_RANGE
        if other.low == 0:
            return UNKNOWN_RANGE

        def all_bounds():
            return itertools.chain(self.bounds(), other.bounds())

        if any(isinf(x) for x in all_bounds()):
            return UNKNOWN_RANGE

        res = [v1 // v2 for v1, v2 in
               itertools.product(self.bounds(), other.bounds())]
        return Interval(min(res), max(res))

    __truediv__ = __div__

    def __add__(self, other):
        """
        Combiner for Addition operation.

        >>> Interval(-12, 5) + Interval(-5, -3)
        Interval(low=-17, high=2)
        """
        if isinstance(other, IntervalTuple):
            return UNKNOWN_RANGE

        sl, sh, ol, oh = self.low, self.high, other.low, other.high

        if isinf(sl) and isinf(ol) and sl * ol < 0:
            return UNKNOWN_RANGE
        if isinf(sh) and isinf(oh) and sh * oh < 0:
            return UNKNOWN_RANGE

        return Interval(sl + ol, sh + oh)

    def __sub__(self, other):
        """
        Combiner for Subtraction operation.

        >>> Interval(1, 5) - Interval(-5, -4)
        Interval(low=5, high=10)
        """
        sl, sh, ol, oh = self.low, self.high, other.low, other.high

        if isinf(sl) and isinf(oh):
            return UNKNOWN_RANGE
        if isinf(sh) and isinf(ol):
            return UNKNOWN_RANGE

        return Interval(sl - oh, sh - ol)

    def __rshift__(range1, range2):
        """
        Combiner for Right shift operation.

        >>> Interval(10, 100) >> Interval(3, 8)
        Interval(low=0, high=12)
        >>> Interval(10, float("inf")) >> Interval(3, 8)
        Interval(low=0, high=inf)
        >>> Interval(-float("inf"), 0) >> Interval(3, 8)
        Interval(low=-inf, high=0)
        >>> Interval(-30, 10) >> Interval(3, float('inf'))
        Interval(low=-4, high=1)
        """
        if range1.low <= 0:
            if isinf(range1.low):
                min_ = range1.low
            else:
                min_ = range1.low >> range2.low
        elif isinf(range2.high):
            min_ = 0
        else:
            min_ = range1.low >> range2.high
        if isinf(range1.high):
            max_ = range1.high
        elif isinf(range2.low):
            max_ = 0
        else:
            max_ = range1.high >> range2.low
        return Interval(min_, max_)

    def __mod__(range1, range2):
        """ Combiner for Modulo operation.

        >>> Interval(-1, 5) % Interval(1, 13)
        Interval(low=0, high=5)
        >>> Interval(-21, 5) % Interval(1, 13)
        Interval(low=0, high=13)
        """
        return Interval(0, min(range2.high,
                               max(abs(range1.high), abs(range1.low))))

    def __pow__(range1, range2):
        """
        Combiner for Power operation.

        >>> Interval(1, 5) ** Interval(-5, -4)
        Interval(low=1.0, high=1.0)
        >>> Interval(-1, 5) ** Interval(-5, 3)
        Interval(low=-1.0, high=125)
        >>> Interval(1, 5) ** Interval(3, 8)
        Interval(low=1, high=390625)
        """
        res = [v1 ** v2 for v1, v2 in
               itertools.product(range1.bounds(), range2.bounds())]
        minres, maxres = min(res), max(res)
        return Interval(type(minres)(numpy.ceil(minres)),
                        type(maxres)(numpy.floor(maxres)))

    def __lshift__(range1, range2):
        """
        Combiner for Left shift operation.

        >>> Interval(1, 5) << Interval(3, 8)
        Interval(low=8, high=1280)
        >>> Interval(1, float("inf")) << Interval(3, 8)
        Interval(low=8, high=inf)
        >>> Interval(-float("inf"), 0) << Interval(3, 8)
        Interval(low=-inf, high=0)
        >>> Interval(-3, 1) << Interval(3, float('inf'))
        Interval(low=-24, high=inf)
        """
        min_inf = isinf(range1.low) or isinf(range2.low)
        max_inf = isinf(range1.high) or isinf(range2.high)
        min_ = -float("inf") if min_inf else (range1.low << range2.low)
        max_ = float("inf") if max_inf else (range1.high << range2.high)
        return Interval(min_, max_)

    def __floordiv__(range1, range2):
        """
        Combiner for Floor divide operation.

        >>> Interval(-1, 5) // Interval(3, 8)
        Interval(low=-1, high=1)
        >>> Interval(-1, 5) // Interval(-5, -4)
        Interval(low=-2, high=0)
        >>> Interval(-1, 5) // Interval(-5, 3)
        Interval(low=-inf, high=inf)
        """
        if range2.low <= 0 and range2.high >= 0:
            return UNKNOWN_RANGE
        if range2.low == 0:
            return UNKNOWN_RANGE
        res = [v1 if isinf(v1) else (v1 // v2) for v1, v2 in
               itertools.product(range1.bounds(), range2.bounds())]
        return Interval(min(res), max(res))

    def __lt__(self, other):
        """
        Combiner for lower than operation.

        >>> Interval(-1, 5) < Interval(6, 7)
        Interval(low=1, high=1)
        >>> Interval(-1, 5) < Interval(5, 7)
        Interval(low=0, high=1)

        >>> Interval(-1, 5) < Interval(-16, -7)
        Interval(low=0, high=0)

        >>> Interval(1, 5) < Interval(3, 7)
        Interval(low=0, high=1)

        """
        if self.high < other.low:
            return Interval(1, 1)
        if self.low >= other.high:
            return Interval(0, 0)
        return Interval(0, 1)

    def __le__(self, other):
        """
        Combiner for lower than or equal operation.

        >>> Interval(-1, 5) <= Interval(6, 7)
        Interval(low=1, high=1)
        >>> Interval(-1, 5) <= Interval(5, 7)
        Interval(low=1, high=1)

        >>> Interval(-1, 5) <= Interval(-16, -7)
        Interval(low=0, high=0)

        >>> Interval(1, 5) <= Interval(3, 7)
        Interval(low=0, high=1)

        """
        if self.high <= other.low:
            return Interval(1, 1)
        if self.low > other.high:
            return Interval(0, 0)
        return Interval(0, 1)

    def __gt__(self, other):
        """
        Combiner for greater than operation.

        >>> Interval(-5, 1) > Interval(-7, -6)
        Interval(low=1, high=1)
        >>> Interval(-5, 1) > Interval(-7, -5)
        Interval(low=0, high=1)

        >>> Interval(-1, 5) > Interval(6, 7)
        Interval(low=0, high=0)

        >>> Interval(1, 5) > Interval(3, 7)
        Interval(low=0, high=1)

        """
        if self.low > other.high:
            return Interval(1, 1)
        if self.high <= other.low:
            return Interval(0, 0)
        return Interval(0, 1)

    def __ge__(self, other):
        """
        Combiner for greater than or equal operation.

        >>> Interval(-5, 1) >= Interval(-7, -6)
        Interval(low=1, high=1)
        >>> Interval(-5, 1) >= Interval(-7, -5)
        Interval(low=1, high=1)

        >>> Interval(-1, 5) >= Interval(6, 7)
        Interval(low=0, high=0)

        >>> Interval(1, 5) >= Interval(3, 7)
        Interval(low=0, high=1)

        """
        if self.low >= other.high:
            return Interval(1, 1)
        if self.high < other.low:
            return Interval(0, 0)
        return Interval(0, 1)

    def __eq__(self, other):
        """
        Combiner for equal operation.

        >>> Interval(-5, 1) == Interval(-7, -6)
        Interval(low=0, high=0)
        >>> Interval(-5, 1) == Interval(-7, -5)
        Interval(low=0, high=1)
        >>> Interval(-1, 5) == Interval(6, 7)
        Interval(low=0, high=0)
        """
        if isinf(self.low):
            return Interval(0, 1)
        elif self.low == self.high == other.low == other.high:
            return Interval(1, 1)
        elif (self < other) or (self > other):
            return Interval(0, 0)
        else:
            return Interval(0, 1)

    def __ne__(self, other):
        """
        Combiner for not equal operation.

        >>> Interval(-5, 1) != Interval(-7, -6)
        Interval(low=1, high=1)
        >>> Interval(-5, 1) != Interval(-7, -5)
        Interval(low=0, high=1)
        >>> Interval(-1, 5) != Interval(6, 7)
        Interval(low=1, high=1)
        """
        if isinf(self.low):
            return Interval(0, 1)
        elif self.low == self.high == other.low == other.high:
            return Interval(1, 1)
        elif (self < other) or (self > other):
            return Interval(1, 1)
        else:
            return Interval(0, 1)

    def __nonzero__(self):
        return not isinf(self.high) and self.low == self.high and self .low > 0

    def __getitem__(self, index):
        return UNKNOWN_RANGE

    __bool__ = __nonzero__


class IntervalTuple(object):

    def __init__(self, values):
        self.values = tuple(values)

    def union(self, other):
        if isinstance(other, Interval):
            return UNKNOWN_TUPLE_RANGE
        return IntervalTuple(x.union(y) for x, y in zip(self.values,
                                                        other.values))

    def intersect(self, other):
        if isinstance(other, Interval):
            return UNKNOWN_TUPLE_RANGE
        return IntervalTuple(x.intersect(y) for x, y in zip(self.values,
                                                            other.values))

    @property
    def high(self):
        return UNKNOWN_RANGE.high

    @property
    def low(self):
        return UNKNOWN_RANGE.low

    def __getitem__(self, index):
        out = None
        low = max(0, index.low)
        high = min(len(self.values) - 1, index.high)
        for i in range(low, high + 1):
            if out is None:
                out = self.values[i]
            else:
                out = out.union(self.values[i])
        return out or UNKNOWN_RANGE

    def widen(self, other):
        if isinstance(other, Interval):
            return UNKNOWN_TUPLE_RANGE
        return IntervalTuple(s.widen(o) for s, o in zip(self.values,
                                                        other.values))

    def __add__(self, other):
        if isinstance(other, Interval):
            return UNKNOWN_TUPLE_RANGE
        return IntervalTuple(self.values + other.values)


UNKNOWN_RANGE = Interval(-float("inf"), float("inf"))
UNKNOWN_TUPLE_RANGE = IntervalTuple([UNKNOWN_RANGE])


def range_values(args):
    """ Function used to compute returned range value of [x]range function. """
    if len(args) == 1:
        return Interval(0, args[0].high)
    elif len(args) == 2:
        return Interval(args[0].low, args[1].high)
    elif len(args) == 3:
        is_neg = args[2].low < 0
        is_pos = args[2].high > 0
        if is_neg and is_pos:
            return UNKNOWN_RANGE
        elif is_neg:
            return Interval(args[1].low, args[0].high - 1)
        else:
            return Interval(args[0].low, args[1].high - 1)


def bool_values(_):
    """ Return the range of a boolean value, i.e. [0, 1]. """
    return Interval(0, 1)


def cmp_values(_):
    """ Return the range of a comparison value, i.e. [-1, 1]. """
    return Interval(-1, 1)


def positive_values(_):
    """ Return a positive range without upper bound. """
    return Interval(0, float("inf"))


def max_values(args):
    """ Return possible range for max function. """
    return Interval(max(x.low for x in args), max(x.high for x in args))


def min_values(args):
    """ Return possible range for min function. """
    return Interval(min(x.low for x in args), min(x.high for x in args))


def ord_values(_):
    """ Return possible range for ord function. """
    return Interval(0, 255)
