
"""
provides a classes of simple units that will be used with AxesDivider
class (or others) to determine the size of each axes. The unit
classes define `get_size` method that returns a tuple of two floats,
meaning relative and absolute sizes, respectively.

Note that this class is nothing more than a simple tuple of two
floats. Take a look at the Divider class to see how these two
values are used.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib.cbook as cbook
from matplotlib.axes import Axes

class _Base(object):
    "Base class"

    def __rmul__(self, other):
        float(other) # just to check if number if given
        return Fraction(other, self)

    def __add__(self, other):
        if isinstance(other, _Base):
            return Add(self, other)
        else:
            float(other)
            other = Fixed(other)
            return Add(self, other)


class Add(_Base):
    def __init__(self, a, b):
        self._a = a
        self._b = b

    def get_size(self, renderer):
        a_rel_size, a_abs_size = self._a.get_size(renderer)
        b_rel_size, b_abs_size = self._b.get_size(renderer)
        return a_rel_size + b_rel_size, a_abs_size + b_abs_size

class AddList(_Base):
    def __init__(self, add_list):
        self._list = add_list

    def get_size(self, renderer):
        sum_rel_size = sum([a.get_size(renderer)[0] for a in self._list])
        sum_abs_size = sum([a.get_size(renderer)[1] for a in self._list])
        return sum_rel_size, sum_abs_size


class Fixed(_Base):
    "Simple fixed size  with absolute part = *fixed_size* and relative part = 0"
    def __init__(self, fixed_size):
        self.fixed_size = fixed_size

    def get_size(self, renderer):
        rel_size = 0.
        abs_size = self.fixed_size
        return rel_size, abs_size


class Scaled(_Base):
    "Simple scaled(?) size with absolute part = 0 and relative part = *scalable_size*"
    def __init__(self, scalable_size):
        self._scalable_size = scalable_size

    def get_size(self, renderer):
        rel_size = self._scalable_size
        abs_size = 0.
        return rel_size, abs_size

Scalable=Scaled

def _get_axes_aspect(ax):
    aspect = ax.get_aspect()
    # when aspec is "auto", consider it as 1.
    if aspect in ('normal', 'auto'):
        aspect = 1.
    elif aspect == "equal":
        aspect = 1
    else:
        aspect = float(aspect)

    return aspect

class AxesX(_Base):
    """
    Scaled size whose relative part corresponds to the data width
    of the *axes* multiplied by the *aspect*.
    """
    def __init__(self, axes, aspect=1., ref_ax=None):
        self._axes = axes
        self._aspect = aspect
        if aspect == "axes" and ref_ax is None:
            raise ValueError("ref_ax must be set when aspect='axes'")
        self._ref_ax = ref_ax

    def get_size(self, renderer):
        l1, l2 = self._axes.get_xlim()
        if self._aspect == "axes":
            ref_aspect = _get_axes_aspect(self._ref_ax)
            aspect = ref_aspect/_get_axes_aspect(self._axes)
        else:
            aspect = self._aspect

        rel_size = abs(l2-l1)*aspect
        abs_size = 0.
        return rel_size, abs_size

class AxesY(_Base):
    """
    Scaled size whose relative part corresponds to the data height
    of the *axes* multiplied by the *aspect*.
    """
    def __init__(self, axes, aspect=1., ref_ax=None):
        self._axes = axes
        self._aspect = aspect
        if aspect == "axes" and ref_ax is None:
            raise ValueError("ref_ax must be set when aspect='axes'")
        self._ref_ax = ref_ax

    def get_size(self, renderer):
        l1, l2 = self._axes.get_ylim()

        if self._aspect == "axes":
            ref_aspect = _get_axes_aspect(self._ref_ax)
            aspect = _get_axes_aspect(self._axes)
        else:
            aspect = self._aspect

        rel_size = abs(l2-l1)*aspect
        abs_size = 0.
        return rel_size, abs_size


class MaxExtent(_Base):
    """
    Size whose absolute part is the largest width (or height) of
    the given *artist_list*.
    """
    def __init__(self, artist_list, w_or_h):
        self._artist_list = artist_list

        if w_or_h not in ["width", "height"]:
            raise ValueError()

        self._w_or_h = w_or_h

    def add_artist(self, a):
        self._artist_list.append(a)

    def get_size(self, renderer):
        rel_size = 0.
        w_list, h_list = [], []
        for a in self._artist_list:
            bb = a.get_window_extent(renderer)
            w_list.append(bb.width)
            h_list.append(bb.height)
        dpi = a.get_figure().get_dpi()
        if self._w_or_h == "width":
            abs_size = max(w_list)/dpi
        elif self._w_or_h == "height":
            abs_size = max(h_list)/dpi

        return rel_size, abs_size


class MaxWidth(_Base):
    """
    Size whose absolute part is the largest width of
    the given *artist_list*.
    """
    def __init__(self, artist_list):
        self._artist_list = artist_list

    def add_artist(self, a):
        self._artist_list.append(a)

    def get_size(self, renderer):
        rel_size = 0.
        w_list = []
        for a in self._artist_list:
            bb = a.get_window_extent(renderer)
            w_list.append(bb.width)
        dpi = a.get_figure().get_dpi()
        abs_size = max(w_list)/dpi

        return rel_size, abs_size



class MaxHeight(_Base):
    """
    Size whose absolute part is the largest height of
    the given *artist_list*.
    """
    def __init__(self, artist_list):
        self._artist_list = artist_list

    def add_artist(self, a):
        self._artist_list.append(a)

    def get_size(self, renderer):
        rel_size = 0.
        h_list = []
        for a in self._artist_list:
            bb = a.get_window_extent(renderer)
            h_list.append(bb.height)
        dpi = a.get_figure().get_dpi()
        abs_size = max(h_list)/dpi

        return rel_size, abs_size


class Fraction(_Base):
    """
    An instance whose size is a *fraction* of the *ref_size*.
    ::

      >>> s = Fraction(0.3, AxesX(ax))

    """
    def __init__(self, fraction, ref_size):
        self._fraction_ref = ref_size
        self._fraction = fraction

    def get_size(self, renderer):
        if self._fraction_ref is None:
            return self._fraction, 0.
        else:
            r, a = self._fraction_ref.get_size(renderer)
            rel_size = r*self._fraction
            abs_size = a*self._fraction
            return rel_size, abs_size

class Padded(_Base):
    """
    Return a instance where the absolute part of *size* is
    increase by the amount of *pad*.
    """
    def __init__(self, size, pad):
        self._size = size
        self._pad = pad

    def get_size(self, renderer):
        r, a = self._size.get_size(renderer)
        rel_size = r
        abs_size = a + self._pad
        return rel_size, abs_size

def from_any(size, fraction_ref=None):
    """
    Creates Fixed unit when the first argument is a float, or a
    Fraction unit if that is a string that ends with %. The second
    argument is only meaningful when Fraction unit is created.::

      >>> a = Size.from_any(1.2) # => Size.Fixed(1.2)
      >>> Size.from_any("50%", a) # => Size.Fraction(0.5, a)

    """
    if cbook.is_numlike(size):
        return Fixed(size)
    elif isinstance(size, six.string_types):
        if size[-1] == "%":
            return Fraction(float(size[:-1]) / 100, fraction_ref)

    raise ValueError("Unknown format")


class SizeFromFunc(_Base):
    def __init__(self, func):
        self._func = func

    def get_size(self, renderer):
        rel_size = 0.

        bb = self._func(renderer)
        dpi = renderer.points_to_pixels(72.)
        abs_size = bb/dpi

        return rel_size, abs_size

class GetExtentHelper(object):
    def _get_left(tight_bbox, axes_bbox):
        return axes_bbox.xmin - tight_bbox.xmin

    def _get_right(tight_bbox, axes_bbox):
        return tight_bbox.xmax - axes_bbox.xmax

    def _get_bottom(tight_bbox, axes_bbox):
        return axes_bbox.ymin - tight_bbox.ymin

    def _get_top(tight_bbox, axes_bbox):
        return tight_bbox.ymax - axes_bbox.ymax

    _get_func_map = dict(left=_get_left,
                         right=_get_right,
                         bottom=_get_bottom,
                         top=_get_top)

    del _get_left, _get_right, _get_bottom, _get_top

    def __init__(self, ax, direction):
        if isinstance(ax, Axes):
            self._ax_list = [ax]
        else:
            self._ax_list = ax

        try:
            self._get_func = self._get_func_map[direction]
        except KeyError:
            raise KeyError("direction must be one of left, right, bottom, top")

    def __call__(self, renderer):
        vl = [self._get_func(ax.get_tightbbox(renderer, False),
                             ax.bbox) for ax in self._ax_list]
        return max(vl)
