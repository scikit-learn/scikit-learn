"""
A module for converting numbers or color arguments to *RGB* or *RGBA*

*RGB* and *RGBA* are sequences of, respectively, 3 or 4 floats in the
range 0-1.

This module includes functions and classes for color specification
conversions, and for mapping numbers to colors in a 1-D array of colors called
a colormap. Colormapping typically involves two steps: a data array is first
mapped onto the range 0-1 using an instance of :class:`Normalize` or of a
subclass; then this number in the 0-1 range is mapped to a color using an
instance of a subclass of :class:`Colormap`.  Two are provided here:
:class:`LinearSegmentedColormap`, which is used to generate all the built-in
colormap instances, but is also useful for making custom colormaps, and
:class:`ListedColormap`, which is used for generating a custom colormap from a
list of color specifications.

The module also provides functions for checking whether an object can be
interpreted as a color (:func:`is_color_like`), for converting such an object
to an RGBA tuple (:func:`to_rgba`) or to an HTML-like hex string in the
`#rrggbb` format (:func:`to_hex`), and a sequence of colors to an `(n, 4)`
RGBA array (:func:`to_rgba_array`).  Caching is used for efficiency.

Matplotlib recognizes the following formats to specify a color:

* an RGB or RGBA tuple of float values in ``[0, 1]`` (e.g., ``(0.1, 0.2, 0.5)``
  or  ``(0.1, 0.2, 0.5, 0.3)``);
* a hex RGB or RGBA string (e.g., ``'#0F0F0F'`` or ``'#0F0F0F0F'``);
* a string representation of a float value in ``[0, 1]`` inclusive for gray
  level (e.g., ``'0.5'``);
* one of ``{'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}``;
* a X11/CSS4 color name;
* a name from the `xkcd color survey <https://xkcd.com/color/rgb/>`__;
  prefixed with ``'xkcd:'`` (e.g., ``'xkcd:sky blue'``);
* one of ``{'tab:blue', 'tab:orange', 'tab:green',
  'tab:red', 'tab:purple', 'tab:brown', 'tab:pink',
  'tab:gray', 'tab:olive', 'tab:cyan'}`` which are the Tableau Colors from the
  'T10' categorical palette (which is the default color cycle);
* a "CN" color spec, i.e. `'C'` followed by a single digit, which is an index
  into the default property cycle (``matplotlib.rcParams['axes.prop_cycle']``);
  the indexing occurs at artist creation time and defaults to black if the
  cycle does not include color.

All string specifications of color, other than "CN", are case-insensitive.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip

from collections import Sized
import itertools
import re
import warnings

import numpy as np
import matplotlib.cbook as cbook
from ._color_data import BASE_COLORS, TABLEAU_COLORS, CSS4_COLORS, XKCD_COLORS


class _ColorMapping(dict):
    def __init__(self, mapping):
        super(_ColorMapping, self).__init__(mapping)
        self.cache = {}

    def __setitem__(self, key, value):
        super(_ColorMapping, self).__setitem__(key, value)
        self.cache.clear()

    def __delitem__(self, key):
        super(_ColorMapping, self).__delitem__(key)
        self.cache.clear()


_colors_full_map = {}
# Set by reverse priority order.
_colors_full_map.update(XKCD_COLORS)
_colors_full_map.update({k.replace('grey', 'gray'): v
                         for k, v in XKCD_COLORS.items()
                         if 'grey' in k})
_colors_full_map.update(CSS4_COLORS)
_colors_full_map.update(TABLEAU_COLORS)
_colors_full_map.update({k.replace('gray', 'grey'): v
                         for k, v in TABLEAU_COLORS.items()
                         if 'gray' in k})
_colors_full_map.update(BASE_COLORS)
_colors_full_map = _ColorMapping(_colors_full_map)


def get_named_colors_mapping():
    """Return the global mapping of names to named colors."""
    return _colors_full_map


def _sanitize_extrema(ex):
    if ex is None:
        return ex
    try:
        ret = np.asscalar(ex)
    except AttributeError:
        ret = float(ex)
    return ret


def _is_nth_color(c):
    """Return whether *c* can be interpreted as an item in the color cycle."""
    return isinstance(c, six.string_types) and re.match(r"\AC[0-9]\Z", c)


def is_color_like(c):
    """Return whether *c* can be interpreted as an RGB(A) color."""
    # Special-case nth color syntax because it cannot be parsed during
    # setup.
    if _is_nth_color(c):
        return True
    try:
        to_rgba(c)
    except ValueError:
        return False
    else:
        return True


def same_color(c1, c2):
    """
    Compare two colors to see if they are the same.

    Parameters
    ----------
    c1, c2 : Matplotlib colors

    Returns
    -------
    bool
        ``True`` if *c1* and *c2* are the same color, otherwise ``False``.
    """
    return (to_rgba_array(c1) == to_rgba_array(c2)).all()


def to_rgba(c, alpha=None):
    """
    Convert *c* to an RGBA color.

    Parameters
    ----------
    c : Matplotlib color

    alpha : scalar, optional
        If *alpha* is not ``None``, it forces the alpha value, except if *c* is
        ``"none"`` (case-insensitive), which always maps to ``(0, 0, 0, 0)``.

    Returns
    -------
    tuple
        Tuple of ``(r, g, b, a)`` scalars.
    """
    # Special-case nth color syntax because it should not be cached.
    if _is_nth_color(c):
        from matplotlib import rcParams
        prop_cycler = rcParams['axes.prop_cycle']
        colors = prop_cycler.by_key().get('color', ['k'])
        c = colors[int(c[1]) % len(colors)]
    try:
        rgba = _colors_full_map.cache[c, alpha]
    except (KeyError, TypeError):  # Not in cache, or unhashable.
        rgba = _to_rgba_no_colorcycle(c, alpha)
        try:
            _colors_full_map.cache[c, alpha] = rgba
        except TypeError:
            pass
    return rgba


def _to_rgba_no_colorcycle(c, alpha=None):
    """Convert *c* to an RGBA color, with no support for color-cycle syntax.

    If *alpha* is not ``None``, it forces the alpha value, except if *c* is
    ``"none"`` (case-insensitive), which always maps to ``(0, 0, 0, 0)``.
    """
    orig_c = c
    if isinstance(c, six.string_types):
        if c.lower() == "none":
            return (0., 0., 0., 0.)
        # Named color.
        try:
            # This may turn c into a non-string, so we check again below.
            c = _colors_full_map[c.lower()]
        except KeyError:
            pass
    if isinstance(c, six.string_types):
        # hex color with no alpha.
        match = re.match(r"\A#[a-fA-F0-9]{6}\Z", c)
        if match:
            return (tuple(int(n, 16) / 255
                          for n in [c[1:3], c[3:5], c[5:7]])
                    + (alpha if alpha is not None else 1.,))
        # hex color with alpha.
        match = re.match(r"\A#[a-fA-F0-9]{8}\Z", c)
        if match:
            color = [int(n, 16) / 255
                     for n in [c[1:3], c[3:5], c[5:7], c[7:9]]]
            if alpha is not None:
                color[-1] = alpha
            return tuple(color)
        # string gray.
        try:
            return (float(c),) * 3 + (alpha if alpha is not None else 1.,)
        except ValueError:
            pass
        raise ValueError("Invalid RGBA argument: {!r}".format(orig_c))
    # tuple color.
    c = np.array(c)
    if not np.can_cast(c.dtype, float, "same_kind") or c.ndim != 1:
        # Test the dtype explicitly as `map(float, ...)`, `np.array(...,
        # float)` and `np.array(...).astype(float)` all convert "0.5" to 0.5.
        # Test dimensionality to reject single floats.
        raise ValueError("Invalid RGBA argument: {!r}".format(orig_c))
    # Return a tuple to prevent the cached value from being modified.
    c = tuple(c.astype(float))
    if len(c) not in [3, 4]:
        raise ValueError("RGBA sequence should have length 3 or 4")
    if len(c) == 3 and alpha is None:
        alpha = 1
    if alpha is not None:
        c = c[:3] + (alpha,)
    if any(elem < 0 or elem > 1 for elem in c):
        raise ValueError("RGBA values should be within 0-1 range")
    return c


def to_rgba_array(c, alpha=None):
    """Convert *c* to a (n, 4) array of RGBA colors.

    If *alpha* is not ``None``, it forces the alpha value.  If *c* is
    ``"none"`` (case-insensitive) or an empty list, an empty array is returned.
    """
    # Special-case inputs that are already arrays, for performance.  (If the
    # array has the wrong kind or shape, raise the error during one-at-a-time
    # conversion.)
    if (isinstance(c, np.ndarray) and c.dtype.kind in "if"
            and c.ndim == 2 and c.shape[1] in [3, 4]):
        if c.shape[1] == 3:
            result = np.column_stack([c, np.zeros(len(c))])
            result[:, -1] = alpha if alpha is not None else 1.
        elif c.shape[1] == 4:
            result = c.copy()
            if alpha is not None:
                result[:, -1] = alpha
        if np.any((result < 0) | (result > 1)):
            raise ValueError("RGBA values should be within 0-1 range")
        return result
    # Handle single values.
    # Note that this occurs *after* handling inputs that are already arrays, as
    # `to_rgba(c, alpha)` (below) is expensive for such inputs, due to the need
    # to format the array in the ValueError message(!).
    if isinstance(c, six.string_types) and c.lower() == "none":
        return np.zeros((0, 4), float)
    try:
        return np.array([to_rgba(c, alpha)], float)
    except (ValueError, TypeError):
        pass
    # Convert one at a time.
    result = np.empty((len(c), 4), float)
    for i, cc in enumerate(c):
        result[i] = to_rgba(cc, alpha)
    return result


def to_rgb(c):
    """Convert *c* to an RGB color, silently dropping the alpha channel."""
    return to_rgba(c)[:3]


def to_hex(c, keep_alpha=False):
    """Convert *c* to a hex color.

    Uses the ``#rrggbb`` format if *keep_alpha* is False (the default),
    ``#rrggbbaa`` otherwise.
    """
    c = to_rgba(c)
    if not keep_alpha:
        c = c[:3]
    return "#" + "".join(format(int(np.round(val * 255)), "02x")
                         for val in c)


### Backwards-compatible color-conversion API


cnames = CSS4_COLORS
hexColorPattern = re.compile(r"\A#[a-fA-F0-9]{6}\Z")
rgb2hex = to_hex
hex2color = to_rgb


class ColorConverter(object):
    """
    Provides methods for converting color specifications to *RGB* or *RGBA*

    Caching is used for more efficient conversion upon repeated calls
    with the same argument.

    Ordinarily only the single instance instantiated in this module,
    *colorConverter*, is needed.
    """

    colors = _colors_full_map
    cache = _colors_full_map.cache

    @staticmethod
    def to_rgb(arg):
        """
        Returns an *RGB* tuple of three floats from 0-1.

        *arg* can be an *RGB* or *RGBA* sequence or a string in any of
        several forms:

            1) a letter from the set 'rgbcmykw'
            2) a hex color string, like '#00FFFF'
            3) a standard name, like 'aqua'
            4) a string representation of a float, like '0.4',
               indicating gray on a 0-1 scale

        if *arg* is *RGBA*, the *A* will simply be discarded.
        """
        return to_rgb(arg)

    @staticmethod
    def to_rgba(arg, alpha=None):
        """
        Returns an *RGBA* tuple of four floats from 0-1.

        For acceptable values of *arg*, see :meth:`to_rgb`.
        In addition, if *arg* is "none" (case-insensitive),
        then (0,0,0,0) will be returned.
        If *arg* is an *RGBA* sequence and *alpha* is not *None*,
        *alpha* will replace the original *A*.
        """
        return to_rgba(arg, alpha)

    @staticmethod
    def to_rgba_array(arg, alpha=None):
        """
        Returns a numpy array of *RGBA* tuples.

        Accepts a single mpl color spec or a sequence of specs.

        Special case to handle "no color": if *c* is "none" (case-insensitive),
        then an empty array will be returned.  Same for an empty list.
        """
        return to_rgba_array(arg, alpha)


colorConverter = ColorConverter()


### End of backwards-compatible color-conversion API


def makeMappingArray(N, data, gamma=1.0):
    """Create an *N* -element 1-d lookup table

    *data* represented by a list of x,y0,y1 mapping correspondences.
    Each element in this list represents how a value between 0 and 1
    (inclusive) represented by x is mapped to a corresponding value
    between 0 and 1 (inclusive). The two values of y are to allow
    for discontinuous mapping functions (say as might be found in a
    sawtooth) where y0 represents the value of y for values of x
    <= to that given, and y1 is the value to be used for x > than
    that given). The list must start with x=0, end with x=1, and
    all values of x must be in increasing order. Values between
    the given mapping points are determined by simple linear interpolation.

    Alternatively, data can be a function mapping values between 0 - 1
    to 0 - 1.

    The function returns an array "result" where ``result[x*(N-1)]``
    gives the closest value for values of x between 0 and 1.
    """

    if callable(data):
        xind = np.linspace(0, 1, N) ** gamma
        lut = np.clip(np.array(data(xind), dtype=float), 0, 1)
        return lut

    try:
        adata = np.array(data)
    except Exception:
        raise TypeError("data must be convertible to an array")
    shape = adata.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError("data must be nx3 format")

    x = adata[:, 0]
    y0 = adata[:, 1]
    y1 = adata[:, 2]

    if x[0] != 0. or x[-1] != 1.0:
        raise ValueError(
            "data mapping points must start with x=0 and end with x=1")
    if (np.diff(x) < 0).any():
        raise ValueError("data mapping points must have x in increasing order")
    # begin generation of lookup table
    x = x * (N - 1)
    lut = np.zeros((N,), float)
    xind = (N - 1) * np.linspace(0, 1, N) ** gamma
    ind = np.searchsorted(x, xind)[1:-1]

    distance = (xind[1:-1] - x[ind - 1]) / (x[ind] - x[ind - 1])
    lut[1:-1] = distance * (y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0] = y1[0]
    lut[-1] = y0[-1]
    # ensure that the lut is confined to values between 0 and 1 by clipping it
    return np.clip(lut, 0.0, 1.0)


class Colormap(object):
    """
    Baseclass for all scalar to RGBA mappings.

    Typically Colormap instances are used to convert data values (floats) from
    the interval ``[0, 1]`` to the RGBA color that the respective Colormap
    represents. For scaling of data into the ``[0, 1]`` interval see
    :class:`matplotlib.colors.Normalize`. It is worth noting that
    :class:`matplotlib.cm.ScalarMappable` subclasses make heavy use of this
    ``data->normalize->map-to-color`` processing chain.

    """
    def __init__(self, name, N=256):
        """
        Parameters
        ----------
        name : str
            The name of the colormap.
        N : int
            The number of rgb quantization levels.

        """
        self.name = name
        self.N = int(N)  # ensure that N is always int
        self._rgba_bad = (0.0, 0.0, 0.0, 0.0)  # If bad, don't paint anything.
        self._rgba_under = None
        self._rgba_over = None
        self._i_under = self.N
        self._i_over = self.N + 1
        self._i_bad = self.N + 2
        self._isinit = False

        #: When this colormap exists on a scalar mappable and colorbar_extend
        #: is not False, colorbar creation will pick up ``colorbar_extend`` as
        #: the default value for the ``extend`` keyword in the
        #: :class:`matplotlib.colorbar.Colorbar` constructor.
        self.colorbar_extend = False

    def __call__(self, X, alpha=None, bytes=False):
        """
        Parameters
        ----------
        X : scalar, ndarray
            The data value(s) to convert to RGBA.
            For floats, X should be in the interval ``[0.0, 1.0]`` to
            return the RGBA values ``X*100`` percent along the Colormap line.
            For integers, X should be in the interval ``[0, Colormap.N)`` to
            return RGBA values *indexed* from the Colormap with index ``X``.
        alpha : float, None
            Alpha must be a scalar between 0 and 1, or None.
        bytes : bool
            If False (default), the returned RGBA values will be floats in the
            interval ``[0, 1]`` otherwise they will be uint8s in the interval
            ``[0, 255]``.

        Returns
        -------
        Tuple of RGBA values if X is scalar, otherwise an array of
        RGBA values with a shape of ``X.shape + (4, )``.

        """
        # See class docstring for arg/kwarg documentation.
        if not self._isinit:
            self._init()
        mask_bad = None
        if not cbook.iterable(X):
            vtype = 'scalar'
            xa = np.array([X])
        else:
            vtype = 'array'
            xma = np.ma.array(X, copy=True)  # Copy here to avoid side effects.
            mask_bad = xma.mask              # Mask will be used below.
            xa = xma.filled()                # Fill to avoid infs, etc.
            del xma

        # Calculations with native byteorder are faster, and avoid a
        # bug that otherwise can occur with putmask when the last
        # argument is a numpy scalar.
        if not xa.dtype.isnative:
            xa = xa.byteswap().newbyteorder()

        if xa.dtype.kind == "f":
            xa *= self.N
            # Negative values are out of range, but astype(int) would truncate
            # them towards zero.
            xa[xa < 0] = -1
            # xa == 1 (== N after multiplication) is not out of range.
            xa[xa == self.N] = self.N - 1
            # Avoid converting large positive values to negative integers.
            np.clip(xa, -1, self.N, out=xa)
            xa = xa.astype(int)
        # Set the over-range indices before the under-range;
        # otherwise the under-range values get converted to over-range.
        xa[xa > self.N - 1] = self._i_over
        xa[xa < 0] = self._i_under
        if mask_bad is not None:
            if mask_bad.shape == xa.shape:
                np.copyto(xa, self._i_bad, where=mask_bad)
            elif mask_bad:
                xa.fill(self._i_bad)
        if bytes:
            lut = (self._lut * 255).astype(np.uint8)
        else:
            lut = self._lut.copy()  # Don't let alpha modify original _lut.

        if alpha is not None:
            alpha = min(alpha, 1.0)  # alpha must be between 0 and 1
            alpha = max(alpha, 0.0)
            if bytes:
                alpha = int(alpha * 255)
            if (lut[-1] == 0).all():
                lut[:-1, -1] = alpha
                # All zeros is taken as a flag for the default bad
                # color, which is no color--fully transparent.  We
                # don't want to override this.
            else:
                lut[:, -1] = alpha
                # If the bad value is set to have a color, then we
                # override its alpha just as for any other value.

        rgba = np.empty(shape=xa.shape + (4,), dtype=lut.dtype)
        lut.take(xa, axis=0, mode='clip', out=rgba)
        if vtype == 'scalar':
            rgba = tuple(rgba[0, :])
        return rgba

    def __copy__(self):
        """Create new object with the same class, update attributes
        """
        cls = self.__class__
        cmapobject = cls.__new__(cls)
        cmapobject.__dict__.update(self.__dict__)
        if self._isinit:
            cmapobject._lut = np.copy(self._lut)
        return cmapobject

    def set_bad(self, color='k', alpha=None):
        """Set color to be used for masked values.
        """
        self._rgba_bad = colorConverter.to_rgba(color, alpha)
        if self._isinit:
            self._set_extremes()

    def set_under(self, color='k', alpha=None):
        """Set color to be used for low out-of-range values.
           Requires norm.clip = False
        """
        self._rgba_under = colorConverter.to_rgba(color, alpha)
        if self._isinit:
            self._set_extremes()

    def set_over(self, color='k', alpha=None):
        """Set color to be used for high out-of-range values.
           Requires norm.clip = False
        """
        self._rgba_over = colorConverter.to_rgba(color, alpha)
        if self._isinit:
            self._set_extremes()

    def _set_extremes(self):
        if self._rgba_under:
            self._lut[self._i_under] = self._rgba_under
        else:
            self._lut[self._i_under] = self._lut[0]
        if self._rgba_over:
            self._lut[self._i_over] = self._rgba_over
        else:
            self._lut[self._i_over] = self._lut[self.N - 1]
        self._lut[self._i_bad] = self._rgba_bad

    def _init(self):
        """Generate the lookup table, self._lut"""
        raise NotImplementedError("Abstract class only")

    def is_gray(self):
        if not self._isinit:
            self._init()
        return (np.all(self._lut[:, 0] == self._lut[:, 1]) and
                np.all(self._lut[:, 0] == self._lut[:, 2]))

    def _resample(self, lutsize):
        """
        Return a new color map with *lutsize* entries.
        """
        raise NotImplementedError()

    def reversed(self, name=None):
        """
        Make a reversed instance of the Colormap.

        .. note :: Function not implemented for base class.

        Parameters
        ----------
        name : str, optional
            The name for the reversed colormap. If it's None the
            name will be the name of the parent colormap + "_r".

        Notes
        -----
        See :meth:`LinearSegmentedColormap.reversed` and
        :meth:`ListedColormap.reversed`
        """
        raise NotImplementedError()


class LinearSegmentedColormap(Colormap):
    """Colormap objects based on lookup tables using linear segments.

    The lookup table is generated using linear interpolation for each
    primary color, with the 0-1 domain divided into any number of
    segments.
    """
    def __init__(self, name, segmentdata, N=256, gamma=1.0):
        """Create color map from linear mapping segments

        segmentdata argument is a dictionary with a red, green and blue
        entries. Each entry should be a list of *x*, *y0*, *y1* tuples,
        forming rows in a table. Entries for alpha are optional.

        Example: suppose you want red to increase from 0 to 1 over
        the bottom half, green to do the same over the middle half,
        and blue over the top half.  Then you would use::

            cdict = {'red':   [(0.0,  0.0, 0.0),
                               (0.5,  1.0, 1.0),
                               (1.0,  1.0, 1.0)],

                     'green': [(0.0,  0.0, 0.0),
                               (0.25, 0.0, 0.0),
                               (0.75, 1.0, 1.0),
                               (1.0,  1.0, 1.0)],

                     'blue':  [(0.0,  0.0, 0.0),
                               (0.5,  0.0, 0.0),
                               (1.0,  1.0, 1.0)]}

        Each row in the table for a given color is a sequence of
        *x*, *y0*, *y1* tuples.  In each sequence, *x* must increase
        monotonically from 0 to 1.  For any input value *z* falling
        between *x[i]* and *x[i+1]*, the output value of a given color
        will be linearly interpolated between *y1[i]* and *y0[i+1]*::

            row i:   x  y0  y1
                           /
                          /
            row i+1: x  y0  y1

        Hence y0 in the first row and y1 in the last row are never used.


        .. seealso::

               :meth:`LinearSegmentedColormap.from_list`
               Static method; factory function for generating a
               smoothly-varying LinearSegmentedColormap.

               :func:`makeMappingArray`
               For information about making a mapping array.
        """
        # True only if all colors in map are identical; needed for contouring.
        self.monochrome = False
        Colormap.__init__(self, name, N)
        self._segmentdata = segmentdata
        self._gamma = gamma

    def _init(self):
        self._lut = np.ones((self.N + 3, 4), float)
        self._lut[:-3, 0] = makeMappingArray(
            self.N, self._segmentdata['red'], self._gamma)
        self._lut[:-3, 1] = makeMappingArray(
            self.N, self._segmentdata['green'], self._gamma)
        self._lut[:-3, 2] = makeMappingArray(
            self.N, self._segmentdata['blue'], self._gamma)
        if 'alpha' in self._segmentdata:
            self._lut[:-3, 3] = makeMappingArray(
                self.N, self._segmentdata['alpha'], 1)
        self._isinit = True
        self._set_extremes()

    def set_gamma(self, gamma):
        """
        Set a new gamma value and regenerate color map.
        """
        self._gamma = gamma
        self._init()

    @staticmethod
    def from_list(name, colors, N=256, gamma=1.0):
        """
        Make a linear segmented colormap with *name* from a sequence
        of *colors* which evenly transitions from colors[0] at val=0
        to colors[-1] at val=1.  *N* is the number of rgb quantization
        levels.
        Alternatively, a list of (value, color) tuples can be given
        to divide the range unevenly.
        """

        if not cbook.iterable(colors):
            raise ValueError('colors must be iterable')

        if (isinstance(colors[0], Sized) and len(colors[0]) == 2
                and not isinstance(colors[0], six.string_types)):
            # List of value, color pairs
            vals, colors = zip(*colors)
        else:
            vals = np.linspace(0, 1, len(colors))

        cdict = dict(red=[], green=[], blue=[], alpha=[])
        for val, color in zip(vals, colors):
            r, g, b, a = colorConverter.to_rgba(color)
            cdict['red'].append((val, r, r))
            cdict['green'].append((val, g, g))
            cdict['blue'].append((val, b, b))
            cdict['alpha'].append((val, a, a))

        return LinearSegmentedColormap(name, cdict, N, gamma)

    def _resample(self, lutsize):
        """
        Return a new color map with *lutsize* entries.
        """
        return LinearSegmentedColormap(self.name, self._segmentdata, lutsize)

    def reversed(self, name=None):
        """
        Make a reversed instance of the Colormap.

        Parameters
        ----------
        name : str, optional
            The name for the reversed colormap. If it's None the
            name will be the name of the parent colormap + "_r".

        Returns
        -------
        LinearSegmentedColormap
            The reversed colormap.
        """
        if name is None:
            name = self.name + "_r"

        # Function factory needed to deal with 'late binding' issue.
        def factory(dat):
            def func_r(x):
                return dat(1.0 - x)
            return func_r

        data_r = dict()
        for key, data in six.iteritems(self._segmentdata):
            if callable(data):
                data_r[key] = factory(data)
            else:
                new_data = [(1.0 - x, y1, y0) for x, y0, y1 in reversed(data)]
                data_r[key] = new_data

        return LinearSegmentedColormap(name, data_r, self.N, self._gamma)


class ListedColormap(Colormap):
    """Colormap object generated from a list of colors.

    This may be most useful when indexing directly into a colormap,
    but it can also be used to generate special colormaps for ordinary
    mapping.
    """
    def __init__(self, colors, name='from_list', N=None):
        """
        Make a colormap from a list of colors.

        *colors*
            a list of matplotlib color specifications,
            or an equivalent Nx3 or Nx4 floating point array
            (*N* rgb or rgba values)
        *name*
            a string to identify the colormap
        *N*
            the number of entries in the map.  The default is *None*,
            in which case there is one colormap entry for each
            element in the list of colors.  If::

                N < len(colors)

            the list will be truncated at *N*.  If::

                N > len(colors)

            the list will be extended by repetition.
        """
        self.monochrome = False  # True only if all colors in map are
                                 # identical; needed for contouring.
        if N is None:
            self.colors = colors
            N = len(colors)
        else:
            if isinstance(colors, six.string_types):
                self.colors = [colors] * N
                self.monochrome = True
            elif cbook.iterable(colors):
                if len(colors) == 1:
                    self.monochrome = True
                self.colors = list(
                    itertools.islice(itertools.cycle(colors), N))
            else:
                try:
                    gray = float(colors)
                except TypeError:
                    pass
                else:
                    self.colors = [gray] * N
                self.monochrome = True
        Colormap.__init__(self, name, N)

    def _init(self):
        rgba = colorConverter.to_rgba_array(self.colors)
        self._lut = np.zeros((self.N + 3, 4), float)
        self._lut[:-3] = rgba
        self._isinit = True
        self._set_extremes()

    def _resample(self, lutsize):
        """
        Return a new color map with *lutsize* entries.
        """
        colors = self(np.linspace(0, 1, lutsize))
        return ListedColormap(colors, name=self.name)

    def reversed(self, name=None):
        """
        Make a reversed instance of the Colormap.

        Parameters
        ----------
        name : str, optional
            The name for the reversed colormap. If it's None the
            name will be the name of the parent colormap + "_r".

        Returns
        -------
        ListedColormap
            A reversed instance of the colormap.
        """
        if name is None:
            name = self.name + "_r"

        colors_r = list(reversed(self.colors))
        return ListedColormap(colors_r, name=name, N=self.N)


class Normalize(object):
    """
    A class which, when called, can normalize data into
    the ``[0.0, 1.0]`` interval.

    """
    def __init__(self, vmin=None, vmax=None, clip=False):
        """
        If *vmin* or *vmax* is not given, they are initialized from the
        minimum and maximum value respectively of the first input
        processed.  That is, *__call__(A)* calls *autoscale_None(A)*.
        If *clip* is *True* and the given value falls outside the range,
        the returned value will be 0 or 1, whichever is closer.
        Returns 0 if::

            vmin==vmax

        Works with scalars or arrays, including masked arrays.  If
        *clip* is *True*, masked values are set to 1; otherwise they
        remain masked.  Clipping silently defeats the purpose of setting
        the over, under, and masked colors in the colormap, so it is
        likely to lead to surprises; therefore the default is
        *clip* = *False*.
        """
        self.vmin = _sanitize_extrema(vmin)
        self.vmax = _sanitize_extrema(vmax)
        self.clip = clip

    @staticmethod
    def process_value(value):
        """
        Homogenize the input *value* for easy and efficient normalization.

        *value* can be a scalar or sequence.

        Returns *result*, *is_scalar*, where *result* is a
        masked array matching *value*.  Float dtypes are preserved;
        integer types with two bytes or smaller are converted to
        np.float32, and larger types are converted to np.float64.
        Preserving float32 when possible, and using in-place operations,
        can greatly improve speed for large arrays.

        Experimental; we may want to add an option to force the
        use of float32.
        """
        is_scalar = not cbook.iterable(value)
        if is_scalar:
            value = [value]
        dtype = np.min_scalar_type(value)
        if np.issubdtype(dtype, np.integer) or dtype.type is np.bool_:
            # bool_/int8/int16 -> float32; int32/int64 -> float64
            dtype = np.promote_types(dtype, np.float32)
        # ensure data passed in as an ndarray subclass are interpreted as
        # an ndarray. See issue #6622.
        mask = np.ma.getmask(value)
        data = np.asarray(np.ma.getdata(value))
        result = np.ma.array(data, mask=mask, dtype=dtype, copy=True)
        return result, is_scalar

    def __call__(self, value, clip=None):
        """
        Normalize *value* data in the ``[vmin, vmax]`` interval into
        the ``[0.0, 1.0]`` interval and return it.  *clip* defaults
        to *self.clip* (which defaults to *False*).  If not already
        initialized, *vmin* and *vmax* are initialized using
        *autoscale_None(value)*.
        """
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        # Convert at least to float, without losing precision.
        (vmin,), _ = self.process_value(self.vmin)
        (vmax,), _ = self.process_value(self.vmax)
        if vmin == vmax:
            result.fill(0)   # Or should it be all masked?  Or 0.5?
        elif vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            resdat -= vmin
            resdat /= (vmax - vmin)
            result = np.ma.array(resdat, mask=result.mask, copy=False)
        # Agg cannot handle float128.  We actually only need 32-bit of
        # precision, but on Windows, `np.dtype(np.longdouble) == np.float64`,
        # so casting to float32 would lose precision on float64s as well.
        if result.dtype == np.longdouble:
            result = result.astype(np.float64)
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        (vmin,), _ = self.process_value(self.vmin)
        (vmax,), _ = self.process_value(self.vmax)

        if cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin + val * (vmax - vmin)
        else:
            return vmin + value * (vmax - vmin)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        A = np.asanyarray(A)
        self.vmin = A.min()
        self.vmax = A.max()

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        A = np.asanyarray(A)
        if self.vmin is None and A.size:
            self.vmin = A.min()
        if self.vmax is None and A.size:
            self.vmax = A.max()

    def scaled(self):
        'return true if vmin and vmax set'
        return (self.vmin is not None and self.vmax is not None)


class LogNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a log scale
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        result = np.ma.masked_less_equal(result, 0, copy=False)

        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin <= 0:
            raise ValueError("values must all be positive")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            # in-place equivalent of above can be much faster
            resdat = result.data
            mask = result.mask
            if mask is np.ma.nomask:
                mask = (resdat <= 0)
            else:
                mask |= resdat <= 0
            np.copyto(resdat, 1, where=mask)
            np.log(resdat, resdat)
            resdat -= np.log(vmin)
            resdat /= (np.log(vmax) - np.log(vmin))
            result = np.ma.array(resdat, mask=mask, copy=False)
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax

        if cbook.iterable(value):
            val = np.ma.asarray(value)
            return vmin * np.ma.power((vmax / vmin), val)
        else:
            return vmin * pow((vmax / vmin), value)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        A = np.ma.masked_less_equal(A, 0, copy=False)
        self.vmin = np.ma.min(A)
        self.vmax = np.ma.max(A)

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        if self.vmin is not None and self.vmax is not None:
            return
        A = np.ma.masked_less_equal(A, 0, copy=False)
        if self.vmin is None and A.size:
            self.vmin = A.min()
        if self.vmax is None and A.size:
            self.vmax = A.max()


class SymLogNorm(Normalize):
    """
    The symmetrical logarithmic scale is logarithmic in both the
    positive and negative directions from the origin.

    Since the values close to zero tend toward infinity, there is a
    need to have a range around zero that is linear.  The parameter
    *linthresh* allows the user to specify the size of this range
    (-*linthresh*, *linthresh*).
    """
    def __init__(self,  linthresh, linscale=1.0,
                 vmin=None, vmax=None, clip=False):
        """
        *linthresh*:
        The range within which the plot is linear (to
        avoid having the plot go to infinity around zero).

        *linscale*:
        This allows the linear range (-*linthresh* to *linthresh*)
        to be stretched relative to the logarithmic range.  Its
        value is the number of decades to use for each half of the
        linear range.  For example, when *linscale* == 1.0 (the
        default), the space used for the positive and negative
        halves of the linear range will be equal to one decade in
        the logarithmic range. Defaults to 1.
        """
        Normalize.__init__(self, vmin, vmax, clip)
        self.linthresh = float(linthresh)
        self._linscale_adj = (linscale / (1.0 - np.e ** -1))
        if vmin is not None and vmax is not None:
            self._transform_vmin_vmax()

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax

        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            # in-place equivalent of above can be much faster
            resdat = self._transform(result.data)
            resdat -= self._lower
            resdat /= (self._upper - self._lower)

        if is_scalar:
            result = result[0]
        return result

    def _transform(self, a):
        """
        Inplace transformation.
        """
        masked = np.abs(a) > self.linthresh
        sign = np.sign(a[masked])
        log = (self._linscale_adj + np.log(np.abs(a[masked]) / self.linthresh))
        log *= sign * self.linthresh
        a[masked] = log
        a[~masked] *= self._linscale_adj
        return a

    def _inv_transform(self, a):
        """
        Inverse inplace Transformation.
        """
        masked = np.abs(a) > (self.linthresh * self._linscale_adj)
        sign = np.sign(a[masked])
        exp = np.exp(sign * a[masked] / self.linthresh - self._linscale_adj)
        exp *= sign * self.linthresh
        a[masked] = exp
        a[~masked] /= self._linscale_adj
        return a

    def _transform_vmin_vmax(self):
        """
        Calculates vmin and vmax in the transformed system.
        """
        vmin, vmax = self.vmin, self.vmax
        arr = np.array([vmax, vmin]).astype(float)
        self._upper, self._lower = self._transform(arr)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        val = np.ma.asarray(value)
        val = val * (self._upper - self._lower) + self._lower
        return self._inv_transform(val)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = np.ma.min(A)
        self.vmax = np.ma.max(A)
        self._transform_vmin_vmax()

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        if self.vmin is not None and self.vmax is not None:
            pass
        A = np.asanyarray(A)
        if self.vmin is None and A.size:
            self.vmin = A.min()
        if self.vmax is None and A.size:
            self.vmax = A.max()
        self._transform_vmin_vmax()


class PowerNorm(Normalize):
    """
    Normalize a given value to the ``[0, 1]`` interval with a power-law
    scaling. This will clip any negative data points to 0.
    """
    def __init__(self, gamma, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.gamma = gamma

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        gamma = self.gamma
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            res_mask = result.data < 0
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                     mask=mask)
            resdat = result.data
            resdat -= vmin
            np.power(resdat, gamma, resdat)
            resdat /= (vmax - vmin) ** gamma

            result = np.ma.array(resdat, mask=result.mask, copy=False)
            result[res_mask] = 0
        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        gamma = self.gamma
        vmin, vmax = self.vmin, self.vmax

        if cbook.iterable(value):
            val = np.ma.asarray(value)
            return np.ma.power(val, 1. / gamma) * (vmax - vmin) + vmin
        else:
            return pow(value, 1. / gamma) * (vmax - vmin) + vmin

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = np.ma.min(A)
        if self.vmin < 0:
            self.vmin = 0
            warnings.warn("Power-law scaling on negative values is "
                          "ill-defined, clamping to 0.")
        self.vmax = np.ma.max(A)

    def autoscale_None(self, A):
        """autoscale only None-valued vmin or vmax."""
        A = np.asanyarray(A)
        if self.vmin is None and A.size:
            self.vmin = A.min()
            if self.vmin < 0:
                self.vmin = 0
                warnings.warn("Power-law scaling on negative values is "
                              "ill-defined, clamping to 0.")
        if self.vmax is None and A.size:
            self.vmax = A.max()


class BoundaryNorm(Normalize):
    """
    Generate a colormap index based on discrete intervals.

    Unlike :class:`Normalize` or :class:`LogNorm`,
    :class:`BoundaryNorm` maps values to integers instead of to the
    interval 0-1.

    Mapping to the 0-1 interval could have been done via
    piece-wise linear interpolation, but using integers seems
    simpler, and reduces the number of conversions back and forth
    between integer and floating point.
    """
    def __init__(self, boundaries, ncolors, clip=False):
        """
        Parameters
        ----------
        boundaries : array-like
            Monotonically increasing sequence of boundaries
        ncolors : int
            Number of colors in the colormap to be used
        clip : bool, optional
            If clip is ``True``, out of range values are mapped to 0 if they
            are below ``boundaries[0]`` or mapped to ncolors - 1 if they are
            above ``boundaries[-1]``.

            If clip is ``False``, out of range values are mapped to -1 if
            they are below ``boundaries[0]`` or mapped to ncolors if they are
            above ``boundaries[-1]``. These are then converted to valid indices
            by :meth:`Colormap.__call__`.

        Notes
        -----
        *boundaries* defines the edges of bins, and data falling within a bin
        is mapped to the color with the same index.

        If the number of bins doesn't equal *ncolors*, the color is chosen
        by linear interpolation of the bin number onto color numbers.
        """
        self.clip = clip
        self.vmin = boundaries[0]
        self.vmax = boundaries[-1]
        self.boundaries = np.asarray(boundaries)
        self.N = len(self.boundaries)
        self.Ncmap = ncolors
        if self.N - 1 == self.Ncmap:
            self._interp = False
        else:
            self._interp = True

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        xx, is_scalar = self.process_value(value)
        mask = np.ma.getmaskarray(xx)
        xx = np.atleast_1d(xx.filled(self.vmax + 1))
        if clip:
            np.clip(xx, self.vmin, self.vmax, out=xx)
            max_col = self.Ncmap - 1
        else:
            max_col = self.Ncmap
        iret = np.zeros(xx.shape, dtype=np.int16)
        for i, b in enumerate(self.boundaries):
            iret[xx >= b] = i
        if self._interp:
            scalefac = (self.Ncmap - 1) / (self.N - 2)
            iret = (iret * scalefac).astype(np.int16)
        iret[xx < self.vmin] = -1
        iret[xx >= self.vmax] = max_col
        ret = np.ma.array(iret, mask=mask)
        if is_scalar:
            ret = int(ret[0])  # assume python scalar
        return ret

    def inverse(self, value):
        """
        Raises
        ------
        ValueError
            BoundaryNorm is not invertible, so calling this method will always
            raise an error
        """
        return ValueError("BoundaryNorm is not invertible")


class NoNorm(Normalize):
    """
    Dummy replacement for Normalize, for the case where we
    want to use indices directly in a
    :class:`~matplotlib.cm.ScalarMappable` .
    """
    def __call__(self, value, clip=None):
        return value

    def inverse(self, value):
        return value


def rgb_to_hsv(arr):
    """
    convert float rgb values (in the range [0, 1]), in a numpy array to hsv
    values.

    Parameters
    ----------
    arr : (..., 3) array-like
       All values must be in the range [0, 1]

    Returns
    -------
    hsv : (..., 3) ndarray
       Colors converted to hsv values in range [0, 1]
    """
    # make sure it is an ndarray
    arr = np.asarray(arr)

    # check length of the last dimension, should be _some_ sort of rgb
    if arr.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {} was found.".format(arr.shape))

    in_ndim = arr.ndim
    if arr.ndim == 1:
        arr = np.array(arr, ndmin=2)

    # make sure we don't have an int image
    arr = arr.astype(np.promote_types(arr.dtype, np.float32))

    out = np.zeros_like(arr)
    arr_max = arr.max(-1)
    ipos = arr_max > 0
    delta = arr.ptp(-1)
    s = np.zeros_like(delta)
    s[ipos] = delta[ipos] / arr_max[ipos]
    ipos = delta > 0
    # red is max
    idx = (arr[..., 0] == arr_max) & ipos
    out[idx, 0] = (arr[idx, 1] - arr[idx, 2]) / delta[idx]
    # green is max
    idx = (arr[..., 1] == arr_max) & ipos
    out[idx, 0] = 2. + (arr[idx, 2] - arr[idx, 0]) / delta[idx]
    # blue is max
    idx = (arr[..., 2] == arr_max) & ipos
    out[idx, 0] = 4. + (arr[idx, 0] - arr[idx, 1]) / delta[idx]

    out[..., 0] = (out[..., 0] / 6.0) % 1.0
    out[..., 1] = s
    out[..., 2] = arr_max

    if in_ndim == 1:
        out.shape = (3,)

    return out


def hsv_to_rgb(hsv):
    """
    convert hsv values in a numpy array to rgb values
    all values assumed to be in range [0, 1]

    Parameters
    ----------
    hsv : (..., 3) array-like
       All values assumed to be in range [0, 1]

    Returns
    -------
    rgb : (..., 3) ndarray
       Colors converted to RGB values in range [0, 1]
    """
    hsv = np.asarray(hsv)

    # check length of the last dimension, should be _some_ sort of rgb
    if hsv.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {shp} was found.".format(shp=hsv.shape))

    # if we got passed a 1D array, try to treat as
    # a single color and reshape as needed
    in_ndim = hsv.ndim
    if in_ndim == 1:
        hsv = np.array(hsv, ndmin=2)

    # make sure we don't have an int image
    hsv = hsv.astype(np.promote_types(hsv.dtype, np.float32))

    h = hsv[..., 0]
    s = hsv[..., 1]
    v = hsv[..., 2]

    r = np.empty_like(h)
    g = np.empty_like(h)
    b = np.empty_like(h)

    i = (h * 6.0).astype(int)
    f = (h * 6.0) - i
    p = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))

    idx = i % 6 == 0
    r[idx] = v[idx]
    g[idx] = t[idx]
    b[idx] = p[idx]

    idx = i == 1
    r[idx] = q[idx]
    g[idx] = v[idx]
    b[idx] = p[idx]

    idx = i == 2
    r[idx] = p[idx]
    g[idx] = v[idx]
    b[idx] = t[idx]

    idx = i == 3
    r[idx] = p[idx]
    g[idx] = q[idx]
    b[idx] = v[idx]

    idx = i == 4
    r[idx] = t[idx]
    g[idx] = p[idx]
    b[idx] = v[idx]

    idx = i == 5
    r[idx] = v[idx]
    g[idx] = p[idx]
    b[idx] = q[idx]

    idx = s == 0
    r[idx] = v[idx]
    g[idx] = v[idx]
    b[idx] = v[idx]

    # `np.stack([r, g, b], axis=-1)` (numpy 1.10).
    rgb = np.concatenate([r[..., None], g[..., None], b[..., None]], -1)

    if in_ndim == 1:
        rgb.shape = (3,)

    return rgb


def _vector_magnitude(arr):
    # things that don't work here:
    #  * np.linalg.norm
    #    - doesn't broadcast in numpy 1.7
    #    - drops the mask from ma.array
    #  * using keepdims - broken on ma.array until 1.11.2
    #  * using sum - discards mask on ma.array unless entire vector is masked

    sum_sq = 0
    for i in range(arr.shape[-1]):
        sum_sq += np.square(arr[..., i, np.newaxis])
    return np.sqrt(sum_sq)


def _vector_dot(a, b):
    # things that don't work here:
    #   * a.dot(b) - fails on masked arrays until 1.10
    #   * np.ma.dot(a, b) - doesn't mask enough things
    #   * np.ma.dot(a, b, strict=True) - returns a maskedarray with no mask
    dot = 0
    for i in range(a.shape[-1]):
        dot += a[..., i] * b[..., i]
    return dot


class LightSource(object):
    """
    Create a light source coming from the specified azimuth and elevation.
    Angles are in degrees, with the azimuth measured
    clockwise from north and elevation up from the zero plane of the surface.

    The :meth:`shade` is used to produce "shaded" rgb values for a data array.
    :meth:`shade_rgb` can be used to combine an rgb image with
    The :meth:`shade_rgb`
    The :meth:`hillshade` produces an illumination map of a surface.
    """
    def __init__(self, azdeg=315, altdeg=45, hsv_min_val=0, hsv_max_val=1,
                 hsv_min_sat=1, hsv_max_sat=0):
        """
        Specify the azimuth (measured clockwise from south) and altitude
        (measured up from the plane of the surface) of the light source
        in degrees.

        Parameters
        ----------
        azdeg : number, optional
            The azimuth (0-360, degrees clockwise from North) of the light
            source. Defaults to 315 degrees (from the northwest).
        altdeg : number, optional
            The altitude (0-90, degrees up from horizontal) of the light
            source.  Defaults to 45 degrees from horizontal.

        Notes
        -----
        For backwards compatibility, the parameters *hsv_min_val*,
        *hsv_max_val*, *hsv_min_sat*, and *hsv_max_sat* may be supplied at
        initialization as well.  However, these parameters will only be used if
        "blend_mode='hsv'" is passed into :meth:`shade` or :meth:`shade_rgb`.
        See the documentation for :meth:`blend_hsv` for more details.
        """
        self.azdeg = azdeg
        self.altdeg = altdeg
        self.hsv_min_val = hsv_min_val
        self.hsv_max_val = hsv_max_val
        self.hsv_min_sat = hsv_min_sat
        self.hsv_max_sat = hsv_max_sat

    @property
    def direction(self):
        """ The unit vector direction towards the light source """

        # Azimuth is in degrees clockwise from North. Convert to radians
        # counterclockwise from East (mathematical notation).
        az = np.radians(90 - self.azdeg)
        alt = np.radians(self.altdeg)

        return np.array([
            np.cos(az) * np.cos(alt),
            np.sin(az) * np.cos(alt),
            np.sin(alt)
        ])

    def hillshade(self, elevation, vert_exag=1, dx=1, dy=1, fraction=1.):
        """
        Calculates the illumination intensity for a surface using the defined
        azimuth and elevation for the light source.

        This computes the normal vectors for the surface, and then passes them
        on to `shade_normals`

        Parameters
        ----------
        elevation : array-like
            A 2d array (or equivalent) of the height values used to generate an
            illumination map
        vert_exag : number, optional
            The amount to exaggerate the elevation values by when calculating
            illumination. This can be used either to correct for differences in
            units between the x-y coordinate system and the elevation
            coordinate system (e.g. decimal degrees vs meters) or to exaggerate
            or de-emphasize topographic effects.
        dx : number, optional
            The x-spacing (columns) of the input *elevation* grid.
        dy : number, optional
            The y-spacing (rows) of the input *elevation* grid.
        fraction : number, optional
            Increases or decreases the contrast of the hillshade.  Values
            greater than one will cause intermediate values to move closer to
            full illumination or shadow (and clipping any values that move
            beyond 0 or 1). Note that this is not visually or mathematically
            the same as vertical exaggeration.
        Returns
        -------
        intensity : ndarray
            A 2d array of illumination values between 0-1, where 0 is
            completely in shadow and 1 is completely illuminated.
        """

        # Because most image and raster GIS data has the first row in the array
        # as the "top" of the image, dy is implicitly negative.  This is
        # consistent to what `imshow` assumes, as well.
        dy = -dy

        # compute the normal vectors from the partial derivatives
        e_dy, e_dx = np.gradient(vert_exag * elevation, dy, dx)

        # .view is to keep subclasses
        normal = np.empty(elevation.shape + (3,)).view(type(elevation))
        normal[..., 0] = -e_dx
        normal[..., 1] = -e_dy
        normal[..., 2] = 1
        normal /= _vector_magnitude(normal)

        return self.shade_normals(normal, fraction)

    def shade_normals(self, normals, fraction=1.):
        """
        Calculates the illumination intensity for the normal vectors of a
        surface using the defined azimuth and elevation for the light source.

        Imagine an artificial sun placed at infinity in some azimuth and
        elevation position illuminating our surface. The parts of the surface
        that slope toward the sun should brighten while those sides facing away
        should become darker.

        Parameters
        ----------
        fraction : number, optional
            Increases or decreases the contrast of the hillshade.  Values
            greater than one will cause intermediate values to move closer to
            full illumination or shadow (and clipping any values that move
            beyond 0 or 1). Note that this is not visually or mathematically
            the same as vertical exaggeration.

        Returns
        -------
        intensity : ndarray
            A 2d array of illumination values between 0-1, where 0 is
            completely in shadow and 1 is completely illuminated.
        """

        intensity = _vector_dot(normals, self.direction)

        # Apply contrast stretch
        imin, imax = intensity.min(), intensity.max()
        intensity *= fraction

        # Rescale to 0-1, keeping range before contrast stretch
        # If constant slope, keep relative scaling (i.e. flat should be 0.5,
        # fully occluded 0, etc.)
        if (imax - imin) > 1e-6:
            # Strictly speaking, this is incorrect. Negative values should be
            # clipped to 0 because they're fully occluded. However, rescaling
            # in this manner is consistent with the previous implementation and
            # visually appears better than a "hard" clip.
            intensity -= imin
            intensity /= (imax - imin)
        intensity = np.clip(intensity, 0, 1, intensity)

        return intensity

    def shade(self, data, cmap, norm=None, blend_mode='overlay', vmin=None,
              vmax=None, vert_exag=1, dx=1, dy=1, fraction=1, **kwargs):
        """
        Combine colormapped data values with an illumination intensity map
        (a.k.a.  "hillshade") of the values.

        Parameters
        ----------
        data : array-like
            A 2d array (or equivalent) of the height values used to generate a
            shaded map.
        cmap : `~matplotlib.colors.Colormap` instance
            The colormap used to color the *data* array. Note that this must be
            a `~matplotlib.colors.Colormap` instance.  For example, rather than
            passing in `cmap='gist_earth'`, use
            `cmap=plt.get_cmap('gist_earth')` instead.
        norm : `~matplotlib.colors.Normalize` instance, optional
            The normalization used to scale values before colormapping. If
            None, the input will be linearly scaled between its min and max.
        blend_mode : {'hsv', 'overlay', 'soft'} or callable, optional
            The type of blending used to combine the colormapped data
            values with the illumination intensity.  Default is
            "overlay".  Note that for most topographic surfaces,
            "overlay" or "soft" appear more visually realistic. If a
            user-defined function is supplied, it is expected to
            combine an MxNx3 RGB array of floats (ranging 0 to 1) with
            an MxNx1 hillshade array (also 0 to 1).  (Call signature
            `func(rgb, illum, **kwargs)`) Additional kwargs supplied
            to this function will be passed on to the *blend_mode*
            function.
        vmin : scalar or None, optional
            The minimum value used in colormapping *data*. If *None* the
            minimum value in *data* is used. If *norm* is specified, then this
            argument will be ignored.
        vmax : scalar or None, optional
            The maximum value used in colormapping *data*. If *None* the
            maximum value in *data* is used. If *norm* is specified, then this
            argument will be ignored.
        vert_exag : number, optional
            The amount to exaggerate the elevation values by when calculating
            illumination. This can be used either to correct for differences in
            units between the x-y coordinate system and the elevation
            coordinate system (e.g. decimal degrees vs meters) or to exaggerate
            or de-emphasize topography.
        dx : number, optional
            The x-spacing (columns) of the input *elevation* grid.
        dy : number, optional
            The y-spacing (rows) of the input *elevation* grid.
        fraction : number, optional
            Increases or decreases the contrast of the hillshade.  Values
            greater than one will cause intermediate values to move closer to
            full illumination or shadow (and clipping any values that move
            beyond 0 or 1). Note that this is not visually or mathematically
            the same as vertical exaggeration.
        Additional kwargs are passed on to the *blend_mode* function.

        Returns
        -------
        rgba : ndarray
            An MxNx4 array of floats ranging between 0-1.
        """
        if vmin is None:
            vmin = data.min()
        if vmax is None:
            vmax = data.max()
        if norm is None:
            norm = Normalize(vmin=vmin, vmax=vmax)

        rgb0 = cmap(norm(data))
        rgb1 = self.shade_rgb(rgb0, elevation=data, blend_mode=blend_mode,
                              vert_exag=vert_exag, dx=dx, dy=dy,
                              fraction=fraction, **kwargs)
        # Don't overwrite the alpha channel, if present.
        rgb0[..., :3] = rgb1[..., :3]
        return rgb0

    def shade_rgb(self, rgb, elevation, fraction=1., blend_mode='hsv',
                  vert_exag=1, dx=1, dy=1, **kwargs):
        """
        Take the input RGB array (ny*nx*3) adjust their color values
        to given the impression of a shaded relief map with a
        specified light source using the elevation (ny*nx).
        A new RGB array ((ny*nx*3)) is returned.

        Parameters
        ----------
        rgb : array-like
            An MxNx3 RGB array, assumed to be in the range of 0 to 1.
        elevation : array-like
            A 2d array (or equivalent) of the height values used to generate a
            shaded map.
        fraction : number
            Increases or decreases the contrast of the hillshade.  Values
            greater than one will cause intermediate values to move closer to
            full illumination or shadow (and clipping any values that move
            beyond 0 or 1). Note that this is not visually or mathematically
            the same as vertical exaggeration.
        blend_mode : {'hsv', 'overlay', 'soft'} or callable, optional
            The type of blending used to combine the colormapped data values
            with the illumination intensity.  For backwards compatibility, this
            defaults to "hsv". Note that for most topographic surfaces,
            "overlay" or "soft" appear more visually realistic. If a
            user-defined function is supplied, it is expected to combine an
            MxNx3 RGB array of floats (ranging 0 to 1) with an MxNx1 hillshade
            array (also 0 to 1).  (Call signature `func(rgb, illum, **kwargs)`)
            Additional kwargs supplied to this function will be passed on to
            the *blend_mode* function.
        vert_exag : number, optional
            The amount to exaggerate the elevation values by when calculating
            illumination. This can be used either to correct for differences in
            units between the x-y coordinate system and the elevation
            coordinate system (e.g. decimal degrees vs meters) or to exaggerate
            or de-emphasize topography.
        dx : number, optional
            The x-spacing (columns) of the input *elevation* grid.
        dy : number, optional
            The y-spacing (rows) of the input *elevation* grid.
        Additional kwargs are passed on to the *blend_mode* function.

        Returns
        -------
        shaded_rgb : ndarray
            An MxNx3 array of floats ranging between 0-1.
        """
        # Calculate the "hillshade" intensity.
        intensity = self.hillshade(elevation, vert_exag, dx, dy, fraction)
        intensity = intensity[..., np.newaxis]

        # Blend the hillshade and rgb data using the specified mode
        lookup = {
                'hsv': self.blend_hsv,
                'soft': self.blend_soft_light,
                'overlay': self.blend_overlay,
                }
        if blend_mode in lookup:
            blend = lookup[blend_mode](rgb, intensity, **kwargs)
        else:
            try:
                blend = blend_mode(rgb, intensity, **kwargs)
            except TypeError:
                raise ValueError('"blend_mode" must be callable or one of {}'
                                 .format(lookup.keys))

        # Only apply result where hillshade intensity isn't masked
        if hasattr(intensity, 'mask'):
            mask = intensity.mask[..., 0]
            for i in range(3):
                blend[..., i][mask] = rgb[..., i][mask]

        return blend

    def blend_hsv(self, rgb, intensity, hsv_max_sat=None, hsv_max_val=None,
                  hsv_min_val=None, hsv_min_sat=None):
        """
        Take the input data array, convert to HSV values in the given colormap,
        then adjust those color values to give the impression of a shaded
        relief map with a specified light source.  RGBA values are returned,
        which can then be used to plot the shaded image with imshow.

        The color of the resulting image will be darkened by moving the (s,v)
        values (in hsv colorspace) toward (hsv_min_sat, hsv_min_val) in the
        shaded regions, or lightened by sliding (s,v) toward (hsv_max_sat
        hsv_max_val) in regions that are illuminated.  The default extremes are
        chose so that completely shaded points are nearly black (s = 1, v = 0)
        and completely illuminated points are nearly white (s = 0, v = 1).

        Parameters
        ----------
        rgb : ndarray
            An MxNx3 RGB array of floats ranging from 0 to 1 (color image).
        intensity : ndarray
            An MxNx1 array of floats ranging from 0 to 1 (grayscale image).
        hsv_max_sat : number, optional
            The maximum saturation value that the *intensity* map can shift the
            output image to. Defaults to 1.
        hsv_min_sat : number, optional
            The minimum saturation value that the *intensity* map can shift the
            output image to. Defaults to 0.
        hsv_max_val : number, optional
            The maximum value ("v" in "hsv") that the *intensity* map can shift
            the output image to. Defaults to 1.
        hsv_min_val: number, optional
            The minimum value ("v" in "hsv") that the *intensity* map can shift
            the output image to. Defaults to 0.

        Returns
        -------
        rgb : ndarray
            An MxNx3 RGB array representing the combined images.
        """
        # Backward compatibility...
        if hsv_max_sat is None:
            hsv_max_sat = self.hsv_max_sat
        if hsv_max_val is None:
            hsv_max_val = self.hsv_max_val
        if hsv_min_sat is None:
            hsv_min_sat = self.hsv_min_sat
        if hsv_min_val is None:
            hsv_min_val = self.hsv_min_val

        # Expects a 2D intensity array scaled between -1 to 1...
        intensity = intensity[..., 0]
        intensity = 2 * intensity - 1

        # convert to rgb, then rgb to hsv
        hsv = rgb_to_hsv(rgb[:, :, 0:3])

        # modify hsv values to simulate illumination.
        hsv[:, :, 1] = np.where(np.logical_and(np.abs(hsv[:, :, 1]) > 1.e-10,
                                               intensity > 0),
                                ((1. - intensity) * hsv[:, :, 1] +
                                 intensity * hsv_max_sat),
                                hsv[:, :, 1])

        hsv[:, :, 2] = np.where(intensity > 0,
                                ((1. - intensity) * hsv[:, :, 2] +
                                 intensity * hsv_max_val),
                                hsv[:, :, 2])

        hsv[:, :, 1] = np.where(np.logical_and(np.abs(hsv[:, :, 1]) > 1.e-10,
                                               intensity < 0),
                                ((1. + intensity) * hsv[:, :, 1] -
                                 intensity * hsv_min_sat),
                                hsv[:, :, 1])
        hsv[:, :, 2] = np.where(intensity < 0,
                                ((1. + intensity) * hsv[:, :, 2] -
                                 intensity * hsv_min_val),
                                hsv[:, :, 2])
        hsv[:, :, 1:] = np.where(hsv[:, :, 1:] < 0., 0, hsv[:, :, 1:])
        hsv[:, :, 1:] = np.where(hsv[:, :, 1:] > 1., 1, hsv[:, :, 1:])
        # convert modified hsv back to rgb.
        return hsv_to_rgb(hsv)

    def blend_soft_light(self, rgb, intensity):
        """
        Combines an rgb image with an intensity map using "soft light"
        blending.  Uses the "pegtop" formula.

        Parameters
        ----------
        rgb : ndarray
            An MxNx3 RGB array of floats ranging from 0 to 1 (color image).
        intensity : ndarray
            An MxNx1 array of floats ranging from 0 to 1 (grayscale image).

        Returns
        -------
        rgb : ndarray
            An MxNx3 RGB array representing the combined images.
        """
        return 2 * intensity * rgb + (1 - 2 * intensity) * rgb**2

    def blend_overlay(self, rgb, intensity):
        """
        Combines an rgb image with an intensity map using "overlay" blending.

        Parameters
        ----------
        rgb : ndarray
            An MxNx3 RGB array of floats ranging from 0 to 1 (color image).
        intensity : ndarray
            An MxNx1 array of floats ranging from 0 to 1 (grayscale image).

        Returns
        -------
        rgb : ndarray
            An MxNx3 RGB array representing the combined images.
        """
        low = 2 * intensity * rgb
        high = 1 - 2 * (1 - intensity) * (1 - rgb)
        return np.where(rgb <= 0.5, low, high)


def from_levels_and_colors(levels, colors, extend='neither'):
    """
    A helper routine to generate a cmap and a norm instance which
    behave similar to contourf's levels and colors arguments.

    Parameters
    ----------
    levels : sequence of numbers
        The quantization levels used to construct the :class:`BoundaryNorm`.
        Values ``v`` are quantizized to level ``i`` if
        ``lev[i] <= v < lev[i+1]``.
    colors : sequence of colors
        The fill color to use for each level. If `extend` is "neither" there
        must be ``n_level - 1`` colors. For an `extend` of "min" or "max" add
        one extra color, and for an `extend` of "both" add two colors.
    extend : {'neither', 'min', 'max', 'both'}, optional
        The behaviour when a value falls out of range of the given levels.
        See :func:`~matplotlib.pyplot.contourf` for details.

    Returns
    -------
    (cmap, norm) : tuple containing a :class:`Colormap` and a \
                   :class:`Normalize` instance
    """
    colors_i0 = 0
    colors_i1 = None

    if extend == 'both':
        colors_i0 = 1
        colors_i1 = -1
        extra_colors = 2
    elif extend == 'min':
        colors_i0 = 1
        extra_colors = 1
    elif extend == 'max':
        colors_i1 = -1
        extra_colors = 1
    elif extend == 'neither':
        extra_colors = 0
    else:
        raise ValueError('Unexpected value for extend: {0!r}'.format(extend))

    n_data_colors = len(levels) - 1
    n_expected_colors = n_data_colors + extra_colors
    if len(colors) != n_expected_colors:
        raise ValueError('With extend == {0!r} and n_levels == {1!r} expected'
                         ' n_colors == {2!r}. Got {3!r}.'
                         ''.format(extend, len(levels), n_expected_colors,
                                   len(colors)))

    cmap = ListedColormap(colors[colors_i0:colors_i1], N=n_data_colors)

    if extend in ['min', 'both']:
        cmap.set_under(colors[0])
    else:
        cmap.set_under('none')

    if extend in ['max', 'both']:
        cmap.set_over(colors[-1])
    else:
        cmap.set_over('none')

    cmap.colorbar_extend = extend

    norm = BoundaryNorm(levels, ncolors=n_data_colors)
    return cmap, norm
