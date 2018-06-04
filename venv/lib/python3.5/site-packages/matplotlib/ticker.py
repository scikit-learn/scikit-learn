"""
Tick locating and formatting
============================

This module contains classes to support completely configurable tick
locating and formatting. Although the locators know nothing about major
or minor ticks, they are used by the Axis class to support major and
minor tick locating and formatting. Generic tick locators and
formatters are provided, as well as domain specific custom ones.

Default Formatter
-----------------

The default formatter identifies when the x-data being plotted is a
small range on top of a large off set. To reduce the chances that the
ticklabels overlap the ticks are labeled as deltas from a fixed offset.
For example::

   ax.plot(np.arange(2000, 2010), range(10))

will have tick of 0-9 with an offset of +2e3. If this is not desired
turn off the use of the offset on the default formatter::

   ax.get_xaxis().get_major_formatter().set_useOffset(False)

set the rcParam ``axes.formatter.useoffset=False`` to turn it off
globally, or set a different formatter.

Tick locating
-------------

The Locator class is the base class for all tick locators. The locators
handle autoscaling of the view limits based on the data limits, and the
choosing of tick locations. A useful semi-automatic tick locator is
`MultipleLocator`. It is initialized with a base, e.g., 10, and it picks
axis limits and ticks that are multiples of that base.

The Locator subclasses defined here are

:class:`AutoLocator`
    `MaxNLocator` with simple defaults.  This is the default tick locator for
    most plotting.

:class:`MaxNLocator`
    Finds up to a max number of intervals with ticks at nice locations.

:class:`LinearLocator`
    Space ticks evenly from min to max.

:class:`LogLocator`
    Space ticks logarithmically from min to max.

:class:`MultipleLocator`
    Ticks and range are a multiple of base; either integer or float.

:class:`FixedLocator`
    Tick locations are fixed.

:class:`IndexLocator`
    Locator for index plots (e.g., where ``x = range(len(y))``).

:class:`NullLocator`
    No ticks.

:class:`SymmetricalLogLocator`
    Locator for use with with the symlog norm; works like `LogLocator` for the
    part outside of the threshold and adds 0 if inside the limits.

:class:`LogitLocator`
    Locator for logit scaling.

:class:`OldAutoLocator`
    Choose a `MultipleLocator` and dynamically reassign it for intelligent
    ticking during navigation.

:class:`AutoMinorLocator`
    Locator for minor ticks when the axis is linear and the
    major ticks are uniformly spaced.  Subdivides the major
    tick interval into a specified number of minor intervals,
    defaulting to 4 or 5 depending on the major interval.


There are a number of locators specialized for date locations - see
the `dates` module.

You can define your own locator by deriving from Locator. You must
override the ``__call__`` method, which returns a sequence of locations,
and you will probably want to override the autoscale method to set the
view limits from the data limits.

If you want to override the default locator, use one of the above or a custom
locator and pass it to the x or y axis instance. The relevant methods are::

  ax.xaxis.set_major_locator(xmajor_locator)
  ax.xaxis.set_minor_locator(xminor_locator)
  ax.yaxis.set_major_locator(ymajor_locator)
  ax.yaxis.set_minor_locator(yminor_locator)

The default minor locator is `NullLocator`, i.e., no minor ticks on by default.

Tick formatting
---------------

Tick formatting is controlled by classes derived from Formatter. The formatter
operates on a single tick value and returns a string to the axis.

:class:`NullFormatter`
    No labels on the ticks.

:class:`IndexFormatter`
    Set the strings from a list of labels.

:class:`FixedFormatter`
    Set the strings manually for the labels.

:class:`FuncFormatter`
    User defined function sets the labels.

:class:`StrMethodFormatter`
    Use string `format` method.

:class:`FormatStrFormatter`
    Use an old-style sprintf format string.

:class:`ScalarFormatter`
    Default formatter for scalars: autopick the format string.

:class:`LogFormatter`
    Formatter for log axes.

:class:`LogFormatterExponent`
    Format values for log axis using ``exponent = log_base(value)``.

:class:`LogFormatterMathtext`
    Format values for log axis using ``exponent = log_base(value)``
    using Math text.

:class:`LogFormatterSciNotation`
    Format values for log axis using scientific notation.

:class:`LogitFormatter`
    Probability formatter.

:class:`EngFormatter`
    Format labels in engineering notation

:class:`PercentFormatter`
    Format labels as a percentage

You can derive your own formatter from the Formatter base class by
simply overriding the ``__call__`` method. The formatter class has
access to the axis view and data limits.

To control the major and minor tick label formats, use one of the
following methods::

  ax.xaxis.set_major_formatter(xmajor_formatter)
  ax.xaxis.set_minor_formatter(xminor_formatter)
  ax.yaxis.set_major_formatter(ymajor_formatter)
  ax.yaxis.set_minor_formatter(yminor_formatter)

See :ref:`sphx_glr_gallery_ticks_and_spines_major_minor_demo.py` for an
example of setting major and minor ticks. See the :mod:`matplotlib.dates`
module for more information and examples of using date locators and formatters.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import itertools
import locale
import math
import numpy as np
from matplotlib import rcParams
from matplotlib import cbook
from matplotlib import transforms as mtransforms
from matplotlib.cbook import mplDeprecation

import warnings


__all__ = ('TickHelper', 'Formatter', 'FixedFormatter',
           'NullFormatter', 'FuncFormatter', 'FormatStrFormatter',
           'StrMethodFormatter', 'ScalarFormatter', 'LogFormatter',
           'LogFormatterExponent', 'LogFormatterMathtext',
           'IndexFormatter', 'LogFormatterSciNotation',
           'LogitFormatter', 'EngFormatter', 'PercentFormatter',
           'Locator', 'IndexLocator', 'FixedLocator', 'NullLocator',
           'LinearLocator', 'LogLocator', 'AutoLocator',
           'MultipleLocator', 'MaxNLocator', 'AutoMinorLocator',
           'SymmetricalLogLocator', 'LogitLocator')


if six.PY3:
    long = int


# Work around numpy/numpy#6127.
def _divmod(x, y):
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(y, np.generic):
        y = y.item()
    return six.moves.builtins.divmod(x, y)


def _mathdefault(s):
    return '\\mathdefault{%s}' % s


class _DummyAxis(object):
    def __init__(self, minpos=0):
        self.dataLim = mtransforms.Bbox.unit()
        self.viewLim = mtransforms.Bbox.unit()
        self._minpos = minpos

    def get_view_interval(self):
        return self.viewLim.intervalx

    def set_view_interval(self, vmin, vmax):
        self.viewLim.intervalx = vmin, vmax

    def get_minpos(self):
        return self._minpos

    def get_data_interval(self):
        return self.dataLim.intervalx

    def set_data_interval(self, vmin, vmax):
        self.dataLim.intervalx = vmin, vmax

    def get_tick_space(self):
        # Just use the long-standing default of nbins==9
        return 9


class TickHelper(object):
    axis = None

    def set_axis(self, axis):
        self.axis = axis

    def create_dummy_axis(self, **kwargs):
        if self.axis is None:
            self.axis = _DummyAxis(**kwargs)

    def set_view_interval(self, vmin, vmax):
        self.axis.set_view_interval(vmin, vmax)

    def set_data_interval(self, vmin, vmax):
        self.axis.set_data_interval(vmin, vmax)

    def set_bounds(self, vmin, vmax):
        self.set_view_interval(vmin, vmax)
        self.set_data_interval(vmin, vmax)


class Formatter(TickHelper):
    """
    Create a string based on a tick value and location.
    """
    # some classes want to see all the locs to help format
    # individual ones
    locs = []

    def __call__(self, x, pos=None):
        """
        Return the format for tick value `x` at position pos.
        ``pos=None`` indicates an unspecified location.
        """
        raise NotImplementedError('Derived must override')

    def format_data(self, value):
        """
        Returns the full string representation of the value with the
        position unspecified.
        """
        return self.__call__(value)

    def format_data_short(self, value):
        """
        Return a short string version of the tick value.

        Defaults to the position-independent long value.
        """
        return self.format_data(value)

    def get_offset(self):
        return ''

    def set_locs(self, locs):
        self.locs = locs

    def fix_minus(self, s):
        """
        Some classes may want to replace a hyphen for minus with the
        proper unicode symbol (U+2212) for typographical correctness.
        The default is to not replace it.

        Note, if you use this method, e.g., in :meth:`format_data` or
        call, you probably don't want to use it for
        :meth:`format_data_short` since the toolbar uses this for
        interactive coord reporting and I doubt we can expect GUIs
        across platforms will handle the unicode correctly.  So for
        now the classes that override :meth:`fix_minus` should have an
        explicit :meth:`format_data_short` method
        """
        return s


class IndexFormatter(Formatter):
    """
    Format the position x to the nearest i-th label where i=int(x+0.5)
    """
    def __init__(self, labels):
        self.labels = labels
        self.n = len(labels)

    def __call__(self, x, pos=None):
        """
        Return the format for tick value `x` at position pos.

        The position is ignored and the value is rounded to the nearest
        integer, which is used to look up the label.
        """
        i = int(x + 0.5)
        if i < 0 or i >= self.n:
            return ''
        else:
            return self.labels[i]


class NullFormatter(Formatter):
    """
    Always return the empty string.
    """
    def __call__(self, x, pos=None):
        """
        Returns an empty string for all inputs.
        """
        return ''


class FixedFormatter(Formatter):
    """
    Return fixed strings for tick labels based only on position, not
    value.
    """
    def __init__(self, seq):
        """
        Set the sequence of strings that will be used for labels.
        """
        self.seq = seq
        self.offset_string = ''

    def __call__(self, x, pos=None):
        """
        Returns the label that matches the position regardless of the
        value.

        For positions ``pos < len(seq)``, return `seq[i]` regardless of
        `x`. Otherwise return empty string. `seq` is the sequence of
        strings that this object was initialized with.
        """
        if pos is None or pos >= len(self.seq):
            return ''
        else:
            return self.seq[pos]

    def get_offset(self):
        return self.offset_string

    def set_offset_string(self, ofs):
        self.offset_string = ofs


class FuncFormatter(Formatter):
    """
    Use a user-defined function for formatting.

    The function should take in two inputs (a tick value ``x`` and a
    position ``pos``), and return a string containing the corresponding
    tick label.
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, x, pos=None):
        """
        Return the value of the user defined function.

        `x` and `pos` are passed through as-is.
        """
        return self.func(x, pos)


class FormatStrFormatter(Formatter):
    """
    Use an old-style ('%' operator) format string to format the tick.

    The format string should have a single variable format (%) in it.
    It will be applied to the value (not the position) of the tick.
    """
    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, x, pos=None):
        """
        Return the formatted label string.

        Only the value `x` is formatted. The position is ignored.
        """
        return self.fmt % x


class StrMethodFormatter(Formatter):
    """
    Use a new-style format string (as used by `str.format()`)
    to format the tick.

    The field used for the value must be labeled `x` and the field used
    for the position must be labeled `pos`.
    """
    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, x, pos=None):
        """
        Return the formatted label string.

        `x` and `pos` are passed to `str.format` as keyword arguments
        with those exact names.
        """
        return self.fmt.format(x=x, pos=pos)


class OldScalarFormatter(Formatter):
    """
    Tick location is a plain old number.
    """

    def __call__(self, x, pos=None):
        """
        Return the format for tick val `x` based on the width of the
        axis.

        The position `pos` is ignored.
        """
        xmin, xmax = self.axis.get_view_interval()
        d = abs(xmax - xmin)

        return self.pprint_val(x, d)

    def pprint_val(self, x, d):
        """
        Formats the value `x` based on the size of the axis range `d`.
        """
        #if the number is not too big and it's an int, format it as an
        #int
        if abs(x) < 1e4 and x == int(x):
            return '%d' % x

        if d < 1e-2:
            fmt = '%1.3e'
        elif d < 1e-1:
            fmt = '%1.3f'
        elif d > 1e5:
            fmt = '%1.1e'
        elif d > 10:
            fmt = '%1.1f'
        elif d > 1:
            fmt = '%1.2f'
        else:
            fmt = '%1.3f'
        s = fmt % x
        tup = s.split('e')
        if len(tup) == 2:
            mantissa = tup[0].rstrip('0').rstrip('.')
            sign = tup[1][0].replace('+', '')
            exponent = tup[1][1:].lstrip('0')
            s = '%se%s%s' % (mantissa, sign, exponent)
        else:
            s = s.rstrip('0').rstrip('.')
        return s


class ScalarFormatter(Formatter):
    """
    Format tick values as a number.

    Tick value is interpreted as a plain old number. If
    ``useOffset==True`` and the data range is much smaller than the data
    average, then an offset will be determined such that the tick labels
    are meaningful. Scientific notation is used for ``data < 10^-n`` or
    ``data >= 10^m``, where ``n`` and ``m`` are the power limits set
    using ``set_powerlimits((n,m))``. The defaults for these are
    controlled by the ``axes.formatter.limits`` rc parameter.
    """
    def __init__(self, useOffset=None, useMathText=None, useLocale=None):
        # useOffset allows plotting small data ranges with large offsets: for
        # example: [1+1e-9,1+2e-9,1+3e-9] useMathText will render the offset
        # and scientific notation in mathtext

        if useOffset is None:
            useOffset = rcParams['axes.formatter.useoffset']
        self._offset_threshold = rcParams['axes.formatter.offset_threshold']
        self.set_useOffset(useOffset)
        self._usetex = rcParams['text.usetex']
        if useMathText is None:
            useMathText = rcParams['axes.formatter.use_mathtext']
        self.set_useMathText(useMathText)
        self.orderOfMagnitude = 0
        self.format = ''
        self._scientific = True
        self._powerlimits = rcParams['axes.formatter.limits']
        if useLocale is None:
            useLocale = rcParams['axes.formatter.use_locale']
        self._useLocale = useLocale

    def get_useOffset(self):
        return self._useOffset

    def set_useOffset(self, val):
        if val in [True, False]:
            self.offset = 0
            self._useOffset = val
        else:
            self._useOffset = False
            self.offset = val

    useOffset = property(fget=get_useOffset, fset=set_useOffset)

    def get_useLocale(self):
        return self._useLocale

    def set_useLocale(self, val):
        if val is None:
            self._useLocale = rcParams['axes.formatter.use_locale']
        else:
            self._useLocale = val

    useLocale = property(fget=get_useLocale, fset=set_useLocale)

    def get_useMathText(self):
        return self._useMathText

    def set_useMathText(self, val):
        if val is None:
            self._useMathText = rcParams['axes.formatter.use_mathtext']
        else:
            self._useMathText = val

    useMathText = property(fget=get_useMathText, fset=set_useMathText)

    def fix_minus(self, s):
        """
        Replace hyphens with a unicode minus.
        """
        if rcParams['text.usetex'] or not rcParams['axes.unicode_minus']:
            return s
        else:
            return s.replace('-', '\N{MINUS SIGN}')

    def __call__(self, x, pos=None):
        """
        Return the format for tick value `x` at position `pos`.
        """
        if len(self.locs) == 0:
            return ''
        else:
            s = self.pprint_val(x)
            return self.fix_minus(s)

    def set_scientific(self, b):
        """
        Turn scientific notation on or off.

        .. seealso:: Method :meth:`set_powerlimits`
        """
        self._scientific = bool(b)

    def set_powerlimits(self, lims):
        """
        Sets size thresholds for scientific notation.

        ``lims`` is a two-element sequence containing the powers of 10
        that determine the switchover threshold. Numbers below
        ``10**lims[0]`` and above ``10**lims[1]`` will be displayed in
        scientific notation.

        For example, ``formatter.set_powerlimits((-3, 4))`` sets the
        pre-2007 default in which scientific notation is used for
        numbers less than 1e-3 or greater than 1e4.

        .. seealso:: Method :meth:`set_scientific`
        """
        if len(lims) != 2:
            raise ValueError("'lims' must be a sequence of length 2")
        self._powerlimits = lims

    def format_data_short(self, value):
        """
        Return a short formatted string representation of a number.
        """
        if self._useLocale:
            return locale.format_string('%-12g', (value,))
        else:
            return '%-12g' % value

    def format_data(self, value):
        """
        Return a formatted string representation of a number.
        """
        if self._useLocale:
            s = locale.format_string('%1.10e', (value,))
        else:
            s = '%1.10e' % value
        s = self._formatSciNotation(s)
        return self.fix_minus(s)

    def get_offset(self):
        """
        Return scientific notation, plus offset.
        """
        if len(self.locs) == 0:
            return ''
        s = ''
        if self.orderOfMagnitude or self.offset:
            offsetStr = ''
            sciNotStr = ''
            if self.offset:
                offsetStr = self.format_data(self.offset)
                if self.offset > 0:
                    offsetStr = '+' + offsetStr
            if self.orderOfMagnitude:
                if self._usetex or self._useMathText:
                    sciNotStr = self.format_data(10 ** self.orderOfMagnitude)
                else:
                    sciNotStr = '1e%d' % self.orderOfMagnitude
            if self._useMathText:
                if sciNotStr != '':
                    sciNotStr = r'\times%s' % _mathdefault(sciNotStr)
                s = ''.join(('$', sciNotStr, _mathdefault(offsetStr), '$'))
            elif self._usetex:
                if sciNotStr != '':
                    sciNotStr = r'\times%s' % sciNotStr
                s = ''.join(('$', sciNotStr, offsetStr, '$'))
            else:
                s = ''.join((sciNotStr, offsetStr))

        return self.fix_minus(s)

    def set_locs(self, locs):
        """
        Set the locations of the ticks.
        """
        self.locs = locs
        if len(self.locs) > 0:
            vmin, vmax = self.axis.get_view_interval()
            d = abs(vmax - vmin)
            if self._useOffset:
                self._compute_offset()
            self._set_orderOfMagnitude(d)
            self._set_format(vmin, vmax)

    def _compute_offset(self):
        locs = self.locs
        if locs is None or not len(locs):
            self.offset = 0
            return
        # Restrict to visible ticks.
        vmin, vmax = sorted(self.axis.get_view_interval())
        locs = np.asarray(locs)
        locs = locs[(vmin <= locs) & (locs <= vmax)]
        if not len(locs):
            self.offset = 0
            return
        lmin, lmax = locs.min(), locs.max()
        # Only use offset if there are at least two ticks and every tick has
        # the same sign.
        if lmin == lmax or lmin <= 0 <= lmax:
            self.offset = 0
            return
        # min, max comparing absolute values (we want division to round towards
        # zero so we work on absolute values).
        abs_min, abs_max = sorted([abs(float(lmin)), abs(float(lmax))])
        sign = math.copysign(1, lmin)
        # What is the smallest power of ten such that abs_min and abs_max are
        # equal up to that precision?
        # Note: Internally using oom instead of 10 ** oom avoids some numerical
        # accuracy issues.
        oom_max = np.ceil(math.log10(abs_max))
        oom = 1 + next(oom for oom in itertools.count(oom_max, -1)
                       if abs_min // 10 ** oom != abs_max // 10 ** oom)
        if (abs_max - abs_min) / 10 ** oom <= 1e-2:
            # Handle the case of straddling a multiple of a large power of ten
            # (relative to the span).
            # What is the smallest power of ten such that abs_min and abs_max
            # are no more than 1 apart at that precision?
            oom = 1 + next(oom for oom in itertools.count(oom_max, -1)
                           if abs_max // 10 ** oom - abs_min // 10 ** oom > 1)
        # Only use offset if it saves at least _offset_threshold digits.
        n = self._offset_threshold - 1
        self.offset = (sign * (abs_max // 10 ** oom) * 10 ** oom
                       if abs_max // 10 ** oom >= 10**n
                       else 0)

    def _set_orderOfMagnitude(self, range):
        # if scientific notation is to be used, find the appropriate exponent
        # if using an numerical offset, find the exponent after applying the
        # offset
        if not self._scientific:
            self.orderOfMagnitude = 0
            return
        locs = np.abs(self.locs)
        if self.offset:
            oom = math.floor(math.log10(range))
        else:
            if locs[0] > locs[-1]:
                val = locs[0]
            else:
                val = locs[-1]
            if val == 0:
                oom = 0
            else:
                oom = math.floor(math.log10(val))
        if oom <= self._powerlimits[0]:
            self.orderOfMagnitude = oom
        elif oom >= self._powerlimits[1]:
            self.orderOfMagnitude = oom
        else:
            self.orderOfMagnitude = 0

    def _set_format(self, vmin, vmax):
        # set the format string to format all the ticklabels
        if len(self.locs) < 2:
            # Temporarily augment the locations with the axis end points.
            _locs = list(self.locs) + [vmin, vmax]
        else:
            _locs = self.locs
        locs = (np.asarray(_locs) - self.offset) / 10. ** self.orderOfMagnitude
        loc_range = np.ptp(locs)
        # Curvilinear coordinates can yield two identical points.
        if loc_range == 0:
            loc_range = np.max(np.abs(locs))
        # Both points might be zero.
        if loc_range == 0:
            loc_range = 1
        if len(self.locs) < 2:
            # We needed the end points only for the loc_range calculation.
            locs = locs[:-2]
        loc_range_oom = int(math.floor(math.log10(loc_range)))
        # first estimate:
        sigfigs = max(0, 3 - loc_range_oom)
        # refined estimate:
        thresh = 1e-3 * 10 ** loc_range_oom
        while sigfigs >= 0:
            if np.abs(locs - np.round(locs, decimals=sigfigs)).max() < thresh:
                sigfigs -= 1
            else:
                break
        sigfigs += 1
        self.format = '%1.' + str(sigfigs) + 'f'
        if self._usetex:
            self.format = '$%s$' % self.format
        elif self._useMathText:
            self.format = '$%s$' % _mathdefault(self.format)

    def pprint_val(self, x):
        xp = (x - self.offset) / (10. ** self.orderOfMagnitude)
        if np.abs(xp) < 1e-8:
            xp = 0
        if self._useLocale:
            return locale.format_string(self.format, (xp,))
        else:
            return self.format % xp

    def _formatSciNotation(self, s):
        # transform 1e+004 into 1e4, for example
        if self._useLocale:
            decimal_point = locale.localeconv()['decimal_point']
            positive_sign = locale.localeconv()['positive_sign']
        else:
            decimal_point = '.'
            positive_sign = '+'
        tup = s.split('e')
        try:
            significand = tup[0].rstrip('0').rstrip(decimal_point)
            sign = tup[1][0].replace(positive_sign, '')
            exponent = tup[1][1:].lstrip('0')
            if self._useMathText or self._usetex:
                if significand == '1' and exponent != '':
                    # reformat 1x10^y as 10^y
                    significand = ''
                if exponent:
                    exponent = '10^{%s%s}' % (sign, exponent)
                if significand and exponent:
                    return r'%s{\times}%s' % (significand, exponent)
                else:
                    return r'%s%s' % (significand, exponent)
            else:
                s = ('%se%s%s' % (significand, sign, exponent)).rstrip('e')
                return s
        except IndexError:
            return s


class LogFormatter(Formatter):
    """
    Base class for formatting ticks on a log or symlog scale.

    It may be instantiated directly, or subclassed.

    Parameters
    ----------
    base : float, optional, default: 10.
        Base of the logarithm used in all calculations.

    labelOnlyBase : bool, optional, default: False
        If True, label ticks only at integer powers of base.
        This is normally True for major ticks and False for
        minor ticks.

    minor_thresholds : (subset, all), optional, default: (1, 0.4)
        If labelOnlyBase is False, these two numbers control
        the labeling of ticks that are not at integer powers of
        base; normally these are the minor ticks. The controlling
        parameter is the log of the axis data range.  In the typical
        case where base is 10 it is the number of decades spanned
        by the axis, so we can call it 'numdec'. If ``numdec <= all``,
        all minor ticks will be labeled.  If ``all < numdec <= subset``,
        then only a subset of minor ticks will be labeled, so as to
        avoid crowding. If ``numdec > subset`` then no minor ticks will
        be labeled.

    linthresh : None or float, optional, default: None
        If a symmetric log scale is in use, its ``linthresh``
        parameter must be supplied here.

    Notes
    -----
    The `set_locs` method must be called to enable the subsetting
    logic controlled by the ``minor_thresholds`` parameter.

    In some cases such as the colorbar, there is no distinction between
    major and minor ticks; the tick locations might be set manually,
    or by a locator that puts ticks at integer powers of base and
    at intermediate locations.  For this situation, disable the
    minor_thresholds logic by using ``minor_thresholds=(np.inf, np.inf)``,
    so that all ticks will be labeled.

    To disable labeling of minor ticks when 'labelOnlyBase' is False,
    use ``minor_thresholds=(0, 0)``.  This is the default for the
    "classic" style.

    Examples
    --------
    To label a subset of minor ticks when the view limits span up
    to 2 decades, and all of the ticks when zoomed in to 0.5 decades
    or less, use ``minor_thresholds=(2, 0.5)``.

    To label all minor ticks when the view limits span up to 1.5
    decades, use ``minor_thresholds=(1.5, 1.5)``.

    """
    def __init__(self, base=10.0, labelOnlyBase=False,
                 minor_thresholds=None,
                 linthresh=None):

        self._base = float(base)
        self.labelOnlyBase = labelOnlyBase
        if minor_thresholds is None:
            if rcParams['_internal.classic_mode']:
                minor_thresholds = (0, 0)
            else:
                minor_thresholds = (1, 0.4)
        self.minor_thresholds = minor_thresholds
        self._sublabels = None
        self._linthresh = linthresh

    def base(self, base):
        """
        change the `base` for labeling.

        .. warning::
           Should always match the base used for :class:`LogLocator`

        """
        self._base = base

    def label_minor(self, labelOnlyBase):
        """
        Switch minor tick labeling on or off.

        Parameters
        ----------
        labelOnlyBase : bool
            If True, label ticks only at integer powers of base.

        """
        self.labelOnlyBase = labelOnlyBase

    def set_locs(self, locs=None):
        """
        Use axis view limits to control which ticks are labeled.

        The ``locs`` parameter is ignored in the present algorithm.

        """
        if np.isinf(self.minor_thresholds[0]):
            self._sublabels = None
            return

        # Handle symlog case:
        linthresh = self._linthresh
        if linthresh is None:
            try:
                linthresh = self.axis.get_transform().linthresh
            except AttributeError:
                pass

        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin, vmax = vmax, vmin

        if linthresh is None and vmin <= 0:
            # It's probably a colorbar with
            # a format kwarg setting a LogFormatter in the manner
            # that worked with 1.5.x, but that doesn't work now.
            self._sublabels = set((1,))  # label powers of base
            return

        b = self._base
        if linthresh is not None:  # symlog
            # Only compute the number of decades in the logarithmic part of the
            # axis
            numdec = 0
            if vmin < -linthresh:
                rhs = min(vmax, -linthresh)
                numdec += math.log(vmin / rhs) / math.log(b)
            if vmax > linthresh:
                lhs = max(vmin, linthresh)
                numdec += math.log(vmax / lhs) / math.log(b)
        else:
            vmin = math.log(vmin) / math.log(b)
            vmax = math.log(vmax) / math.log(b)
            numdec = abs(vmax - vmin)

        if numdec > self.minor_thresholds[0]:
            # Label only bases
            self._sublabels = {1}
        elif numdec > self.minor_thresholds[1]:
            # Add labels between bases at log-spaced coefficients;
            # include base powers in case the locations include
            # "major" and "minor" points, as in colorbar.
            c = np.logspace(0, 1, int(b)//2 + 1, base=b)
            self._sublabels = set(np.round(c))
            # For base 10, this yields (1, 2, 3, 4, 6, 10).
        else:
            # Label all integer multiples of base**n.
            self._sublabels = set(np.arange(1, b + 1))

    def _num_to_string(self, x, vmin, vmax):
        if x > 10000:
            s = '%1.0e' % x
        elif x < 1:
            s = '%1.0e' % x
        else:
            s = self.pprint_val(x, vmax - vmin)
        return s

    def __call__(self, x, pos=None):
        """
        Return the format for tick val `x`.
        """
        if x == 0.0:  # Symlog
            return '0'

        x = abs(x)
        b = self._base
        # only label the decades
        fx = math.log(x) / math.log(b)
        is_x_decade = is_close_to_int(fx)
        exponent = np.round(fx) if is_x_decade else np.floor(fx)
        coeff = np.round(x / b ** exponent)

        if self.labelOnlyBase and not is_x_decade:
            return ''
        if self._sublabels is not None and coeff not in self._sublabels:
            return ''

        vmin, vmax = self.axis.get_view_interval()
        vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander=0.05)
        s = self._num_to_string(x, vmin, vmax)
        return self.fix_minus(s)

    def format_data(self, value):
        b = self.labelOnlyBase
        self.labelOnlyBase = False
        value = cbook.strip_math(self.__call__(value))
        self.labelOnlyBase = b
        return value

    def format_data_short(self, value):
        """
        Return a short formatted string representation of a number.
        """
        return '%-12g' % value

    def pprint_val(self, x, d):
        #if the number is not too big and it's an int, format it as an
        #int
        if abs(x) < 1e4 and x == int(x):
            return '%d' % x

        if d < 1e-2:
            fmt = '%1.3e'
        elif d < 1e-1:
            fmt = '%1.3f'
        elif d > 1e5:
            fmt = '%1.1e'
        elif d > 10:
            fmt = '%1.1f'
        elif d > 1:
            fmt = '%1.2f'
        else:
            fmt = '%1.3f'
        s = fmt % x

        tup = s.split('e')
        if len(tup) == 2:
            mantissa = tup[0].rstrip('0').rstrip('.')
            exponent = int(tup[1])
            if exponent:
                s = '%se%d' % (mantissa, exponent)
            else:
                s = mantissa
        else:
            s = s.rstrip('0').rstrip('.')
        return s


class LogFormatterExponent(LogFormatter):
    """
    Format values for log axis using ``exponent = log_base(value)``.
    """
    def _num_to_string(self, x, vmin, vmax):
        fx = math.log(x) / math.log(self._base)
        if abs(fx) > 10000:
            s = '%1.0g' % fx
        elif abs(fx) < 1:
            s = '%1.0g' % fx
        else:
            fd = math.log(vmax - vmin) / math.log(self._base)
            s = self.pprint_val(fx, fd)
        return s


class LogFormatterMathtext(LogFormatter):
    """
    Format values for log axis using ``exponent = log_base(value)``.
    """

    def _non_decade_format(self, sign_string, base, fx, usetex):
        'Return string for non-decade locations'
        if usetex:
            return (r'$%s%s^{%.2f}$') % (sign_string, base, fx)
        else:
            return ('$%s$' % _mathdefault('%s%s^{%.2f}' %
                                          (sign_string, base, fx)))

    def __call__(self, x, pos=None):
        """
        Return the format for tick value `x`.

        The position `pos` is ignored.
        """
        usetex = rcParams['text.usetex']
        min_exp = rcParams['axes.formatter.min_exponent']

        if x == 0:  # Symlog
            if usetex:
                return '$0$'
            else:
                return '$%s$' % _mathdefault('0')

        sign_string = '-' if x < 0 else ''
        x = abs(x)
        b = self._base

        # only label the decades
        fx = math.log(x) / math.log(b)
        is_x_decade = is_close_to_int(fx)
        exponent = np.round(fx) if is_x_decade else np.floor(fx)
        coeff = np.round(x / b ** exponent)
        if is_x_decade:
            fx = nearest_long(fx)

        if self.labelOnlyBase and not is_x_decade:
            return ''
        if self._sublabels is not None and coeff not in self._sublabels:
            return ''

        # use string formatting of the base if it is not an integer
        if b % 1 == 0.0:
            base = '%d' % b
        else:
            base = '%s' % b

        if np.abs(fx) < min_exp:
            if usetex:
                return r'${0}{1:g}$'.format(sign_string, x)
            else:
                return '${0}$'.format(_mathdefault(
                    '{0}{1:g}'.format(sign_string, x)))
        elif not is_x_decade:
            return self._non_decade_format(sign_string, base, fx, usetex)
        else:
            if usetex:
                return (r'$%s%s^{%d}$') % (sign_string,
                                           base,
                                           nearest_long(fx))
            else:
                return ('$%s$' % _mathdefault(
                    '%s%s^{%d}' %
                    (sign_string, base, nearest_long(fx))))


class LogFormatterSciNotation(LogFormatterMathtext):
    """
    Format values following scientific notation in a logarithmic axis
    """

    def _non_decade_format(self, sign_string, base, fx, usetex):
        'Return string for non-decade locations'
        b = float(base)
        exponent = math.floor(fx)
        coeff = b ** fx / b ** exponent
        if is_close_to_int(coeff):
            coeff = nearest_long(coeff)
        if usetex:
            return (r'$%s%g\times%s^{%d}$') % \
                                        (sign_string, coeff, base, exponent)
        else:
            return ('$%s$' % _mathdefault(r'%s%g\times%s^{%d}' %
                                        (sign_string, coeff, base, exponent)))


class LogitFormatter(Formatter):
    """
    Probability formatter (using Math text).
    """
    def __call__(self, x, pos=None):
        s = ''
        if 0.01 <= x <= 0.99:
            s = '{:.2f}'.format(x)
        elif x < 0.01:
            if is_decade(x):
                s = '$10^{{{:.0f}}}$'.format(np.log10(x))
            else:
                s = '${:.5f}$'.format(x)
        else:  # x > 0.99
            if is_decade(1-x):
                s = '$1-10^{{{:.0f}}}$'.format(np.log10(1-x))
            else:
                s = '$1-{:.5f}$'.format(1-x)
        return s

    def format_data_short(self, value):
        'return a short formatted string representation of a number'
        return '%-12g' % value


class EngFormatter(Formatter):
    """
    Formats axis values using engineering prefixes to represent powers
    of 1000, plus a specified unit, e.g., 10 MHz instead of 1e7.
    """

    # The SI engineering prefixes
    ENG_PREFIXES = {
        -24: "y",
        -21: "z",
        -18: "a",
        -15: "f",
        -12: "p",
         -9: "n",
         -6: "\N{GREEK SMALL LETTER MU}",
         -3: "m",
          0: "",
          3: "k",
          6: "M",
          9: "G",
         12: "T",
         15: "P",
         18: "E",
         21: "Z",
         24: "Y"
    }

    def __init__(self, unit="", places=None, sep=" "):
        """
        Parameters
        ----------
        unit : str (default: "")
            Unit symbol to use, suitable for use with single-letter
            representations of powers of 1000. For example, 'Hz' or 'm'.

        places : int (default: None)
            Precision with which to display the number, specified in
            digits after the decimal point (there will be between one
            and three digits before the decimal point). If it is None,
            the formatting falls back to the floating point format '%g',
            which displays up to 6 *significant* digits, i.e. the equivalent
            value for *places* varies between 0 and 5 (inclusive).

        sep : str (default: " ")
            Separator used between the value and the prefix/unit. For
            example, one get '3.14 mV' if ``sep`` is " " (default) and
            '3.14mV' if ``sep`` is "". Besides the default behavior, some
            other useful options may be:

            * ``sep=""`` to append directly the prefix/unit to the value;
            * ``sep="\\N{THIN SPACE}"`` (``U+2009``);
            * ``sep="\\N{NARROW NO-BREAK SPACE}"`` (``U+202F``);
            * ``sep="\\N{NO-BREAK SPACE}"`` (``U+00A0``).
        """
        self.unit = unit
        self.places = places
        self.sep = sep

    def __call__(self, x, pos=None):
        s = "%s%s" % (self.format_eng(x), self.unit)
        # Remove the trailing separator when there is neither prefix nor unit
        if len(self.sep) > 0 and s.endswith(self.sep):
            s = s[:-len(self.sep)]
        return self.fix_minus(s)

    def format_eng(self, num):
        """
        Formats a number in engineering notation, appending a letter
        representing the power of 1000 of the original number.
        Some examples:

        >>> format_eng(0)       # for self.places = 0
        '0'

        >>> format_eng(1000000) # for self.places = 1
        '1.0 M'

        >>> format_eng("-1e-6") # for self.places = 2
        u'-1.00 \N{GREEK SMALL LETTER MU}'

        `num` may be a numeric value or a string that can be converted
        to a numeric value with ``float(num)``.
        """
        if isinstance(num, six.string_types):
            warnings.warn(
                "Passing a string as *num* argument is deprecated since"
                "Matplotlib 2.1, and is expected to be removed in 2.3.",
                mplDeprecation)

        dnum = float(num)
        sign = 1
        fmt = "g" if self.places is None else ".{:d}f".format(self.places)

        if dnum < 0:
            sign = -1
            dnum = -dnum

        if dnum != 0:
            pow10 = int(math.floor(math.log10(dnum) / 3) * 3)
        else:
            pow10 = 0
            # Force dnum to zero, to avoid inconsistencies like
            # format_eng(-0) = "0" and format_eng(0.0) = "0"
            # but format_eng(-0.0) = "-0.0"
            dnum = 0.0

        pow10 = np.clip(pow10, min(self.ENG_PREFIXES), max(self.ENG_PREFIXES))

        mant = sign * dnum / (10.0 ** pow10)
        # Taking care of the cases like 999.9..., which
        # may be rounded to 1000 instead of 1 k.  Beware
        # of the corner case of values that are beyond
        # the range of SI prefixes (i.e. > 'Y').
        _fmant = float("{mant:{fmt}}".format(mant=mant, fmt=fmt))
        if _fmant >= 1000 and pow10 != max(self.ENG_PREFIXES):
            mant /= 1000
            pow10 += 3

        prefix = self.ENG_PREFIXES[int(pow10)]

        formatted = "{mant:{fmt}}{sep}{prefix}".format(
            mant=mant, sep=self.sep, prefix=prefix, fmt=fmt)

        return formatted


class PercentFormatter(Formatter):
    """
    Format numbers as a percentage.

    How the number is converted into a percentage is determined by the
    `xmax` parameter. `xmax` is the data value that corresponds to 100%.
    Percentages are computed as ``x / xmax * 100``. So if the data is
    already scaled to be percentages, `xmax` will be 100. Another common
    situation is where `xmax` is 1.0.

    `symbol` is a string which will be appended to the label. It may be
    `None` or empty to indicate that no symbol should be used. LaTeX
    special characters are escaped in `symbol` whenever latex mode is
    enabled, unless `is_latex` is `True`.

    `decimals` is the number of decimal places to place after the point.
    If it is set to `None` (the default), the number will be computed
    automatically.
    """
    def __init__(self, xmax=100, decimals=None, symbol='%', is_latex=False):
        self.xmax = xmax + 0.0
        self.decimals = decimals
        self._symbol = symbol
        self._is_latex = is_latex

    def __call__(self, x, pos=None):
        """
        Formats the tick as a percentage with the appropriate scaling.
        """
        ax_min, ax_max = self.axis.get_view_interval()
        display_range = abs(ax_max - ax_min)

        return self.fix_minus(self.format_pct(x, display_range))

    def format_pct(self, x, display_range):
        """
        Formats the number as a percentage number with the correct
        number of decimals and adds the percent symbol, if any.

        If `self.decimals` is `None`, the number of digits after the
        decimal point is set based on the `display_range` of the axis
        as follows:

        +---------------+----------+------------------------+
        | display_range | decimals |          sample        |
        +---------------+----------+------------------------+
        | >50           |     0    | ``x = 34.5`` => 35%    |
        +---------------+----------+------------------------+
        | >5            |     1    | ``x = 34.5`` => 34.5%  |
        +---------------+----------+------------------------+
        | >0.5          |     2    | ``x = 34.5`` => 34.50% |
        +---------------+----------+------------------------+
        |      ...      |    ...   |          ...           |
        +---------------+----------+------------------------+

        This method will not be very good for tiny axis ranges or
        extremely large ones. It assumes that the values on the chart
        are percentages displayed on a reasonable scale.
        """
        x = self.convert_to_pct(x)
        if self.decimals is None:
            # conversion works because display_range is a difference
            scaled_range = self.convert_to_pct(display_range)
            if scaled_range <= 0:
                decimals = 0
            else:
                # Luckily Python's built-in ceil rounds to +inf, not away from
                # zero. This is very important since the equation for decimals
                # starts out as `scaled_range > 0.5 * 10**(2 - decimals)`
                # and ends up with `decimals > 2 - log10(2 * scaled_range)`.
                decimals = math.ceil(2.0 - math.log10(2.0 * scaled_range))
                if decimals > 5:
                    decimals = 5
                elif decimals < 0:
                    decimals = 0
        else:
            decimals = self.decimals
        s = '{x:0.{decimals}f}'.format(x=x, decimals=int(decimals))

        return s + self.symbol

    def convert_to_pct(self, x):
        return 100.0 * (x / self.xmax)

    @property
    def symbol(self):
        """
        The configured percent symbol as a string.

        If LaTeX is enabled via :rc:`text.usetex`, the special characters
        ``{'#', '$', '%', '&', '~', '_', '^', '\\', '{', '}'}`` are
        automatically escaped in the string.
        """
        symbol = self._symbol
        if not symbol:
            symbol = ''
        elif rcParams['text.usetex'] and not self._is_latex:
            # Source: http://www.personal.ceu.hu/tex/specchar.htm
            # Backslash must be first for this to work correctly since
            # it keeps getting added in
            for spec in r'\#$%&~_^{}':
                symbol = symbol.replace(spec, '\\' + spec)
        return symbol

    @symbol.setter
    def symbol(self, symbol):
        self._symbol = symbol


class Locator(TickHelper):
    """
    Determine the tick locations;

    Note, you should not use the same locator between different
    :class:`~matplotlib.axis.Axis` because the locator stores references to
    the Axis data and view limits
    """

    # Some automatic tick locators can generate so many ticks they
    # kill the machine when you try and render them.
    # This parameter is set to cause locators to raise an error if too
    # many ticks are generated.
    MAXTICKS = 1000

    def tick_values(self, vmin, vmax):
        """
        Return the values of the located ticks given **vmin** and **vmax**.

        .. note::
            To get tick locations with the vmin and vmax values defined
            automatically for the associated :attr:`axis` simply call
            the Locator instance::

                >>> print((type(loc)))
                <type 'Locator'>
                >>> print((loc()))
                [1, 2, 3, 4]

        """
        raise NotImplementedError('Derived must override')

    def set_params(self, **kwargs):
        """
        Do nothing, and rase a warning. Any locator class not supporting the
        set_params() function will call this.
        """
        warnings.warn("'set_params()' not defined for locator of type " +
                      str(type(self)))

    def __call__(self):
        """Return the locations of the ticks"""
        # note: some locators return data limits, other return view limits,
        # hence there is no *one* interface to call self.tick_values.
        raise NotImplementedError('Derived must override')

    def raise_if_exceeds(self, locs):
        """raise a RuntimeError if Locator attempts to create more than
           MAXTICKS locs"""
        if len(locs) >= self.MAXTICKS:
            raise RuntimeError("Locator attempting to generate {} ticks from "
                               "{} to {}: exceeds Locator.MAXTICKS".format(
                                   len(locs), locs[0], locs[-1]))
        return locs

    def view_limits(self, vmin, vmax):
        """
        select a scale for the range from vmin to vmax

        Normally this method is overridden by subclasses to
        change locator behaviour.
        """
        return mtransforms.nonsingular(vmin, vmax)

    def autoscale(self):
        """autoscale the view limits"""
        return self.view_limits(*self.axis.get_view_interval())

    def pan(self, numsteps):
        """Pan numticks (can be positive or negative)"""
        ticks = self()
        numticks = len(ticks)

        vmin, vmax = self.axis.get_view_interval()
        vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander=0.05)
        if numticks > 2:
            step = numsteps * abs(ticks[0] - ticks[1])
        else:
            d = abs(vmax - vmin)
            step = numsteps * d / 6.

        vmin += step
        vmax += step
        self.axis.set_view_interval(vmin, vmax, ignore=True)

    def zoom(self, direction):
        "Zoom in/out on axis; if direction is >0 zoom in, else zoom out"

        vmin, vmax = self.axis.get_view_interval()
        vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander=0.05)
        interval = abs(vmax - vmin)
        step = 0.1 * interval * direction
        self.axis.set_view_interval(vmin + step, vmax - step, ignore=True)

    def refresh(self):
        """refresh internal information based on current lim"""
        pass


class IndexLocator(Locator):
    """
    Place a tick on every multiple of some base number of points
    plotted, e.g., on every 5th point.  It is assumed that you are doing
    index plotting; i.e., the axis is 0, len(data).  This is mainly
    useful for x ticks.
    """
    def __init__(self, base, offset):
        'place ticks on the i-th data points where (i-offset)%base==0'
        self._base = base
        self.offset = offset

    def set_params(self, base=None, offset=None):
        """Set parameters within this locator"""
        if base is not None:
            self._base = base
        if offset is not None:
            self.offset = offset

    def __call__(self):
        """Return the locations of the ticks"""
        dmin, dmax = self.axis.get_data_interval()
        return self.tick_values(dmin, dmax)

    def tick_values(self, vmin, vmax):
        return self.raise_if_exceeds(
            np.arange(vmin + self.offset, vmax + 1, self._base))


class FixedLocator(Locator):
    """
    Tick locations are fixed.  If nbins is not None,
    the array of possible positions will be subsampled to
    keep the number of ticks <= nbins +1.
    The subsampling will be done so as to include the smallest
    absolute value; for example, if zero is included in the
    array of possibilities, then it is guaranteed to be one of
    the chosen ticks.
    """

    def __init__(self, locs, nbins=None):
        self.locs = np.asarray(locs)
        self.nbins = nbins
        if self.nbins is not None:
            self.nbins = max(self.nbins, 2)

    def set_params(self, nbins=None):
        """Set parameters within this locator."""
        if nbins is not None:
            self.nbins = nbins

    def __call__(self):
        return self.tick_values(None, None)

    def tick_values(self, vmin, vmax):
        """"
        Return the locations of the ticks.

        .. note::

            Because the values are fixed, vmin and vmax are not used in this
            method.

        """
        if self.nbins is None:
            return self.locs
        step = max(int(np.ceil(len(self.locs) / self.nbins)), 1)
        ticks = self.locs[::step]
        for i in range(1, step):
            ticks1 = self.locs[i::step]
            if np.abs(ticks1).min() < np.abs(ticks).min():
                ticks = ticks1
        return self.raise_if_exceeds(ticks)


class NullLocator(Locator):
    """
    No ticks
    """

    def __call__(self):
        return self.tick_values(None, None)

    def tick_values(self, vmin, vmax):
        """"
        Return the locations of the ticks.

        .. note::

            Because the values are Null, vmin and vmax are not used in this
            method.
        """
        return []


class LinearLocator(Locator):
    """
    Determine the tick locations

    The first time this function is called it will try to set the
    number of ticks to make a nice tick partitioning.  Thereafter the
    number of ticks will be fixed so that interactive navigation will
    be nice

    """
    def __init__(self, numticks=None, presets=None):
        """
        Use presets to set locs based on lom.  A dict mapping vmin, vmax->locs
        """
        self.numticks = numticks
        if presets is None:
            self.presets = {}
        else:
            self.presets = presets

    def set_params(self, numticks=None, presets=None):
        """Set parameters within this locator."""
        if presets is not None:
            self.presets = presets
        if numticks is not None:
            self.numticks = numticks

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander=0.05)
        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if (vmin, vmax) in self.presets:
            return self.presets[(vmin, vmax)]

        if self.numticks is None:
            self._set_numticks()

        if self.numticks == 0:
            return []
        ticklocs = np.linspace(vmin, vmax, self.numticks)

        return self.raise_if_exceeds(ticklocs)

    def _set_numticks(self):
        self.numticks = 11  # todo; be smart here; this is just for dev

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if vmin == vmax:
            vmin -= 1
            vmax += 1

        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            exponent, remainder = _divmod(
                math.log10(vmax - vmin), math.log10(max(self.numticks - 1, 1)))
            exponent -= (remainder < .5)
            scale = max(self.numticks - 1, 1) ** (-exponent)
            vmin = math.floor(scale * vmin) / scale
            vmax = math.ceil(scale * vmax) / scale

        return mtransforms.nonsingular(vmin, vmax)


def closeto(x, y):
    if abs(x - y) < 1e-10:
        return True
    else:
        return False


class Base(object):
    'this solution has some hacks to deal with floating point inaccuracies'
    def __init__(self, base):
        if base <= 0:
            raise ValueError("'base' must be positive")
        self._base = base

    def lt(self, x):
        'return the largest multiple of base < x'
        d, m = _divmod(x, self._base)
        if closeto(m, 0) and not closeto(m / self._base, 1):
            return (d - 1) * self._base
        return d * self._base

    def le(self, x):
        'return the largest multiple of base <= x'
        d, m = _divmod(x, self._base)
        if closeto(m / self._base, 1):  # was closeto(m, self._base)
            #looks like floating point error
            return (d + 1) * self._base
        return d * self._base

    def gt(self, x):
        'return the smallest multiple of base > x'
        d, m = _divmod(x, self._base)
        if closeto(m / self._base, 1):
            #looks like floating point error
            return (d + 2) * self._base
        return (d + 1) * self._base

    def ge(self, x):
        'return the smallest multiple of base >= x'
        d, m = _divmod(x, self._base)
        if closeto(m, 0) and not closeto(m / self._base, 1):
            return d * self._base
        return (d + 1) * self._base

    def get_base(self):
        return self._base


class MultipleLocator(Locator):
    """
    Set a tick on every integer that is multiple of base in the
    view interval
    """

    def __init__(self, base=1.0):
        self._base = Base(base)

    def set_params(self, base):
        """Set parameters within this locator."""
        if base is not None:
            self._base = base

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        vmin = self._base.ge(vmin)
        base = self._base.get_base()
        n = (vmax - vmin + 0.001 * base) // base
        locs = vmin - base + np.arange(n + 3) * base
        return self.raise_if_exceeds(locs)

    def view_limits(self, dmin, dmax):
        """
        Set the view limits to the nearest multiples of base that
        contain the data
        """
        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            vmin = self._base.le(dmin)
            vmax = self._base.ge(dmax)
            if vmin == vmax:
                vmin -= 1
                vmax += 1
        else:
            vmin = dmin
            vmax = dmax

        return mtransforms.nonsingular(vmin, vmax)


def scale_range(vmin, vmax, n=1, threshold=100):
    dv = abs(vmax - vmin)  # > 0 as nonsingular is called before.
    meanv = (vmax + vmin) / 2
    if abs(meanv) / dv < threshold:
        offset = 0
    else:
        offset = math.copysign(10 ** (math.log10(abs(meanv)) // 1), meanv)
    scale = 10 ** (math.log10(dv / n) // 1)
    return scale, offset


class MaxNLocator(Locator):
    """
    Select no more than N intervals at nice locations.
    """
    default_params = dict(nbins=10,
                          steps=None,
                          integer=False,
                          symmetric=False,
                          prune=None,
                          min_n_ticks=2)

    def __init__(self, *args, **kwargs):
        """
        Keyword args:

        *nbins*
            Maximum number of intervals; one less than max number of
            ticks.  If the string `'auto'`, the number of bins will be
            automatically determined based on the length of the axis.

        *steps*
            Sequence of nice numbers starting with 1 and ending with 10;
            e.g., [1, 2, 4, 5, 10], where the values are acceptable
            tick multiples.  i.e. for the example, 20, 40, 60 would be
            an acceptable set of ticks, as would 0.4, 0.6, 0.8, because
            they are multiples of 2.  However, 30, 60, 90 would not
            be allowed because 3 does not appear in the list of steps.

        *integer*
            If True, ticks will take only integer values, provided
            at least `min_n_ticks` integers are found within the
            view limits.

        *symmetric*
            If True, autoscaling will result in a range symmetric
            about zero.

        *prune*
            ['lower' | 'upper' | 'both' | None]
            Remove edge ticks -- useful for stacked or ganged plots where
            the upper tick of one axes overlaps with the lower tick of the
            axes above it, primarily when :rc:`axes.autolimit_mode` is
            ``'round_numbers'``.  If ``prune=='lower'``, the smallest tick will
            be removed.  If ``prune == 'upper'``, the largest tick will be
            removed.  If ``prune == 'both'``, the largest and smallest ticks
            will be removed.  If ``prune == None``, no ticks will be removed.

        *min_n_ticks*
            Relax `nbins` and `integer` constraints if necessary to
            obtain this minimum number of ticks.

        """
        if args:
            kwargs['nbins'] = args[0]
            if len(args) > 1:
                raise ValueError(
                    "Keywords are required for all arguments except 'nbins'")
        self.set_params(**self.default_params)
        self.set_params(**kwargs)

    @staticmethod
    def _validate_steps(steps):
        if not np.iterable(steps):
            raise ValueError('steps argument must be a sequence of numbers '
                             'from 1 to 10')
        steps = np.asarray(steps)
        if np.any(np.diff(steps) <= 0):
            raise ValueError('steps argument must be uniformly increasing')
        if steps[-1] > 10 or steps[0] < 1:
            warnings.warn('Steps argument should be a sequence of numbers\n'
                          'increasing from 1 to 10, inclusive. Behavior with\n'
                          'values outside this range is undefined, and will\n'
                          'raise a ValueError in future versions of mpl.')
        if steps[0] != 1:
            steps = np.hstack((1, steps))
        if steps[-1] != 10:
            steps = np.hstack((steps, 10))
        return steps

    @staticmethod
    def _staircase(steps):
        # Make an extended staircase within which the needed
        # step will be found.  This is probably much larger
        # than necessary.
        flights = (0.1 * steps[:-1], steps, 10 * steps[1])
        return np.hstack(flights)

    def set_params(self, **kwargs):
        """Set parameters within this locator."""
        if 'nbins' in kwargs:
            self._nbins = kwargs['nbins']
            if self._nbins != 'auto':
                self._nbins = int(self._nbins)
        if 'symmetric' in kwargs:
            self._symmetric = kwargs['symmetric']
        if 'prune' in kwargs:
            prune = kwargs['prune']
            if prune is not None and prune not in ['upper', 'lower', 'both']:
                raise ValueError(
                    "prune must be 'upper', 'lower', 'both', or None")
            self._prune = prune
        if 'min_n_ticks' in kwargs:
            self._min_n_ticks = max(1, kwargs['min_n_ticks'])
        if 'steps' in kwargs:
            steps = kwargs['steps']
            if steps is None:
                self._steps = np.array([1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10])
            else:
                self._steps = self._validate_steps(steps)
            self._extended_steps = self._staircase(self._steps)
        if 'integer' in kwargs:
            self._integer = kwargs['integer']

    def _raw_ticks(self, vmin, vmax):
        if self._nbins == 'auto':
            if self.axis is not None:
                nbins = np.clip(self.axis.get_tick_space(),
                                max(1, self._min_n_ticks - 1), 9)
            else:
                nbins = 9
        else:
            nbins = self._nbins

        scale, offset = scale_range(vmin, vmax, nbins)
        _vmin = vmin - offset
        _vmax = vmax - offset
        raw_step = (vmax - vmin) / nbins
        steps = self._extended_steps * scale
        if self._integer:
            # For steps > 1, keep only integer values.
            igood = (steps < 1) | (np.abs(steps - np.round(steps)) < 0.001)
            steps = steps[igood]

        istep = np.nonzero(steps >= raw_step)[0][0]

        # Classic round_numbers mode may require a larger step.
        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            for istep in range(istep, len(steps)):
                step = steps[istep]
                best_vmin = (_vmin // step) * step
                best_vmax = best_vmin + step * nbins
                if (best_vmax >= _vmax):
                    break

        # This is an upper limit; move to smaller steps if necessary.
        for i in range(istep):
            step = steps[istep - i]
            if (self._integer and
                    np.floor(_vmax) - np.ceil(_vmin) >= self._min_n_ticks - 1):
                step = max(1, step)
            best_vmin = (_vmin // step) * step

            low = np.round(Base(step).le(_vmin - best_vmin) / step)
            high = np.round(Base(step).ge(_vmax - best_vmin) / step)
            ticks = np.arange(low, high + 1) * step + best_vmin + offset
            nticks = ((ticks <= vmax) & (ticks >= vmin)).sum()
            if nticks >= self._min_n_ticks:
                break
        return ticks

    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        if self._symmetric:
            vmax = max(abs(vmin), abs(vmax))
            vmin = -vmax
        vmin, vmax = mtransforms.nonsingular(
            vmin, vmax, expander=1e-13, tiny=1e-14)
        locs = self._raw_ticks(vmin, vmax)

        prune = self._prune
        if prune == 'lower':
            locs = locs[1:]
        elif prune == 'upper':
            locs = locs[:-1]
        elif prune == 'both':
            locs = locs[1:-1]
        return self.raise_if_exceeds(locs)

    def view_limits(self, dmin, dmax):
        if self._symmetric:
            dmax = max(abs(dmin), abs(dmax))
            dmin = -dmax

        dmin, dmax = mtransforms.nonsingular(
            dmin, dmax, expander=1e-12, tiny=1e-13)

        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            return self._raw_ticks(dmin, dmax)[[0, -1]]
        else:
            return dmin, dmax


def decade_down(x, base=10):
    'floor x to the nearest lower decade'
    if x == 0.0:
        return -base
    lx = np.floor(np.log(x) / np.log(base))
    return base ** lx


def decade_up(x, base=10):
    'ceil x to the nearest higher decade'
    if x == 0.0:
        return base
    lx = np.ceil(np.log(x) / np.log(base))
    return base ** lx


def nearest_long(x):
    if x == 0:
        return long(0)
    elif x > 0:
        return long(x + 0.5)
    else:
        return long(x - 0.5)


def is_decade(x, base=10):
    if not np.isfinite(x):
        return False
    if x == 0.0:
        return True
    lx = np.log(np.abs(x)) / np.log(base)
    return is_close_to_int(lx)


def is_close_to_int(x):
    if not np.isfinite(x):
        return False
    return abs(x - nearest_long(x)) < 1e-10


class LogLocator(Locator):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, subs=(1.0,), numdecs=4, numticks=None):
        """
        Place ticks on the locations : subs[j] * base**i

        Parameters
        ----------
        subs : None, string, or sequence of float, optional, default (1.0,)
            Gives the multiples of integer powers of the base at which
            to place ticks.  The default places ticks only at
            integer powers of the base.
            The permitted string values are ``'auto'`` and ``'all'``,
            both of which use an algorithm based on the axis view
            limits to determine whether and how to put ticks between
            integer powers of the base.  With ``'auto'``, ticks are
            placed only between integer powers; with ``'all'``, the
            integer powers are included.  A value of None is
            equivalent to ``'auto'``.

        """
        if numticks is None:
            if rcParams['_internal.classic_mode']:
                numticks = 15
            else:
                numticks = 'auto'
        self.base(base)
        self.subs(subs)
        self.numdecs = numdecs
        self.numticks = numticks

    def set_params(self, base=None, subs=None, numdecs=None, numticks=None):
        """Set parameters within this locator."""
        if base is not None:
            self.base(base)
        if subs is not None:
            self.subs(subs)
        if numdecs is not None:
            self.numdecs = numdecs
        if numticks is not None:
            self.numticks = numticks

    # FIXME: these base and subs functions are contrary to our
    # usual and desired API.

    def base(self, base):
        """
        set the base of the log scaling (major tick every base**i, i integer)
        """
        self._base = float(base)

    def subs(self, subs):
        """
        set the minor ticks for the log scaling every base**i*subs[j]
        """
        if subs is None:  # consistency with previous bad API
            self._subs = 'auto'
        elif isinstance(subs, six.string_types):
            if subs not in ('all', 'auto'):
                raise ValueError("A subs string must be 'all' or 'auto'; "
                                 "found '%s'." % subs)
            self._subs = subs
        else:
            self._subs = np.asarray(subs, dtype=float)

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        if self.numticks == 'auto':
            if self.axis is not None:
                numticks = np.clip(self.axis.get_tick_space(), 2, 9)
            else:
                numticks = 9
        else:
            numticks = self.numticks

        b = self._base
        # dummy axis has no axes attribute
        if hasattr(self.axis, 'axes') and self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            decades = np.arange(vmax - self.numdecs, vmax)
            ticklocs = b ** decades

            return ticklocs

        if vmin <= 0.0:
            if self.axis is not None:
                vmin = self.axis.get_minpos()

            if vmin <= 0.0 or not np.isfinite(vmin):
                raise ValueError(
                    "Data has no positive values, and therefore can not be "
                    "log-scaled.")

        vmin = math.log(vmin) / math.log(b)
        vmax = math.log(vmax) / math.log(b)

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        numdec = math.floor(vmax) - math.ceil(vmin)

        if isinstance(self._subs, six.string_types):
            _first = 2.0 if self._subs == 'auto' else 1.0
            if numdec > 10 or b < 3:
                if self._subs == 'auto':
                    return np.array([])  # no minor or major ticks
                else:
                    subs = np.array([1.0])  # major ticks
            else:
                subs = np.arange(_first, b)
        else:
            subs = self._subs

        stride = 1

        if rcParams['_internal.classic_mode']:
            # Leave the bug left over from the PY2-PY3 transition.
            while numdec / stride + 1 > numticks:
                stride += 1
        else:
            while numdec // stride + 1 > numticks:
                stride += 1

        # Does subs include anything other than 1?
        have_subs = len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0)

        decades = np.arange(math.floor(vmin) - stride,
                            math.ceil(vmax) + 2 * stride, stride)

        if hasattr(self, '_transform'):
            ticklocs = self._transform.inverted().transform(decades)
            if have_subs:
                if stride == 1:
                    ticklocs = np.ravel(np.outer(subs, ticklocs))
                else:
                    ticklocs = []
        else:
            if have_subs:
                ticklocs = []
                if stride == 1:
                    for decadeStart in b ** decades:
                        ticklocs.extend(subs * decadeStart)
            else:
                ticklocs = b ** decades

        return self.raise_if_exceeds(np.asarray(ticklocs))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
        b = self._base

        vmin, vmax = self.nonsingular(vmin, vmax)

        if self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            vmin = b ** (vmax - self.numdecs)

        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            if not is_decade(vmin, self._base):
                vmin = decade_down(vmin, self._base)
            if not is_decade(vmax, self._base):
                vmax = decade_up(vmax, self._base)

        return vmin, vmax

    def nonsingular(self, vmin, vmax):
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            return 1, 10  # initial range, no data plotted yet

        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if vmax <= 0:
            warnings.warn(
                "Data has no positive values, and therefore cannot be "
                "log-scaled.")
            return 1, 10

        minpos = self.axis.get_minpos()
        if not np.isfinite(minpos):
            minpos = 1e-300  # This should never take effect.
        if vmin <= 0:
            vmin = minpos
        if vmin == vmax:
            vmin = decade_down(vmin, self._base)
            vmax = decade_up(vmax, self._base)
        return vmin, vmax


class SymmetricalLogLocator(Locator):
    """
    Determine the tick locations for symmetric log axes
    """

    def __init__(self, transform=None, subs=None, linthresh=None, base=None):
        """
        place ticks on the location= base**i*subs[j]
        """
        if transform is not None:
            self._base = transform.base
            self._linthresh = transform.linthresh
        elif linthresh is not None and base is not None:
            self._base = base
            self._linthresh = linthresh
        else:
            raise ValueError("Either transform, or both linthresh "
                             "and base, must be provided.")
        if subs is None:
            self._subs = [1.0]
        else:
            self._subs = subs
        self.numticks = 15

    def set_params(self, subs=None, numticks=None):
        """Set parameters within this locator."""
        if numticks is not None:
            self.numticks = numticks
        if subs is not None:
            self._subs = subs

    def __call__(self):
        'Return the locations of the ticks'
        # Note, these are untransformed coordinates
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        b = self._base
        t = self._linthresh

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        # The domain is divided into three sections, only some of
        # which may actually be present.
        #
        # <======== -t ==0== t ========>
        # aaaaaaaaa    bbbbb   ccccccccc
        #
        # a) and c) will have ticks at integral log positions.  The
        # number of ticks needs to be reduced if there are more
        # than self.numticks of them.
        #
        # b) has a tick at 0 and only 0 (we assume t is a small
        # number, and the linear segment is just an implementation
        # detail and not interesting.)
        #
        # We could also add ticks at t, but that seems to usually be
        # uninteresting.
        #
        # "simple" mode is when the range falls entirely within (-t,
        # t) -- it should just display (vmin, 0, vmax)

        has_a = has_b = has_c = False
        if vmin < -t:
            has_a = True
            if vmax > -t:
                has_b = True
                if vmax > t:
                    has_c = True
        elif vmin < 0:
            if vmax > 0:
                has_b = True
                if vmax > t:
                    has_c = True
            else:
                return [vmin, vmax]
        elif vmin < t:
            if vmax > t:
                has_b = True
                has_c = True
            else:
                return [vmin, vmax]
        else:
            has_c = True

        def get_log_range(lo, hi):
            lo = np.floor(np.log(lo) / np.log(b))
            hi = np.ceil(np.log(hi) / np.log(b))
            return lo, hi

        # First, calculate all the ranges, so we can determine striding
        if has_a:
            if has_b:
                a_range = get_log_range(t, -vmin + 1)
            else:
                a_range = get_log_range(-vmax, -vmin + 1)
        else:
            a_range = (0, 0)

        if has_c:
            if has_b:
                c_range = get_log_range(t, vmax + 1)
            else:
                c_range = get_log_range(vmin, vmax + 1)
        else:
            c_range = (0, 0)

        total_ticks = (a_range[1] - a_range[0]) + (c_range[1] - c_range[0])
        if has_b:
            total_ticks += 1
        stride = max(total_ticks // (self.numticks - 1), 1)

        decades = []
        if has_a:
            decades.extend(-1 * (b ** (np.arange(a_range[0], a_range[1],
                                                 stride)[::-1])))

        if has_b:
            decades.append(0.0)

        if has_c:
            decades.extend(b ** (np.arange(c_range[0], c_range[1], stride)))

        # Add the subticks if requested
        if self._subs is None:
            subs = np.arange(2.0, b)
        else:
            subs = np.asarray(self._subs)

        if len(subs) > 1 or subs[0] != 1.0:
            ticklocs = []
            for decade in decades:
                if decade == 0:
                    ticklocs.append(decade)
                else:
                    ticklocs.extend(subs * decade)
        else:
            ticklocs = decades

        return self.raise_if_exceeds(np.array(ticklocs))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
        b = self._base
        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if rcParams['axes.autolimit_mode'] == 'round_numbers':
            if not is_decade(abs(vmin), b):
                if vmin < 0:
                    vmin = -decade_up(-vmin, b)
                else:
                    vmin = decade_down(vmin, b)
            if not is_decade(abs(vmax), b):
                if vmax < 0:
                    vmax = -decade_down(-vmax, b)
                else:
                    vmax = decade_up(vmax, b)

            if vmin == vmax:
                if vmin < 0:
                    vmin = -decade_up(-vmin, b)
                    vmax = -decade_down(-vmax, b)
                else:
                    vmin = decade_down(vmin, b)
                    vmax = decade_up(vmax, b)

        result = mtransforms.nonsingular(vmin, vmax)
        return result


class LogitLocator(Locator):
    """
    Determine the tick locations for logit axes
    """

    def __init__(self, minor=False):
        """
        place ticks on the logit locations
        """
        self.minor = minor

    def set_params(self, minor=None):
        """Set parameters within this locator."""
        if minor is not None:
            self.minor = minor

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        # dummy axis has no axes attribute
        if hasattr(self.axis, 'axes') and self.axis.axes.name == 'polar':
            raise NotImplementedError('Polar axis cannot be logit scaled yet')

        vmin, vmax = self.nonsingular(vmin, vmax)
        vmin = np.log10(vmin / (1 - vmin))
        vmax = np.log10(vmax / (1 - vmax))

        decade_min = np.floor(vmin)
        decade_max = np.ceil(vmax)

        # major ticks
        if not self.minor:
            ticklocs = []
            if (decade_min <= -1):
                expo = np.arange(decade_min, min(0, decade_max + 1))
                ticklocs.extend(list(10**expo))
            if (decade_min <= 0) and (decade_max >= 0):
                ticklocs.append(0.5)
            if (decade_max >= 1):
                expo = -np.arange(max(1, decade_min), decade_max + 1)
                ticklocs.extend(list(1 - 10**expo))

        # minor ticks
        else:
            ticklocs = []
            if (decade_min <= -2):
                expo = np.arange(decade_min, min(-1, decade_max))
                newticks = np.outer(np.arange(2, 10), 10**expo).ravel()
                ticklocs.extend(list(newticks))
            if (decade_min <= 0) and (decade_max >= 0):
                ticklocs.extend([0.2, 0.3, 0.4, 0.6, 0.7, 0.8])
            if (decade_max >= 2):
                expo = -np.arange(max(2, decade_min), decade_max + 1)
                newticks = 1 - np.outer(np.arange(2, 10), 10**expo).ravel()
                ticklocs.extend(list(newticks))

        return self.raise_if_exceeds(np.array(ticklocs))

    def nonsingular(self, vmin, vmax):
        initial_range = (1e-7, 1 - 1e-7)
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            return initial_range  # no data plotted yet

        if vmin > vmax:
            vmin, vmax = vmax, vmin

        # what to do if a window beyond ]0, 1[ is chosen
        if self.axis is not None:
            minpos = self.axis.get_minpos()
            if not np.isfinite(minpos):
                return initial_range  # again, no data plotted
        else:
            minpos = 1e-7  # should not occur in normal use

        # NOTE: for vmax, we should query a property similar to get_minpos, but
        # related to the maximal, less-than-one data point. Unfortunately,
        # Bbox._minpos is defined very deep in the BBox and updated with data,
        # so for now we use 1 - minpos as a substitute.

        if vmin <= 0:
            vmin = minpos
        if vmax >= 1:
            vmax = 1 - minpos
        if vmin == vmax:
            return 0.1 * vmin, 1 - 0.1 * vmin

        return vmin, vmax


class AutoLocator(MaxNLocator):
    """
    Dynamically find major tick positions. This is actually a subclass
    of `~matplotlib.ticker.MaxNLocator`, with parameters *nbins = 'auto'*
    and *steps = [1, 2, 2.5, 5, 10]*.
    """
    def __init__(self):
        """
        To know the values of the non-public parameters, please have a
        look to the defaults of `~matplotlib.ticker.MaxNLocator`.
        """
        if rcParams['_internal.classic_mode']:
            nbins = 9
            steps = [1, 2, 5, 10]
        else:
            nbins = 'auto'
            steps = [1, 2, 2.5, 5, 10]
        MaxNLocator.__init__(self, nbins=nbins, steps=steps)


class AutoMinorLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks. The scale must be linear with major ticks evenly spaced.
    """
    def __init__(self, n=None):
        """
        *n* is the number of subdivisions of the interval between
        major ticks; e.g., n=2 will place a single minor tick midway
        between major ticks.

        If *n* is omitted or None, it will be set to 5 or 4.
        """
        self.ndivs = n

    def __call__(self):
        'Return the locations of the ticks'
        if self.axis.get_scale() == 'log':
            warnings.warn('AutoMinorLocator does not work with logarithmic '
                          'scale')
            return []

        majorlocs = self.axis.get_majorticklocs()
        try:
            majorstep = majorlocs[1] - majorlocs[0]
        except IndexError:
            # Need at least two major ticks to find minor tick locations
            # TODO: Figure out a way to still be able to display minor
            # ticks without two major ticks visible. For now, just display
            # no ticks at all.
            return []

        if self.ndivs is None:
            x = int(np.round(10 ** (np.log10(majorstep) % 1)))
            if x in [1, 5, 10]:
                ndivs = 5
            else:
                ndivs = 4
        else:
            ndivs = self.ndivs

        minorstep = majorstep / ndivs

        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin, vmax = vmax, vmin

        t0 = majorlocs[0]
        tmin = ((vmin - t0) // minorstep + 1) * minorstep
        tmax = ((vmax - t0) // minorstep + 1) * minorstep
        locs = np.arange(tmin, tmax, minorstep) + t0
        cond = np.abs((locs - t0) % majorstep) > minorstep / 10.0
        locs = locs.compress(cond)

        return self.raise_if_exceeds(np.array(locs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


class OldAutoLocator(Locator):
    """
    On autoscale this class picks the best MultipleLocator to set the
    view limits and the tick locs.

    """
    def __init__(self):
        self._locator = LinearLocator()

    def __call__(self):
        'Return the locations of the ticks'
        self.refresh()
        return self.raise_if_exceeds(self._locator())

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))

    def refresh(self):
        'refresh internal information based on current lim'
        vmin, vmax = self.axis.get_view_interval()
        vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander=0.05)
        d = abs(vmax - vmin)
        self._locator = self.get_locator(d)

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'

        d = abs(vmax - vmin)
        self._locator = self.get_locator(d)
        return self._locator.view_limits(vmin, vmax)

    def get_locator(self, d):
        'pick the best locator based on a distance'
        d = abs(d)
        if d <= 0:
            locator = MultipleLocator(0.2)
        else:

            try:
                ld = math.log10(d)
            except OverflowError:
                raise RuntimeError('AutoLocator illegal data interval range')

            fld = math.floor(ld)
            base = 10 ** fld

            #if ld==fld:  base = 10**(fld-1)
            #else:        base = 10**fld

            if d >= 5 * base:
                ticksize = base
            elif d >= 2 * base:
                ticksize = base / 2.0
            else:
                ticksize = base / 5.0
            locator = MultipleLocator(ticksize)

        return locator
