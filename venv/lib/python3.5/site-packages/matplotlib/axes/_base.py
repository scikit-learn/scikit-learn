from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict

import six
from six.moves import xrange

import itertools
import warnings
import math
from operator import attrgetter

import numpy as np

import matplotlib

from matplotlib import cbook
from matplotlib.cbook import (_check_1d, _string_to_bool, iterable,
                              index_of, get_label)
from matplotlib import docstring
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.artist as martist
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker
import matplotlib.axis as maxis
import matplotlib.scale as mscale
import matplotlib.spines as mspines
import matplotlib.font_manager as font_manager
import matplotlib.text as mtext
import matplotlib.image as mimage
from matplotlib.offsetbox import OffsetBox
from matplotlib.artist import allow_rasterization
from matplotlib.legend import Legend

from matplotlib.rcsetup import cycler
from matplotlib.rcsetup import validate_axisbelow

rcParams = matplotlib.rcParams

is_string_like = cbook.is_string_like
is_sequence_of_strings = cbook.is_sequence_of_strings

_hold_msg = """axes.hold is deprecated.
    See the API Changes document (http://matplotlib.org/api/api_changes.html)
    for more details."""


def _process_plot_format(fmt):
    """
    Process a MATLAB style color/line style format string.  Return a
    (*linestyle*, *color*) tuple as a result of the processing.  Default
    values are ('-', 'b').  Example format strings include:

    * 'ko': black circles
    * '.b': blue dots
    * 'r--': red dashed lines
    * 'C2--': the third color in the color cycle, dashed lines

    .. seealso::

        :func:`~matplotlib.Line2D.lineStyles` and
        :func:`~matplotlib.pyplot.colors`
            for all possible styles and color format string.
    """

    linestyle = None
    marker = None
    color = None

    # Is fmt just a colorspec?
    try:
        color = mcolors.to_rgba(fmt)

        # We need to differentiate grayscale '1.0' from tri_down marker '1'
        try:
            fmtint = str(int(fmt))
        except ValueError:
            return linestyle, marker, color  # Yes
        else:
            if fmt != fmtint:
                # user definitely doesn't want tri_down marker
                return linestyle, marker, color  # Yes
            else:
                # ignore converted color
                color = None
    except ValueError:
        pass  # No, not just a color.

    # handle the multi char special cases and strip them from the
    # string
    if fmt.find('--') >= 0:
        linestyle = '--'
        fmt = fmt.replace('--', '')
    if fmt.find('-.') >= 0:
        linestyle = '-.'
        fmt = fmt.replace('-.', '')
    if fmt.find(' ') >= 0:
        linestyle = 'None'
        fmt = fmt.replace(' ', '')

    chars = [c for c in fmt]

    i = 0
    while i < len(chars):
        c = chars[i]
        if c in mlines.lineStyles:
            if linestyle is not None:
                raise ValueError(
                    'Illegal format string "%s"; two linestyle symbols' % fmt)
            linestyle = c
        elif c in mlines.lineMarkers:
            if marker is not None:
                raise ValueError(
                    'Illegal format string "%s"; two marker symbols' % fmt)
            marker = c
        elif c in mcolors.get_named_colors_mapping():
            if color is not None:
                raise ValueError(
                    'Illegal format string "%s"; two color symbols' % fmt)
            color = c
        elif c == 'C' and i < len(chars) - 1:
            color_cycle_number = int(chars[i + 1])
            color = mcolors.to_rgba("C{}".format(color_cycle_number))
            i += 1
        else:
            raise ValueError(
                'Unrecognized character %c in format string' % c)
        i += 1

    if linestyle is None and marker is None:
        linestyle = rcParams['lines.linestyle']
    if linestyle is None:
        linestyle = 'None'
    if marker is None:
        marker = 'None'

    return linestyle, marker, color


class _process_plot_var_args(object):
    """
    Process variable length arguments to the plot command, so that
    plot commands like the following are supported::

      plot(t, s)
      plot(t1, s1, t2, s2)
      plot(t1, s1, 'ko', t2, s2)
      plot(t1, s1, 'ko', t2, s2, 'r--', t3, e3)

    an arbitrary number of *x*, *y*, *fmt* are allowed
    """
    def __init__(self, axes, command='plot'):
        self.axes = axes
        self.command = command
        self.set_prop_cycle()

    def __getstate__(self):
        # note: it is not possible to pickle a itertools.cycle instance
        return {'axes': self.axes, 'command': self.command}

    def __setstate__(self, state):
        self.__dict__ = state.copy()
        self.set_prop_cycle()

    def set_prop_cycle(self, *args, **kwargs):
        if not (args or kwargs) or (len(args) == 1 and args[0] is None):
            prop_cycler = rcParams['axes.prop_cycle']
        else:
            prop_cycler = cycler(*args, **kwargs)

        self.prop_cycler = itertools.cycle(prop_cycler)
        # This should make a copy
        self._prop_keys = prop_cycler.keys

    def __call__(self, *args, **kwargs):
        if self.axes.xaxis is not None and self.axes.yaxis is not None:
            xunits = kwargs.pop('xunits', self.axes.xaxis.units)

            if self.axes.name == 'polar':
                xunits = kwargs.pop('thetaunits', xunits)

            yunits = kwargs.pop('yunits', self.axes.yaxis.units)

            if self.axes.name == 'polar':
                yunits = kwargs.pop('runits', yunits)

            if xunits != self.axes.xaxis.units:
                self.axes.xaxis.set_units(xunits)

            if yunits != self.axes.yaxis.units:
                self.axes.yaxis.set_units(yunits)

        ret = self._grab_next_args(*args, **kwargs)
        return ret

    def get_next_color(self):
        """Return the next color in the cycle."""
        if 'color' not in self._prop_keys:
            return 'k'
        return next(self.prop_cycler)['color']

    def set_lineprops(self, line, **kwargs):
        assert self.command == 'plot', 'set_lineprops only works with "plot"'
        line.set(**kwargs)

    def set_patchprops(self, fill_poly, **kwargs):
        assert self.command == 'fill', 'set_patchprops only works with "fill"'
        fill_poly.set(**kwargs)

    def _xy_from_xy(self, x, y):
        if self.axes.xaxis is not None and self.axes.yaxis is not None:
            bx = self.axes.xaxis.update_units(x)
            by = self.axes.yaxis.update_units(y)

            if self.command != 'plot':
                # the Line2D class can handle unitized data, with
                # support for post hoc unit changes etc.  Other mpl
                # artists, e.g., Polygon which _process_plot_var_args
                # also serves on calls to fill, cannot.  So this is a
                # hack to say: if you are not "plot", which is
                # creating Line2D, then convert the data now to
                # floats.  If you are plot, pass the raw data through
                # to Line2D which will handle the conversion.  So
                # polygons will not support post hoc conversions of
                # the unit type since they are not storing the orig
                # data.  Hopefully we can rationalize this at a later
                # date - JDH
                if bx:
                    x = self.axes.convert_xunits(x)
                if by:
                    y = self.axes.convert_yunits(y)

        # like asanyarray, but converts scalar to array, and doesn't change
        # existing compatible sequences
        x = _check_1d(x)
        y = _check_1d(y)
        if x.shape[0] != y.shape[0]:
            raise ValueError("x and y must have same first dimension, but "
                             "have shapes {} and {}".format(x.shape, y.shape))
        if x.ndim > 2 or y.ndim > 2:
            raise ValueError("x and y can be no greater than 2-D, but have "
                             "shapes {} and {}".format(x.shape, y.shape))

        if x.ndim == 1:
            x = x[:, np.newaxis]
        if y.ndim == 1:
            y = y[:, np.newaxis]
        return x, y

    def _getdefaults(self, ignore, *kwargs):
        """
        Only advance the cycler if the cycler has information that
        is not specified in any of the supplied tuple of dicts.
        Ignore any keys specified in the `ignore` set.

        Returns a copy of defaults dictionary if there are any
        keys that are not found in any of the supplied dictionaries.
        If the supplied dictionaries have non-None values for
        everything the property cycler has, then just return
        an empty dictionary. Ignored keys are excluded from the
        returned dictionary.

        """
        prop_keys = self._prop_keys
        if ignore is None:
            ignore = set()
        prop_keys = prop_keys - ignore

        if any(all(kw.get(k, None) is None for kw in kwargs)
               for k in prop_keys):
            # Need to copy this dictionary or else the next time around
            # in the cycle, the dictionary could be missing entries.
            default_dict = next(self.prop_cycler).copy()
            for p in ignore:
                default_dict.pop(p, None)
        else:
            default_dict = {}
        return default_dict

    def _setdefaults(self, defaults, *kwargs):
        """
        Given a defaults dictionary, and any other dictionaries,
        update those other dictionaries with information in defaults if
        none of the other dictionaries contains that information.

        """
        for k in defaults:
            if all(kw.get(k, None) is None for kw in kwargs):
                for kw in kwargs:
                    kw[k] = defaults[k]

    def _makeline(self, x, y, kw, kwargs):
        kw = kw.copy()  # Don't modify the original kw.
        kw.update(kwargs)
        default_dict = self._getdefaults(None, kw)
        self._setdefaults(default_dict, kw)
        seg = mlines.Line2D(x, y, **kw)
        return seg

    def _makefill(self, x, y, kw, kwargs):
        kw = kw.copy()  # Don't modify the original kw.
        kwargs = kwargs.copy()

        # Ignore 'marker'-related properties as they aren't Polygon
        # properties, but they are Line2D properties, and so they are
        # likely to appear in the default cycler construction.
        # This is done here to the defaults dictionary as opposed to the
        # other two dictionaries because we do want to capture when a
        # *user* explicitly specifies a marker which should be an error.
        # We also want to prevent advancing the cycler if there are no
        # defaults needed after ignoring the given properties.
        ignores = {'marker', 'markersize', 'markeredgecolor',
                   'markerfacecolor', 'markeredgewidth'}
        # Also ignore anything provided by *kwargs*.
        for k, v in six.iteritems(kwargs):
            if v is not None:
                ignores.add(k)

        # Only using the first dictionary to use as basis
        # for getting defaults for back-compat reasons.
        # Doing it with both seems to mess things up in
        # various places (probably due to logic bugs elsewhere).
        default_dict = self._getdefaults(ignores, kw)
        self._setdefaults(default_dict, kw)

        # Looks like we don't want "color" to be interpreted to
        # mean both facecolor and edgecolor for some reason.
        # So the "kw" dictionary is thrown out, and only its
        # 'color' value is kept and translated as a 'facecolor'.
        # This design should probably be revisited as it increases
        # complexity.
        facecolor = kw.get('color', None)

        # Throw out 'color' as it is now handled as a facecolor
        default_dict.pop('color', None)

        # To get other properties set from the cycler
        # modify the kwargs dictionary.
        self._setdefaults(default_dict, kwargs)

        seg = mpatches.Polygon(np.hstack((x[:, np.newaxis],
                                          y[:, np.newaxis])),
                               facecolor=facecolor,
                               fill=kwargs.get('fill', True),
                               closed=kw['closed'])
        self.set_patchprops(seg, **kwargs)
        return seg

    def _plot_args(self, tup, kwargs):
        ret = []
        if len(tup) > 1 and isinstance(tup[-1], six.string_types):
            linestyle, marker, color = _process_plot_format(tup[-1])
            tup = tup[:-1]
        elif len(tup) == 3:
            raise ValueError('third arg must be a format string')
        else:
            linestyle, marker, color = None, None, None

        # Don't allow any None value; These will be up-converted
        # to one element array of None which causes problems
        # downstream.
        if any(v is None for v in tup):
            raise ValueError("x and y must not be None")

        kw = {}
        for k, v in zip(('linestyle', 'marker', 'color'),
                        (linestyle, marker, color)):
            if v is not None:
                kw[k] = v

        if 'label' not in kwargs or kwargs['label'] is None:
            kwargs['label'] = get_label(tup[-1], None)

        if len(tup) == 2:
            x = _check_1d(tup[0])
            y = _check_1d(tup[-1])
        else:
            x, y = index_of(tup[-1])

        x, y = self._xy_from_xy(x, y)

        if self.command == 'plot':
            func = self._makeline
        else:
            kw['closed'] = kwargs.get('closed', True)
            func = self._makefill

        ncx, ncy = x.shape[1], y.shape[1]
        if ncx > 1 and ncy > 1 and ncx != ncy:
            cbook.warn_deprecated("2.2", "cycling among columns of inputs "
                                  "with non-matching shapes is deprecated.")
        for j in xrange(max(ncx, ncy)):
            seg = func(x[:, j % ncx], y[:, j % ncy], kw, kwargs)
            ret.append(seg)
        return ret

    def _grab_next_args(self, *args, **kwargs):
        while args:
            this, args = args[:2], args[2:]
            if args and isinstance(args[0], six.string_types):
                this += args[0],
                args = args[1:]
            for seg in self._plot_args(this, kwargs):
                yield seg


class _AxesBase(martist.Artist):
    """
    """
    name = "rectilinear"

    _shared_x_axes = cbook.Grouper()
    _shared_y_axes = cbook.Grouper()
    _twinned_axes = cbook.Grouper()

    def __str__(self):
        return "{0}({1[0]:g},{1[1]:g};{1[2]:g}x{1[3]:g})".format(
            type(self).__name__, self._position.bounds)

    def __init__(self, fig, rect,
                 facecolor=None,  # defaults to rc axes.facecolor
                 frameon=True,
                 sharex=None,  # use Axes instance's xaxis info
                 sharey=None,  # use Axes instance's yaxis info
                 label='',
                 xscale=None,
                 yscale=None,
                 **kwargs
                 ):
        """
        Build an :class:`Axes` instance in
        :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* in
        :class:`~matplotlib.figure.Figure` coordinates

        Optional keyword arguments:

          ================   =========================================
          Keyword            Description
          ================   =========================================
          *adjustable*       [ 'box' | 'datalim' ]
          *alpha*            float: the alpha transparency (can be None)
          *anchor*           [ 'C', 'SW', 'S', 'SE', 'E', 'NE', 'N',
                               'NW', 'W' ]
          *aspect*           [ 'auto' | 'equal' | aspect_ratio ]
          *autoscale_on*     bool; whether to autoscale the *viewlim*
          *axisbelow*        [ bool | 'line' ] draw the grids
                             and ticks below or above most other artists,
                             or below lines but above patches
          *cursor_props*     a (*float*, *color*) tuple
          *figure*           a :class:`~matplotlib.figure.Figure`
                             instance
          *frame_on*         bool; whether to draw the axes frame
          *label*            the axes label
          *navigate*         bool
          *navigate_mode*    [ 'PAN' | 'ZOOM' | None ] the navigation
                             toolbar button status
          *position*         [left, bottom, width, height] in
                             class:`~matplotlib.figure.Figure` coords
          *sharex*           an class:`~matplotlib.axes.Axes` instance
                             to share the x-axis with
          *sharey*           an class:`~matplotlib.axes.Axes` instance
                             to share the y-axis with
          *title*            the title string
          *visible*          bool, whether the axes is visible
          *xlabel*           the xlabel
          *xlim*             (*xmin*, *xmax*) view limits
          *xscale*           [%(scale)s]
          *xticklabels*      sequence of strings
          *xticks*           sequence of floats
          *ylabel*           the ylabel strings
          *ylim*             (*ymin*, *ymax*) view limits
          *yscale*           [%(scale)s]
          *yticklabels*      sequence of strings
          *yticks*           sequence of floats
          ================   =========================================
        """ % {'scale': ' | '.join(
            [repr(x) for x in mscale.get_scale_names()])}
        martist.Artist.__init__(self)
        if isinstance(rect, mtransforms.Bbox):
            self._position = rect
        else:
            self._position = mtransforms.Bbox.from_bounds(*rect)
        if self._position.width < 0 or self._position.height < 0:
            raise ValueError('Width and height specified must be non-negative')
        self._originalPosition = self._position.frozen()
        # self.set_axes(self)
        self.axes = self
        self._aspect = 'auto'
        self._adjustable = 'box'
        self._anchor = 'C'
        self._sharex = sharex
        self._sharey = sharey
        if sharex is not None:
            self._shared_x_axes.join(self, sharex)
        if sharey is not None:
            self._shared_y_axes.join(self, sharey)
        self.set_label(label)
        self.set_figure(fig)

        self.set_axes_locator(kwargs.get("axes_locator", None))

        self.spines = self._gen_axes_spines()

        # this call may differ for non-sep axes, e.g., polar
        self._init_axis()
        if facecolor is None:
            facecolor = rcParams['axes.facecolor']
        self._facecolor = facecolor
        self._frameon = frameon
        self._axisbelow = rcParams['axes.axisbelow']

        self._rasterization_zorder = None

        self._hold = rcParams['axes.hold']
        if self._hold is None:
            self._hold = True

        self._connected = {}  # a dict from events to (id, func)
        self.cla()

        # funcs used to format x and y - fall back on major formatters
        self.fmt_xdata = None
        self.fmt_ydata = None

        self._cachedRenderer = None
        self.set_navigate(True)
        self.set_navigate_mode(None)

        if xscale:
            self.set_xscale(xscale)
        if yscale:
            self.set_yscale(yscale)

        if len(kwargs):
            self.update(kwargs)

        if self.xaxis is not None:
            self._xcid = self.xaxis.callbacks.connect(
                'units finalize', lambda: self._on_units_changed(scalex=True))

        if self.yaxis is not None:
            self._ycid = self.yaxis.callbacks.connect(
                'units finalize', lambda: self._on_units_changed(scaley=True))

        self.tick_params(
            top=rcParams['xtick.top'] and rcParams['xtick.minor.top'],
            bottom=rcParams['xtick.bottom'] and rcParams['xtick.minor.bottom'],
            labeltop=(rcParams['xtick.labeltop'] and
                      rcParams['xtick.minor.top']),
            labelbottom=(rcParams['xtick.labelbottom'] and
                         rcParams['xtick.minor.bottom']),
            left=rcParams['ytick.left'] and rcParams['ytick.minor.left'],
            right=rcParams['ytick.right'] and rcParams['ytick.minor.right'],
            labelleft=(rcParams['ytick.labelleft'] and
                       rcParams['ytick.minor.left']),
            labelright=(rcParams['ytick.labelright'] and
                        rcParams['ytick.minor.right']),
            which='minor')

        self.tick_params(
            top=rcParams['xtick.top'] and rcParams['xtick.major.top'],
            bottom=rcParams['xtick.bottom'] and rcParams['xtick.major.bottom'],
            labeltop=(rcParams['xtick.labeltop'] and
                      rcParams['xtick.major.top']),
            labelbottom=(rcParams['xtick.labelbottom'] and
                         rcParams['xtick.major.bottom']),
            left=rcParams['ytick.left'] and rcParams['ytick.major.left'],
            right=rcParams['ytick.right'] and rcParams['ytick.major.right'],
            labelleft=(rcParams['ytick.labelleft'] and
                       rcParams['ytick.major.left']),
            labelright=(rcParams['ytick.labelright'] and
                        rcParams['ytick.major.right']),
            which='major')

        self._layoutbox = None
        self._poslayoutbox = None

    def __getstate__(self):
        # The renderer should be re-created by the figure, and then cached at
        # that point.
        state = super(_AxesBase, self).__getstate__()
        state['_cachedRenderer'] = None
        state.pop('_layoutbox')
        state.pop('_poslayoutbox')

        return state

    def __setstate__(self, state):
        self.__dict__ = state
        # put the _remove_method back on all artists contained within the axes
        for container_name in ['lines', 'collections', 'tables', 'patches',
                               'texts', 'images']:
            container = getattr(self, container_name)
            for artist in container:
                artist._remove_method = container.remove
        self._stale = True
        self._layoutbox = None
        self._poslayoutbox = None

    def get_window_extent(self, *args, **kwargs):
        """
        get the axes bounding box in display space; *args* and
        *kwargs* are empty
        """
        bbox = self.bbox
        x_pad = self.xaxis.get_tick_padding()
        y_pad = self.yaxis.get_tick_padding()
        return mtransforms.Bbox([[bbox.x0 - x_pad, bbox.y0 - y_pad],
                                 [bbox.x1 + x_pad, bbox.y1 + y_pad]])

    def _init_axis(self):
        "move this out of __init__ because non-separable axes don't use it"
        self.xaxis = maxis.XAxis(self)
        self.spines['bottom'].register_axis(self.xaxis)
        self.spines['top'].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines['left'].register_axis(self.yaxis)
        self.spines['right'].register_axis(self.yaxis)
        self._update_transScale()

    def set_figure(self, fig):
        """
        Set the `.Figure` for this `.Axes`.

        .. ACCEPTS: `.Figure`

        Parameters
        ----------
        fig : `.Figure`
        """
        martist.Artist.set_figure(self, fig)

        self.bbox = mtransforms.TransformedBbox(self._position,
                                                fig.transFigure)
        # these will be updated later as data is added
        self.dataLim = mtransforms.Bbox.null()
        self.viewLim = mtransforms.Bbox.unit()
        self.transScale = mtransforms.TransformWrapper(
            mtransforms.IdentityTransform())

        self._set_lim_and_transforms()

    def _set_lim_and_transforms(self):
        """
        set the *_xaxis_transform*, *_yaxis_transform*,
        *transScale*, *transData*, *transLimits* and *transAxes*
        transformations.

        .. note::

            This method is primarily used by rectilinear projections
            of the :class:`~matplotlib.axes.Axes` class, and is meant
            to be overridden by new kinds of projection axes that need
            different transformations and limits. (See
            :class:`~matplotlib.projections.polar.PolarAxes` for an
            example.

        """
        self.transAxes = mtransforms.BboxTransformTo(self.bbox)

        # Transforms the x and y axis separately by a scale factor.
        # It is assumed that this part will have non-linear components
        # (e.g., for a log scale).
        self.transScale = mtransforms.TransformWrapper(
            mtransforms.IdentityTransform())

        # An affine transformation on the data, generally to limit the
        # range of the axes
        self.transLimits = mtransforms.BboxTransformFrom(
            mtransforms.TransformedBbox(self.viewLim, self.transScale))

        # The parentheses are important for efficiency here -- they
        # group the last two (which are usually affines) separately
        # from the first (which, with log-scaling can be non-affine).
        self.transData = self.transScale + (self.transLimits + self.transAxes)

        self._xaxis_transform = mtransforms.blended_transform_factory(
            self.transData, self.transAxes)
        self._yaxis_transform = mtransforms.blended_transform_factory(
            self.transAxes, self.transData)

    def get_xaxis_transform(self, which='grid'):
        """
        Get the transformation used for drawing x-axis labels, ticks
        and gridlines.  The x-direction is in data coordinates and the
        y-direction is in axis coordinates.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        if which == 'grid':
            return self._xaxis_transform
        elif which == 'tick1':
            # for cartesian projection, this is bottom spine
            return self.spines['bottom'].get_spine_transform()
        elif which == 'tick2':
            # for cartesian projection, this is top spine
            return self.spines['top'].get_spine_transform()
        else:
            raise ValueError('unknown value for which')

    def get_xaxis_text1_transform(self, pad_points):
        """
        Get the transformation used for drawing x-axis labels, which
        will add the given amount of padding (in points) between the
        axes and the label.  The x-direction is in data coordinates
        and the y-direction is in axis coordinates.  Returns a
        3-tuple of the form::

          (transform, valign, halign)

        where *valign* and *halign* are requested alignments for the
        text.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        labels_align = matplotlib.rcParams["xtick.alignment"]

        return (self.get_xaxis_transform(which='tick1') +
                mtransforms.ScaledTranslation(0, -1 * pad_points / 72.0,
                                              self.figure.dpi_scale_trans),
                "top", labels_align)

    def get_xaxis_text2_transform(self, pad_points):
        """
        Get the transformation used for drawing the secondary x-axis
        labels, which will add the given amount of padding (in points)
        between the axes and the label.  The x-direction is in data
        coordinates and the y-direction is in axis coordinates.
        Returns a 3-tuple of the form::

          (transform, valign, halign)

        where *valign* and *halign* are requested alignments for the
        text.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        labels_align = matplotlib.rcParams["xtick.alignment"]
        return (self.get_xaxis_transform(which='tick2') +
                mtransforms.ScaledTranslation(0, pad_points / 72.0,
                                              self.figure.dpi_scale_trans),
                "bottom", labels_align)

    def get_yaxis_transform(self, which='grid'):
        """
        Get the transformation used for drawing y-axis labels, ticks
        and gridlines.  The x-direction is in axis coordinates and the
        y-direction is in data coordinates.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        if which == 'grid':
            return self._yaxis_transform
        elif which == 'tick1':
            # for cartesian projection, this is bottom spine
            return self.spines['left'].get_spine_transform()
        elif which == 'tick2':
            # for cartesian projection, this is top spine
            return self.spines['right'].get_spine_transform()
        else:
            raise ValueError('unknown value for which')

    def get_yaxis_text1_transform(self, pad_points):
        """
        Get the transformation used for drawing y-axis labels, which
        will add the given amount of padding (in points) between the
        axes and the label.  The x-direction is in axis coordinates
        and the y-direction is in data coordinates.  Returns a 3-tuple
        of the form::

          (transform, valign, halign)

        where *valign* and *halign* are requested alignments for the
        text.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        labels_align = matplotlib.rcParams["ytick.alignment"]
        return (self.get_yaxis_transform(which='tick1') +
                mtransforms.ScaledTranslation(-1 * pad_points / 72.0, 0,
                                              self.figure.dpi_scale_trans),
                labels_align, "right")

    def get_yaxis_text2_transform(self, pad_points):
        """
        Get the transformation used for drawing the secondary y-axis
        labels, which will add the given amount of padding (in points)
        between the axes and the label.  The x-direction is in axis
        coordinates and the y-direction is in data coordinates.
        Returns a 3-tuple of the form::

          (transform, valign, halign)

        where *valign* and *halign* are requested alignments for the
        text.

        .. note::

            This transformation is primarily used by the
            :class:`~matplotlib.axis.Axis` class, and is meant to be
            overridden by new kinds of projections that may need to
            place axis elements in different locations.

        """
        labels_align = matplotlib.rcParams["ytick.alignment"]

        return (self.get_yaxis_transform(which='tick2') +
                mtransforms.ScaledTranslation(pad_points / 72.0, 0,
                                              self.figure.dpi_scale_trans),
                labels_align, "left")

    def _update_transScale(self):
        self.transScale.set(
            mtransforms.blended_transform_factory(
                self.xaxis.get_transform(), self.yaxis.get_transform()))
        if hasattr(self, "lines"):
            for line in self.lines:
                try:
                    line._transformed_path.invalidate()
                except AttributeError:
                    pass

    def get_position(self, original=False):
        """
        Get a copy of the axes rectangle as a `.Bbox`.

        Parameters
        ----------
        original : bool
            If ``True``, return the original position. Otherwise return the
            active position. For an explanation of the positions see
            `.set_position`.

        Returns
        -------
        pos : `.Bbox`

        """
        if original:
            return self._originalPosition.frozen()
        else:
            return self._position.frozen()

    def set_position(self, pos, which='both'):
        """
        Set the axes position.

        Axes have two position attributes. The 'original' position is the
        position allocated for the Axes. The 'active' position is the
        position the Axes is actually drawn at. These positions are usually
        the same unless a fixed aspect is set to the Axes. See `.set_aspect`
        for details.

        Parameters
        ----------
        pos : [left, bottom, width, height] or `~matplotlib.transforms.Bbox`
            The new position of the in `.Figure` coordinates.

        which : ['both' | 'active' | 'original'], optional
            Determines which position variables to change.

        """
        self._set_position(pos, which='both')
        # because this is being called externally to the library we
        # zero the constrained layout parts.
        self._layoutbox = None
        self._poslayoutbox = None

    def _set_position(self, pos, which='both'):
        """
        private version of set_position.  Call this internally
        to get the same functionality of `get_position`, but not
        to take the axis out of the constrained_layout
        hierarchy.
        """
        if not isinstance(pos, mtransforms.BboxBase):
            pos = mtransforms.Bbox.from_bounds(*pos)
        for ax in self._twinned_axes.get_siblings(self):
            if which in ('both', 'active'):
                ax._position.set(pos)
            if which in ('both', 'original'):
                ax._originalPosition.set(pos)
        self.stale = True

    def reset_position(self):
        """
        Reset the active position to the original position.

        This resets the a possible position change due to aspect constraints.
        For an explanation of the positions see `.set_position`.
        """
        for ax in self._twinned_axes.get_siblings(self):
            pos = ax.get_position(original=True)
            ax.set_position(pos, which='active')

    def set_axes_locator(self, locator):
        """
        Set the axes locator.

        .. ACCEPTS: a callable object which takes an axes instance and
           renderer and returns a bbox.

        Parameters
        ----------
        locator : callable
            A locator function, which takes an axes and a renderer and returns
            a bbox.
        """
        self._axes_locator = locator
        self.stale = True

    def get_axes_locator(self):
        """
        Return the axes_locator.
        """
        return self._axes_locator

    def _set_artist_props(self, a):
        """set the boilerplate props for artists added to axes"""
        a.set_figure(self.figure)
        if not a.is_transform_set():
            a.set_transform(self.transData)

        a.axes = self
        if a.mouseover:
            self.mouseover_set.add(a)

    def _gen_axes_patch(self):
        """
        Returns the patch used to draw the background of the axes.  It
        is also used as the clipping path for any data elements on the
        axes.

        In the standard axes, this is a rectangle, but in other
        projections it may not be.

        .. note::

            Intended to be overridden by new projection types.

        """
        return mpatches.Rectangle((0.0, 0.0), 1.0, 1.0)

    def _gen_axes_spines(self, locations=None, offset=0.0, units='inches'):
        """
        Returns a dict whose keys are spine names and values are
        Line2D or Patch instances. Each element is used to draw a
        spine of the axes.

        In the standard axes, this is a single line segment, but in
        other projections it may not be.

        .. note::

            Intended to be overridden by new projection types.

        """
        return OrderedDict([
            ('left', mspines.Spine.linear_spine(self, 'left')),
            ('right', mspines.Spine.linear_spine(self, 'right')),
            ('bottom', mspines.Spine.linear_spine(self, 'bottom')),
            ('top', mspines.Spine.linear_spine(self, 'top'))])

    def cla(self):
        """Clear the current axes."""
        # Note: this is called by Axes.__init__()

        # stash the current visibility state
        if hasattr(self, 'patch'):
            patch_visible = self.patch.get_visible()
        else:
            patch_visible = True

        xaxis_visible = self.xaxis.get_visible()
        yaxis_visible = self.yaxis.get_visible()

        self.xaxis.cla()
        self.yaxis.cla()

        for name, spine in six.iteritems(self.spines):
            spine.cla()

        self.ignore_existing_data_limits = True
        self.callbacks = cbook.CallbackRegistry()

        if self._sharex is not None:
            # major and minor are axis.Ticker class instances with
            # locator and formatter attributes
            self.xaxis.major = self._sharex.xaxis.major
            self.xaxis.minor = self._sharex.xaxis.minor
            x0, x1 = self._sharex.get_xlim()
            self.set_xlim(x0, x1, emit=False, auto=None)
            self.xaxis._scale = mscale.scale_factory(
                    self._sharex.xaxis.get_scale(), self.xaxis)
        else:
            self.xaxis._set_scale('linear')
            try:
                self.set_xlim(0, 1)
            except TypeError:
                pass

        if self._sharey is not None:
            self.yaxis.major = self._sharey.yaxis.major
            self.yaxis.minor = self._sharey.yaxis.minor
            y0, y1 = self._sharey.get_ylim()
            self.set_ylim(y0, y1, emit=False, auto=None)
            self.yaxis._scale = mscale.scale_factory(
                    self._sharey.yaxis.get_scale(), self.yaxis)
        else:
            self.yaxis._set_scale('linear')
            try:
                self.set_ylim(0, 1)
            except TypeError:
                pass
        # update the minor locator for x and y axis based on rcParams
        if (rcParams['xtick.minor.visible']):
            self.xaxis.set_minor_locator(mticker.AutoMinorLocator())

        if (rcParams['ytick.minor.visible']):
            self.yaxis.set_minor_locator(mticker.AutoMinorLocator())

        self._autoscaleXon = True
        self._autoscaleYon = True
        self._xmargin = rcParams['axes.xmargin']
        self._ymargin = rcParams['axes.ymargin']
        self._tight = None
        self._use_sticky_edges = True
        self._update_transScale()  # needed?

        self._get_lines = _process_plot_var_args(self)
        self._get_patches_for_fill = _process_plot_var_args(self, 'fill')

        self._gridOn = rcParams['axes.grid']
        self.lines = []
        self.patches = []
        self.texts = []
        self.tables = []
        self.artists = []
        self.images = []
        self.mouseover_set = set()
        self._current_image = None  # strictly for pyplot via _sci, _gci
        self.legend_ = None
        self.collections = []  # collection.Collection instances
        self.containers = []

        self.grid(False)  # Disable grid on init to use rcParameter
        self.grid(self._gridOn, which=rcParams['axes.grid.which'],
                  axis=rcParams['axes.grid.axis'])
        props = font_manager.FontProperties(
            size=rcParams['axes.titlesize'],
            weight=rcParams['axes.titleweight'])

        self.title = mtext.Text(
            x=0.5, y=1.0, text='',
            fontproperties=props,
            verticalalignment='baseline',
            horizontalalignment='center',
            )
        self._left_title = mtext.Text(
            x=0.0, y=1.0, text='',
            fontproperties=props.copy(),
            verticalalignment='baseline',
            horizontalalignment='left', )
        self._right_title = mtext.Text(
            x=1.0, y=1.0, text='',
            fontproperties=props.copy(),
            verticalalignment='baseline',
            horizontalalignment='right',
            )
        title_offset_points = rcParams['axes.titlepad']
        # refactor this out so it can be called in ax.set_title if
        # pad argument used...
        self._set_title_offset_trans(title_offset_points)

        for _title in (self.title, self._left_title, self._right_title):
            self._set_artist_props(_title)

        # The patch draws the background of the axes.  We want this to be below
        # the other artists.  We use the frame to draw the edges so we are
        # setting the edgecolor to None.
        self.patch = self._gen_axes_patch()
        self.patch.set_figure(self.figure)
        self.patch.set_facecolor(self._facecolor)
        self.patch.set_edgecolor('None')
        self.patch.set_linewidth(0)
        self.patch.set_transform(self.transAxes)

        self.set_axis_on()

        self.xaxis.set_clip_path(self.patch)
        self.yaxis.set_clip_path(self.patch)

        self._shared_x_axes.clean()
        self._shared_y_axes.clean()
        if self._sharex:
            self.xaxis.set_visible(xaxis_visible)
            self.patch.set_visible(patch_visible)

        if self._sharey:
            self.yaxis.set_visible(yaxis_visible)
            self.patch.set_visible(patch_visible)

        self.stale = True

    @property
    @cbook.deprecated("2.1", alternative="Axes.patch")
    def axesPatch(self):
        return self.patch

    def clear(self):
        """Clear the axes."""
        self.cla()

    def get_facecolor(self):
        """Get the Axes facecolor."""
        return self.patch.get_facecolor()
    get_fc = get_facecolor

    def set_facecolor(self, color):
        """Set the Axes facecolor.

        .. ACCEPTS: color

        Parameters
        ----------
        color : color
        """
        self._facecolor = color
        return self.patch.set_facecolor(color)
    set_fc = set_facecolor

    def _set_title_offset_trans(self, title_offset_points):
        """
        Set the offset for the title either from rcParams['axes.titlepad']
        or from set_title kwarg ``pad``.
        """
        self.titleOffsetTrans = mtransforms.ScaledTranslation(
                0.0, title_offset_points / 72.0,
                self.figure.dpi_scale_trans)
        for _title in (self.title, self._left_title, self._right_title):
            _title.set_transform(self.transAxes + self.titleOffsetTrans)
            _title.set_clip_box(None)

    def set_prop_cycle(self, *args, **kwargs):
        """
        Set the property cycle for any future plot commands on this Axes.

        set_prop_cycle(arg)
        set_prop_cycle(label, itr)
        set_prop_cycle(label1=itr1[, label2=itr2[, ...]])

        Form 1 simply sets given `Cycler` object.

        Form 2 creates and sets  a `Cycler` from a label and an iterable.

        Form 3 composes and sets  a `Cycler` as an inner product of the
        pairs of keyword arguments. In other words, all of the
        iterables are cycled simultaneously, as if through zip().

        Parameters
        ----------
        arg : Cycler
            Set the given Cycler.
            Can also be `None` to reset to the cycle defined by the
            current style.

        label : str
            The property key. Must be a valid `Artist` property.
            For example, 'color' or 'linestyle'. Aliases are allowed,
            such as 'c' for 'color' and 'lw' for 'linewidth'.

        itr : iterable
            Finite-length iterable of the property values. These values
            are validated and will raise a ValueError if invalid.

        See Also
        --------
            :func:`cycler`      Convenience function for creating your
                                own cyclers.

        """
        if args and kwargs:
            raise TypeError("Cannot supply both positional and keyword "
                            "arguments to this method.")
        if len(args) == 1 and args[0] is None:
            prop_cycle = None
        else:
            prop_cycle = cycler(*args, **kwargs)
        self._get_lines.set_prop_cycle(prop_cycle)
        self._get_patches_for_fill.set_prop_cycle(prop_cycle)

    @cbook.deprecated('1.5', alternative='`.set_prop_cycle`')
    def set_color_cycle(self, clist):
        """
        Set the color cycle for any future plot commands on this Axes.

        Parameters
        ----------
        clist
            A list of mpl color specifiers.
        """
        if clist is None:
            # Calling set_color_cycle() or set_prop_cycle() with None
            # effectively resets the cycle, but you can't do
            # set_prop_cycle('color', None). So we are special-casing this.
            self.set_prop_cycle(None)
        else:
            self.set_prop_cycle('color', clist)

    @cbook.deprecated("2.0")
    def ishold(self):
        """return the HOLD status of the axes

        The `hold` mechanism is deprecated and will be removed in
        v3.0.
        """

        return self._hold

    @cbook.deprecated("2.0", message=_hold_msg)
    def hold(self, b=None):
        """
        Set the hold state.

        The ``hold`` mechanism is deprecated and will be removed in
        v3.0.  The behavior will remain consistent with the
        long-time default value of True.

        If *hold* is *None* (default), toggle the *hold* state.  Else
        set the *hold* state to boolean value *b*.

        Examples::

          # toggle hold
          hold()

          # turn hold on
          hold(True)

          # turn hold off
          hold(False)

        When hold is *True*, subsequent plot commands will be added to
        the current axes.  When hold is *False*, the current axes and
        figure will be cleared on the next plot command

        """
        if b is None:
            self._hold = not self._hold
        else:
            self._hold = b

    def get_aspect(self):
        return self._aspect

    def set_aspect(self, aspect, adjustable=None, anchor=None, share=False):
        """
        Set the aspect of the axis scaling, i.e. the ratio of y-unit to x-unit.

        Parameters
        ----------
        aspect : ['auto' | 'equal'] or num
            Possible values:

            ========   ================================================
            value      description
            ========   ================================================
            'auto'     automatic; fill the position rectangle with data
            'equal'    same scaling from data to plot units for x and y
             num       a circle will be stretched such that the height
                       is num times the width. aspect=1 is the same as
                       aspect='equal'.
            ========   ================================================

        adjustable : None or ['box' | 'datalim'], optional
            If not ``None``, this defines which parameter will be adjusted to
            meet the required aspect. See `.set_adjustable` for further
            details.

        anchor : None or str or 2-tuple of float, optional
            If not ``None``, this defines where the Axes will be drawn if there
            is extra space due to aspect constraints. The most common way to
            to specify the anchor are abbreviations of cardinal directions:

            =====   =====================
            value   description
            =====   =====================
            'C'     centered
            'SW'    lower left corner
            'S'     middle of bottom edge
            'SE'    lower right corner
            etc.
            =====   =====================

            See `.set_anchor` for further details.

        share : bool, optional
            If ``True``, apply the settings to all shared Axes.
            Default is ``False``.

        See Also
        --------
        matplotlib.axes.Axes.set_adjustable
            defining the parameter to adjust in order to meet the required
            aspect.
        matplotlib.axes.Axes.set_anchor
            defining the position in case of extra space.
        """
        if not (isinstance(aspect, six.string_types)
                and aspect in ('equal', 'auto')):
            aspect = float(aspect)  # raise ValueError if necessary
        if share:
            axes = set(self._shared_x_axes.get_siblings(self)
                       + self._shared_y_axes.get_siblings(self))
        else:
            axes = [self]
        for ax in axes:
            ax._aspect = aspect

        if adjustable is None:
            adjustable = self._adjustable
        self.set_adjustable(adjustable, share=share)  # Handle sharing.

        if anchor is not None:
            self.set_anchor(anchor, share=share)
        self.stale = True

    def get_adjustable(self):
        return self._adjustable

    def set_adjustable(self, adjustable, share=False):
        """
        Define which parameter the Axes will change to achieve a given aspect.

        Parameters
        ----------
        adjustable : ['box' | 'datalim']
            If 'box', change the physical dimensions of the Axes.
            If 'datalim', change the ``x`` or ``y`` data limits.

        share : bool, optional
            If ``True``, apply the settings to all shared Axes.
            Default is ``False``.

        .. ACCEPTS: [ 'box' | 'datalim']

        See Also
        --------
        matplotlib.axes.Axes.set_aspect
            for a description of aspect handling.

        Notes
        -----
        Shared Axes (of which twinned Axes are a special case)
        impose restrictions on how aspect ratios can be imposed.
        For twinned Axes, use 'datalim'.  For Axes that share both
        x and y, use 'box'.  Otherwise, either 'datalim' or 'box'
        may be used.  These limitations are partly a requirement
        to avoid over-specification, and partly a result of the
        particular implementation we are currently using, in
        which the adjustments for aspect ratios are done sequentially
        and independently on each Axes as it is drawn.
        """
        if adjustable == 'box-forced':
            warnings.warn("The 'box-forced' keyword argument is deprecated"
                          " since 2.2.", cbook.mplDeprecation)
        if adjustable not in ('box', 'datalim', 'box-forced'):
            raise ValueError("argument must be 'box', or 'datalim'")
        if share:
            axes = set(self._shared_x_axes.get_siblings(self)
                       + self._shared_y_axes.get_siblings(self))
        else:
            axes = [self]
        for ax in axes:
            ax._adjustable = adjustable
        self.stale = True

    def get_anchor(self):
        """
        Get the anchor location.

        See Also
        --------
        matplotlib.axes.Axes.set_anchor
            for a description of the anchor.
        matplotlib.axes.Axes.set_aspect
            for a description of aspect handling.
        """
        return self._anchor

    def set_anchor(self, anchor, share=False):
        """
        Define the anchor location.

        The actual drawing area (active position) of the Axes may be smaller
        than the Bbox (original position) when a fixed aspect is required. The
        anchor defines where the drawing area will be located within the
        available space.

        .. ACCEPTS: [ 'C' | 'SW' | 'S' | 'SE' | 'E' | 'NE' | 'N' | 'NW' | 'W' ]

        Parameters
        ----------
        anchor : str or 2-tuple of floats
            The anchor position may be either:

            - a sequence (*cx*, *cy*). *cx* and *cy* may range from 0
              to 1, where 0 is left or bottom and 1 is right or top.

            - a string using cardinal directions as abbreviation:

              - 'C' for centered
              - 'S' (south) for bottom-center
              - 'SW' (south west) for bottom-left
              - etc.

              Here is an overview of the possible positions:

              +------+------+------+
              | 'NW' | 'N'  | 'NE' |
              +------+------+------+
              | 'W'  | 'C'  | 'E'  |
              +------+------+------+
              | 'SW' | 'S'  | 'SE' |
              +------+------+------+

        share : bool, optional
            If ``True``, apply the settings to all shared Axes.
            Default is ``False``.

        See Also
        --------
        matplotlib.axes.Axes.set_aspect
            for a description of aspect handling.
        """
        if not (anchor in mtransforms.Bbox.coefs or len(anchor) == 2):
            raise ValueError('argument must be among %s' %
                             ', '.join(mtransforms.Bbox.coefs))
        if share:
            axes = set(self._shared_x_axes.get_siblings(self)
                       + self._shared_y_axes.get_siblings(self))
        else:
            axes = [self]
        for ax in axes:
            ax._anchor = anchor

        self.stale = True

    def get_data_ratio(self):
        """
        Returns the aspect ratio of the raw data.

        This method is intended to be overridden by new projection
        types.
        """
        xmin, xmax = self.get_xbound()
        ymin, ymax = self.get_ybound()

        xsize = max(abs(xmax - xmin), 1e-30)
        ysize = max(abs(ymax - ymin), 1e-30)

        return ysize / xsize

    def get_data_ratio_log(self):
        """
        Returns the aspect ratio of the raw data in log scale.
        Will be used when both axis scales are in log.
        """
        xmin, xmax = self.get_xbound()
        ymin, ymax = self.get_ybound()

        xsize = max(abs(math.log10(xmax) - math.log10(xmin)), 1e-30)
        ysize = max(abs(math.log10(ymax) - math.log10(ymin)), 1e-30)

        return ysize / xsize

    def apply_aspect(self, position=None):
        """
        Adjust the Axes for a specified data aspect ratio.

        Depending on `.get_adjustable` this will modify either the Axes box
        (position) or the view limits. In the former case, `.get_anchor`
        will affect the position.

        Notes
        -----
        This is called automatically when each Axes is drawn.  You may need
        to call it yourself if you need to update the Axes position and/or
        view limits before the Figure is drawn.

        See Also
        --------
        matplotlib.axes.Axes.set_aspect
            for a description of aspect ratio handling.
        matplotlib.axes.Axes.set_adjustable
            defining the parameter to adjust in order to meet the required
            aspect.
        matplotlib.axes.Axes.set_anchor
            defining the position in case of extra space.
        """
        if position is None:
            position = self.get_position(original=True)

        aspect = self.get_aspect()

        if self.name != 'polar':
            xscale, yscale = self.get_xscale(), self.get_yscale()
            if xscale == "linear" and yscale == "linear":
                aspect_scale_mode = "linear"
            elif xscale == "log" and yscale == "log":
                aspect_scale_mode = "log"
            elif ((xscale == "linear" and yscale == "log") or
                  (xscale == "log" and yscale == "linear")):
                if aspect != "auto":
                    warnings.warn(
                        'aspect is not supported for Axes with xscale=%s, '
                        'yscale=%s' % (xscale, yscale))
                    aspect = "auto"
            else:  # some custom projections have their own scales.
                pass
        else:
            aspect_scale_mode = "linear"

        if aspect == 'auto':
            self._set_position(position, which='active')
            return

        if aspect == 'equal':
            A = 1
        else:
            A = aspect

        figW, figH = self.get_figure().get_size_inches()
        fig_aspect = figH / figW
        if self._adjustable in ['box', 'box-forced']:
            if self in self._twinned_axes:
                raise RuntimeError("Adjustable 'box' is not allowed in a"
                                   " twinned Axes.  Use 'datalim' instead.")
            if aspect_scale_mode == "log":
                box_aspect = A * self.get_data_ratio_log()
            else:
                box_aspect = A * self.get_data_ratio()
            pb = position.frozen()
            pb1 = pb.shrunk_to_aspect(box_aspect, pb, fig_aspect)
            self._set_position(pb1.anchored(self.get_anchor(), pb), 'active')
            return

        # reset active to original in case it had been changed
        # by prior use of 'box'
        self._set_position(position, which='active')

        xmin, xmax = self.get_xbound()
        ymin, ymax = self.get_ybound()

        if aspect_scale_mode == "log":
            xmin, xmax = math.log10(xmin), math.log10(xmax)
            ymin, ymax = math.log10(ymin), math.log10(ymax)

        xsize = max(abs(xmax - xmin), 1e-30)
        ysize = max(abs(ymax - ymin), 1e-30)

        l, b, w, h = position.bounds
        box_aspect = fig_aspect * (h / w)
        data_ratio = box_aspect / A

        y_expander = (data_ratio * xsize / ysize - 1.0)
        # If y_expander > 0, the dy/dx viewLim ratio needs to increase
        if abs(y_expander) < 0.005:
            return

        if aspect_scale_mode == "log":
            dL = self.dataLim
            dL_width = math.log10(dL.x1) - math.log10(dL.x0)
            dL_height = math.log10(dL.y1) - math.log10(dL.y0)
            xr = 1.05 * dL_width
            yr = 1.05 * dL_height
        else:
            dL = self.dataLim
            xr = 1.05 * dL.width
            yr = 1.05 * dL.height

        xmarg = xsize - xr
        ymarg = ysize - yr
        Ysize = data_ratio * xsize
        Xsize = ysize / data_ratio
        Xmarg = Xsize - xr
        Ymarg = Ysize - yr
        # Setting these targets to, e.g., 0.05*xr does not seem to
        # help.
        xm = 0
        ym = 0

        shared_x = self in self._shared_x_axes
        shared_y = self in self._shared_y_axes
        # Not sure whether we need this check:
        if shared_x and shared_y:
            raise RuntimeError("adjustable='datalim' is not allowed when both"
                               " axes are shared.")

        # If y is shared, then we are only allowed to change x, etc.
        if shared_y:
            adjust_y = False
        else:
            if xmarg > xm and ymarg > ym:
                adjy = ((Ymarg > 0 and y_expander < 0) or
                        (Xmarg < 0 and y_expander > 0))
            else:
                adjy = y_expander > 0
            adjust_y = shared_x or adjy  # (Ymarg > xmarg)

        if adjust_y:
            yc = 0.5 * (ymin + ymax)
            y0 = yc - Ysize / 2.0
            y1 = yc + Ysize / 2.0
            if aspect_scale_mode == "log":
                self.set_ybound((10. ** y0, 10. ** y1))
            else:
                self.set_ybound((y0, y1))
        else:
            xc = 0.5 * (xmin + xmax)
            x0 = xc - Xsize / 2.0
            x1 = xc + Xsize / 2.0
            if aspect_scale_mode == "log":
                self.set_xbound((10. ** x0, 10. ** x1))
            else:
                self.set_xbound((x0, x1))

    def axis(self, *v, **kwargs):
        """Set axis properties.

        Valid signatures::

          xmin, xmax, ymin, ymax = axis()
          xmin, xmax, ymin, ymax = axis(list_arg)
          xmin, xmax, ymin, ymax = axis(string_arg)
          xmin, xmax, ymin, ymax = axis(**kwargs)

        Parameters
        ----------
        v : list of float or {'on', 'off', 'equal', 'tight', 'scaled',\
            'normal', 'auto', 'image', 'square'}
            Optional positional argument

            Axis data limits set from a list; or a command relating to axes:

                ========== ================================================
                Value      Description
                ========== ================================================
                'on'       Toggle axis lines and labels on
                'off'      Toggle axis lines and labels off
                'equal'    Equal scaling by changing limits
                'scaled'   Equal scaling by changing box dimensions
                'tight'    Limits set such that all data is shown
                'auto'     Automatic scaling, fill rectangle with data
                'normal'   Same as 'auto'; deprecated
                'image'    'scaled' with axis limits equal to data limits
                'square'   Square plot; similar to 'scaled', but initially\
                           forcing xmax-xmin = ymax-ymin
                ========== ================================================

        emit : bool, optional
            Passed to set_{x,y}lim functions, if observers
            are notified of axis limit change

        xmin, ymin, xmax, ymax : float, optional
            The axis limits to be set

        Returns
        -------
        xmin, xmax, ymin, ymax : float
            The axis limits

        """

        if len(v) == 0 and len(kwargs) == 0:
            xmin, xmax = self.get_xlim()
            ymin, ymax = self.get_ylim()
            return xmin, xmax, ymin, ymax

        emit = kwargs.get('emit', True)

        if len(v) == 1 and isinstance(v[0], six.string_types):
            s = v[0].lower()
            if s == 'on':
                self.set_axis_on()
            elif s == 'off':
                self.set_axis_off()
            elif s in ('equal', 'tight', 'scaled', 'normal',
                       'auto', 'image', 'square'):
                self.set_autoscale_on(True)
                self.set_aspect('auto')
                self.autoscale_view(tight=False)
                # self.apply_aspect()
                if s == 'equal':
                    self.set_aspect('equal', adjustable='datalim')
                elif s == 'scaled':
                    self.set_aspect('equal', adjustable='box', anchor='C')
                    self.set_autoscale_on(False)  # Req. by Mark Bakker
                elif s == 'tight':
                    self.autoscale_view(tight=True)
                    self.set_autoscale_on(False)
                elif s == 'image':
                    self.autoscale_view(tight=True)
                    self.set_autoscale_on(False)
                    self.set_aspect('equal', adjustable='box', anchor='C')
                elif s == 'square':
                    self.set_aspect('equal', adjustable='box', anchor='C')
                    self.set_autoscale_on(False)
                    xlim = self.get_xlim()
                    ylim = self.get_ylim()
                    edge_size = max(np.diff(xlim), np.diff(ylim))
                    self.set_xlim([xlim[0], xlim[0] + edge_size],
                                  emit=emit, auto=False)
                    self.set_ylim([ylim[0], ylim[0] + edge_size],
                                  emit=emit, auto=False)
            else:
                raise ValueError('Unrecognized string %s to axis; '
                                 'try on or off' % s)
            xmin, xmax = self.get_xlim()
            ymin, ymax = self.get_ylim()
            return xmin, xmax, ymin, ymax

        try:
            v[0]
        except IndexError:
            xmin = kwargs.get('xmin', None)
            xmax = kwargs.get('xmax', None)
            auto = False  # turn off autoscaling, unless...
            if xmin is None and xmax is None:
                auto = None  # leave autoscaling state alone
            xmin, xmax = self.set_xlim(xmin, xmax, emit=emit, auto=auto)

            ymin = kwargs.get('ymin', None)
            ymax = kwargs.get('ymax', None)
            auto = False  # turn off autoscaling, unless...
            if ymin is None and ymax is None:
                auto = None  # leave autoscaling state alone
            ymin, ymax = self.set_ylim(ymin, ymax, emit=emit, auto=auto)
            return xmin, xmax, ymin, ymax

        v = v[0]
        if len(v) != 4:
            raise ValueError('v must contain [xmin xmax ymin ymax]')

        self.set_xlim([v[0], v[1]], emit=emit, auto=False)
        self.set_ylim([v[2], v[3]], emit=emit, auto=False)

        return v

    def get_legend(self):
        """Return the `Legend` instance, or None if no legend is defined."""
        return self.legend_

    def get_images(self):
        """return a list of Axes images contained by the Axes"""
        return cbook.silent_list('AxesImage', self.images)

    def get_lines(self):
        """Return a list of lines contained by the Axes"""
        return cbook.silent_list('Line2D', self.lines)

    def get_xaxis(self):
        """Return the XAxis instance."""
        return self.xaxis

    def get_xgridlines(self):
        """Get the x grid lines as a list of `Line2D` instances."""
        return cbook.silent_list('Line2D xgridline',
                                 self.xaxis.get_gridlines())

    def get_xticklines(self):
        """Get the x tick lines as a list of `Line2D` instances."""
        return cbook.silent_list('Line2D xtickline',
                                 self.xaxis.get_ticklines())

    def get_yaxis(self):
        """Return the YAxis instance."""
        return self.yaxis

    def get_ygridlines(self):
        """Get the y grid lines as a list of `Line2D` instances."""
        return cbook.silent_list('Line2D ygridline',
                                 self.yaxis.get_gridlines())

    def get_yticklines(self):
        """Get the y tick lines as a list of `Line2D` instances."""
        return cbook.silent_list('Line2D ytickline',
                                 self.yaxis.get_ticklines())

    # Adding and tracking artists

    def _sci(self, im):
        """
        helper for :func:`~matplotlib.pyplot.sci`;
        do not use elsewhere.
        """
        if isinstance(im, matplotlib.contour.ContourSet):
            if im.collections[0] not in self.collections:
                raise ValueError(
                    "ContourSet must be in current Axes")
        elif im not in self.images and im not in self.collections:
            raise ValueError(
                "Argument must be an image, collection, or ContourSet in "
                "this Axes")
        self._current_image = im

    def _gci(self):
        """
        Helper for :func:`~matplotlib.pyplot.gci`;
        do not use elsewhere.
        """
        return self._current_image

    def has_data(self):
        """
        Return *True* if any artists have been added to axes.

        This should not be used to determine whether the *dataLim*
        need to be updated, and may not actually be useful for
        anything.
        """
        return (
            len(self.collections) +
            len(self.images) +
            len(self.lines) +
            len(self.patches)) > 0

    def add_artist(self, a):
        """Add any :class:`~matplotlib.artist.Artist` to the axes.

        Use `add_artist` only for artists for which there is no dedicated
        "add" method; and if necessary, use a method such as `update_datalim`
        to manually update the dataLim if the artist is to be included in
        autoscaling.

        Returns the artist.
        """
        a.axes = self
        self.artists.append(a)
        self._set_artist_props(a)
        a.set_clip_path(self.patch)
        a._remove_method = lambda h: self.artists.remove(h)
        self.stale = True
        return a

    def add_collection(self, collection, autolim=True):
        """
        Add a :class:`~matplotlib.collections.Collection` instance
        to the axes.

        Returns the collection.
        """
        label = collection.get_label()
        if not label:
            collection.set_label('_collection%d' % len(self.collections))
        self.collections.append(collection)
        self._set_artist_props(collection)

        if collection.get_clip_path() is None:
            collection.set_clip_path(self.patch)

        if autolim:
            self.update_datalim(collection.get_datalim(self.transData))

        collection._remove_method = lambda h: self.collections.remove(h)
        self.stale = True
        return collection

    def add_image(self, image):
        """
        Add a :class:`~matplotlib.image.AxesImage` to the axes.

        Returns the image.
        """
        self._set_artist_props(image)
        if not image.get_label():
            image.set_label('_image%d' % len(self.images))
        self.images.append(image)
        image._remove_method = lambda h: self.images.remove(h)
        self.stale = True
        return image

    def _update_image_limits(self, image):
        xmin, xmax, ymin, ymax = image.get_extent()
        self.axes.update_datalim(((xmin, ymin), (xmax, ymax)))

    def add_line(self, line):
        """
        Add a :class:`~matplotlib.lines.Line2D` to the list of plot
        lines

        Returns the line.
        """
        self._set_artist_props(line)
        if line.get_clip_path() is None:
            line.set_clip_path(self.patch)

        self._update_line_limits(line)
        if not line.get_label():
            line.set_label('_line%d' % len(self.lines))
        self.lines.append(line)
        line._remove_method = lambda h: self.lines.remove(h)
        self.stale = True
        return line

    def _add_text(self, txt):
        """

        """
        self._set_artist_props(txt)
        self.texts.append(txt)
        txt._remove_method = lambda h: self.texts.remove(h)
        self.stale = True
        return txt

    def _update_line_limits(self, line):
        """
        Figures out the data limit of the given line, updating self.dataLim.
        """
        path = line.get_path()
        if path.vertices.size == 0:
            return

        line_trans = line.get_transform()

        if line_trans == self.transData:
            data_path = path

        elif any(line_trans.contains_branch_seperately(self.transData)):
            # identify the transform to go from line's coordinates
            # to data coordinates
            trans_to_data = line_trans - self.transData

            # if transData is affine we can use the cached non-affine component
            # of line's path. (since the non-affine part of line_trans is
            # entirely encapsulated in trans_to_data).
            if self.transData.is_affine:
                line_trans_path = line._get_transformed_path()
                na_path, _ = line_trans_path.get_transformed_path_and_affine()
                data_path = trans_to_data.transform_path_affine(na_path)
            else:
                data_path = trans_to_data.transform_path(path)
        else:
            # for backwards compatibility we update the dataLim with the
            # coordinate range of the given path, even though the coordinate
            # systems are completely different. This may occur in situations
            # such as when ax.transAxes is passed through for absolute
            # positioning.
            data_path = path

        if data_path.vertices.size > 0:
            updatex, updatey = line_trans.contains_branch_seperately(
                self.transData)
            self.dataLim.update_from_path(data_path,
                                          self.ignore_existing_data_limits,
                                          updatex=updatex,
                                          updatey=updatey)
            self.ignore_existing_data_limits = False

    def add_patch(self, p):
        """
        Add a :class:`~matplotlib.patches.Patch` *p* to the list of
        axes patches; the clipbox will be set to the Axes clipping
        box.  If the transform is not set, it will be set to
        :attr:`transData`.

        Returns the patch.
        """

        self._set_artist_props(p)
        if p.get_clip_path() is None:
            p.set_clip_path(self.patch)
        self._update_patch_limits(p)
        self.patches.append(p)
        p._remove_method = lambda h: self.patches.remove(h)
        return p

    def _update_patch_limits(self, patch):
        """update the data limits for patch *p*"""
        # hist can add zero height Rectangles, which is useful to keep
        # the bins, counts and patches lined up, but it throws off log
        # scaling.  We'll ignore rects with zero height or width in
        # the auto-scaling

        # cannot check for '==0' since unitized data may not compare to zero
        # issue #2150 - we update the limits if patch has non zero width
        # or height.
        if (isinstance(patch, mpatches.Rectangle) and
                ((not patch.get_width()) and (not patch.get_height()))):
            return
        vertices = patch.get_path().vertices
        if vertices.size > 0:
            xys = patch.get_patch_transform().transform(vertices)
            if patch.get_data_transform() != self.transData:
                patch_to_data = (patch.get_data_transform() -
                                 self.transData)
                xys = patch_to_data.transform(xys)

            updatex, updatey = patch.get_transform().\
                contains_branch_seperately(self.transData)
            self.update_datalim(xys, updatex=updatex,
                                updatey=updatey)

    def add_table(self, tab):
        """
        Add a :class:`~matplotlib.table.Table` instance to the
        list of axes tables

        Parameters
        ----------
        tab: `matplotlib.table.Table`
            Table instance

        Returns
        -------
        `matplotlib.table.Table`: the table.
        """
        self._set_artist_props(tab)
        self.tables.append(tab)
        tab.set_clip_path(self.patch)
        tab._remove_method = lambda h: self.tables.remove(h)
        return tab

    def add_container(self, container):
        """
        Add a :class:`~matplotlib.container.Container` instance
        to the axes.

        Returns the collection.
        """
        label = container.get_label()
        if not label:
            container.set_label('_container%d' % len(self.containers))
        self.containers.append(container)
        container.set_remove_method(lambda h: self.containers.remove(h))
        return container

    def _on_units_changed(self, scalex=False, scaley=False):
        """
        Callback for processing changes to axis units.

        Currently forces updates of data limits and view limits.
        """
        self.relim()
        self.autoscale_view(scalex=scalex, scaley=scaley)

    def relim(self, visible_only=False):
        """
        Recompute the data limits based on current artists. If you want to
        exclude invisible artists from the calculation, set
        ``visible_only=True``

        At present, :class:`~matplotlib.collections.Collection`
        instances are not supported.
        """
        # Collections are deliberately not supported (yet); see
        # the TODO note in artists.py.
        self.dataLim.ignore(True)
        self.dataLim.set_points(mtransforms.Bbox.null().get_points())
        self.ignore_existing_data_limits = True

        for line in self.lines:
            if not visible_only or line.get_visible():
                self._update_line_limits(line)

        for p in self.patches:
            if not visible_only or p.get_visible():
                self._update_patch_limits(p)

        for image in self.images:
            if not visible_only or image.get_visible():
                self._update_image_limits(image)

    def update_datalim(self, xys, updatex=True, updatey=True):
        """
        Update the data lim bbox with seq of xy tups or equiv. 2-D array
        """
        # if no data is set currently, the bbox will ignore its
        # limits and set the bound to be the bounds of the xydata.
        # Otherwise, it will compute the bounds of it's current data
        # and the data in xydata
        xys = np.asarray(xys)
        if not len(xys):
            return
        self.dataLim.update_from_data_xy(xys, self.ignore_existing_data_limits,
                                         updatex=updatex, updatey=updatey)
        self.ignore_existing_data_limits = False

    def update_datalim_bounds(self, bounds):
        """
        Update the datalim to include the given
        :class:`~matplotlib.transforms.Bbox` *bounds*
        """
        self.dataLim.set(mtransforms.Bbox.union([self.dataLim, bounds]))

    def _process_unit_info(self, xdata=None, ydata=None, kwargs=None):
        """Look for unit *kwargs* and update the axis instances as necessary"""

        if self.xaxis is None or self.yaxis is None:
            return

        if xdata is not None:
            # we only need to update if there is nothing set yet.
            if not self.xaxis.have_units():
                self.xaxis.update_units(xdata)

        if ydata is not None:
            # we only need to update if there is nothing set yet.
            if not self.yaxis.have_units():
                self.yaxis.update_units(ydata)

        # process kwargs 2nd since these will override default units
        if kwargs is not None:
            xunits = kwargs.pop('xunits', self.xaxis.units)
            if self.name == 'polar':
                xunits = kwargs.pop('thetaunits', xunits)
            if xunits != self.xaxis.units:
                self.xaxis.set_units(xunits)
                # If the units being set imply a different converter,
                # we need to update.
                if xdata is not None:
                    self.xaxis.update_units(xdata)

            yunits = kwargs.pop('yunits', self.yaxis.units)
            if self.name == 'polar':
                yunits = kwargs.pop('runits', yunits)
            if yunits != self.yaxis.units:
                self.yaxis.set_units(yunits)
                # If the units being set imply a different converter,
                # we need to update.
                if ydata is not None:
                    self.yaxis.update_units(ydata)
        return kwargs

    def in_axes(self, mouseevent):
        """
        Return *True* if the given *mouseevent* (in display coords)
        is in the Axes
        """
        return self.patch.contains(mouseevent)[0]

    def get_autoscale_on(self):
        """
        Get whether autoscaling is applied for both axes on plot commands
        """
        return self._autoscaleXon and self._autoscaleYon

    def get_autoscalex_on(self):
        """
        Get whether autoscaling for the x-axis is applied on plot commands
        """
        return self._autoscaleXon

    def get_autoscaley_on(self):
        """
        Get whether autoscaling for the y-axis is applied on plot commands
        """
        return self._autoscaleYon

    def set_autoscale_on(self, b):
        """
        Set whether autoscaling is applied on plot commands

        .. ACCEPTS: bool

        Parameters
        ----------
        b : bool
        """
        self._autoscaleXon = b
        self._autoscaleYon = b

    def set_autoscalex_on(self, b):
        """
        Set whether autoscaling for the x-axis is applied on plot commands

        .. ACCEPTS: bool

        Parameters
        ----------
        b : bool
        """
        self._autoscaleXon = b

    def set_autoscaley_on(self, b):
        """
        Set whether autoscaling for the y-axis is applied on plot commands

        .. ACCEPTS: bool

        Parameters
        ----------
        b : bool
        """
        self._autoscaleYon = b

    @property
    def use_sticky_edges(self):
        """
        When autoscaling, whether to obey all `Artist.sticky_edges`.

        Default is ``True``.

        Setting this to ``False`` ensures that the specified margins
        will be applied, even if the plot includes an image, for
        example, which would otherwise force a view limit to coincide
        with its data limit.

        The changing this property does not change the plot until
        `autoscale` or `autoscale_view` is called.
        """
        return self._use_sticky_edges

    @use_sticky_edges.setter
    def use_sticky_edges(self, b):
        self._use_sticky_edges = bool(b)
        # No effect until next autoscaling, which will mark the axes as stale.

    def set_xmargin(self, m):
        """
        Set padding of X data limits prior to autoscaling.

        *m* times the data interval will be added to each
        end of that interval before it is used in autoscaling.
        For example, if your data is in the range [0, 2], a factor of
        ``m = 0.1`` will result in a range [-0.2, 2.2].

        Negative values -0.5 < m < 0 will result in clipping of the data range.
        I.e. for a data range [0, 2], a factor of ``m = -0.1`` will result in
        a range [0.2, 1.8].

        .. ACCEPTS: float greater than -0.5

        Parameters
        ----------
        m : float greater than -0.5
        """
        if m <= -0.5:
            raise ValueError("margin must be greater than -0.5")
        self._xmargin = m
        self.stale = True

    def set_ymargin(self, m):
        """
        Set padding of Y data limits prior to autoscaling.

        *m* times the data interval will be added to each
        end of that interval before it is used in autoscaling.
        For example, if your data is in the range [0, 2], a factor of
        ``m = 0.1`` will result in a range [-0.2, 2.2].

        Negative values -0.5 < m < 0 will result in clipping of the data range.
        I.e. for a data range [0, 2], a factor of ``m = -0.1`` will result in
        a range [0.2, 1.8].

        .. ACCEPTS: float greater than -0.5

        Parameters
        ----------
        m : float greater than -0.5
        """
        if m <= -0.5:
            raise ValueError("margin must be greater than -0.5")
        self._ymargin = m
        self.stale = True

    def margins(self, *args, **kw):
        """
        Set or retrieve autoscaling margins.

        signatures::

            margins()

        returns xmargin, ymargin

        ::

            margins(margin)

            margins(xmargin, ymargin)

            margins(x=xmargin, y=ymargin)

            margins(..., tight=False)

        All three forms above set the xmargin and ymargin parameters.
        All keyword parameters are optional.  A single argument
        specifies both xmargin and ymargin. The padding added to the end of
        each interval is *margin* times the data interval. The *margin* must
        be a float in the range [0, 1].

        The *tight* parameter is passed to :meth:`autoscale_view`
        , which is executed after a margin is changed; the default here is
        *True*, on the assumption that when margins are specified, no
        additional padding to match tick marks is usually desired.  Setting
        *tight* to *None* will preserve the previous setting.

        Specifying any margin changes only the autoscaling; for example,
        if *xmargin* is not None, then *xmargin* times the X data
        interval will be added to each end of that interval before
        it is used in autoscaling.

        """
        if not args and not kw:
            return self._xmargin, self._ymargin

        tight = kw.pop('tight', True)
        mx = kw.pop('x', None)
        my = kw.pop('y', None)
        if len(args) == 1:
            mx = my = args[0]
        elif len(args) == 2:
            mx, my = args
        elif len(args) > 2:
            raise ValueError("more than two arguments were supplied")
        if mx is not None:
            self.set_xmargin(mx)
        if my is not None:
            self.set_ymargin(my)

        scalex = (mx is not None)
        scaley = (my is not None)

        self.autoscale_view(tight=tight, scalex=scalex, scaley=scaley)

    def set_rasterization_zorder(self, z):
        """
        Parameters
        ----------
        z : float or None
            zorder below which artists are rasterized.  ``None`` means that
            artists do not get rasterized based on zorder.

            .. ACCEPTS: float or None
        """
        self._rasterization_zorder = z
        self.stale = True

    def get_rasterization_zorder(self):
        """Return the zorder value below which artists will be rasterized."""
        return self._rasterization_zorder

    def autoscale(self, enable=True, axis='both', tight=None):
        """
        Autoscale the axis view to the data (toggle).

        Convenience method for simple axis view autoscaling.
        It turns autoscaling on or off, and then,
        if autoscaling for either axis is on, it performs
        the autoscaling on the specified axis or axes.

        Parameters
        ----------
        enable : bool or None, optional
            True (default) turns autoscaling on, False turns it off.
            None leaves the autoscaling state unchanged.

        axis : ['both' | 'x' | 'y'], optional
            which axis to operate on; default is 'both'

        tight: bool or None, optional
            If True, set view limits to data limits;
            if False, let the locator and margins expand the view limits;
            if None, use tight scaling if the only artist is an image,
            otherwise treat *tight* as False.
            The *tight* setting is retained for future autoscaling
            until it is explicitly changed.

        """
        if enable is None:
            scalex = True
            scaley = True
        else:
            scalex = False
            scaley = False
            if axis in ['x', 'both']:
                self._autoscaleXon = bool(enable)
                scalex = self._autoscaleXon
            if axis in ['y', 'both']:
                self._autoscaleYon = bool(enable)
                scaley = self._autoscaleYon
        if tight and scalex:
            self._xmargin = 0
        if tight and scaley:
            self._ymargin = 0
        self.autoscale_view(tight=tight, scalex=scalex, scaley=scaley)

    def autoscale_view(self, tight=None, scalex=True, scaley=True):
        """
        Autoscale the view limits using the data limits.

        You can selectively autoscale only a single axis, e.g., the xaxis by
        setting *scaley* to *False*.  The autoscaling preserves any
        axis direction reversal that has already been done.

        If *tight* is *False*, the axis major locator will be used
        to expand the view limits if rcParams['axes.autolimit_mode']
        is 'round_numbers'.  Note that any margins that are in effect
        will be applied first, regardless of whether *tight* is
        *True* or *False*.  Specifying *tight* as *True* or *False*
        saves the setting as a private attribute of the Axes; specifying
        it as *None* (the default) applies the previously saved value.

        The data limits are not updated automatically when artist data are
        changed after the artist has been added to an Axes instance.  In that
        case, use :meth:`matplotlib.axes.Axes.relim` prior to calling
        autoscale_view.
        """
        if tight is not None:
            self._tight = bool(tight)

        if self.use_sticky_edges and (self._xmargin or self._ymargin):
            stickies = [artist.sticky_edges for artist in self.get_children()]
            x_stickies = sum([sticky.x for sticky in stickies], [])
            y_stickies = sum([sticky.y for sticky in stickies], [])
            if self.get_xscale().lower() == 'log':
                x_stickies = [xs for xs in x_stickies if xs > 0]
            if self.get_yscale().lower() == 'log':
                y_stickies = [ys for ys in y_stickies if ys > 0]
        else:  # Small optimization.
            x_stickies, y_stickies = [], []

        def handle_single_axis(scale, autoscaleon, shared_axes, interval,
                               minpos, axis, margin, stickies, set_bound):

            if not (scale and autoscaleon):
                return  # nothing to do...

            shared = shared_axes.get_siblings(self)
            dl = [ax.dataLim for ax in shared]
            # ignore non-finite data limits if good limits exist
            finite_dl = [d for d in dl if np.isfinite(d).all()]
            if len(finite_dl):
                # if finite limits exist for atleast one axis (and the
                # other is infinite), restore the finite limits
                x_finite = [d for d in dl
                            if (np.isfinite(d.intervalx).all() and
                                (d not in finite_dl))]
                y_finite = [d for d in dl
                            if (np.isfinite(d.intervaly).all() and
                                (d not in finite_dl))]

                dl = finite_dl
                dl.extend(x_finite)
                dl.extend(y_finite)

            bb = mtransforms.BboxBase.union(dl)
            x0, x1 = getattr(bb, interval)
            locator = axis.get_major_locator()
            try:
                # e.g., DateLocator has its own nonsingular()
                x0, x1 = locator.nonsingular(x0, x1)
            except AttributeError:
                # Default nonsingular for, e.g., MaxNLocator
                x0, x1 = mtransforms.nonsingular(
                    x0, x1, increasing=False, expander=0.05)

            # Add the margin in figure space and then transform back, to handle
            # non-linear scales.
            minpos = getattr(bb, minpos)
            transform = axis.get_transform()
            inverse_trans = transform.inverted()
            # We cannot use exact equality due to floating point issues e.g.
            # with streamplot.
            do_lower_margin = not np.any(np.isclose(x0, stickies))
            do_upper_margin = not np.any(np.isclose(x1, stickies))
            x0, x1 = axis._scale.limit_range_for_scale(x0, x1, minpos)
            x0t, x1t = transform.transform([x0, x1])
            delta = (x1t - x0t) * margin
            if do_lower_margin:
                x0t -= delta
            if do_upper_margin:
                x1t += delta
            x0, x1 = inverse_trans.transform([x0t, x1t])

            if not self._tight:
                x0, x1 = locator.view_limits(x0, x1)
            set_bound(x0, x1)
            # End of definition of internal function 'handle_single_axis'.

        handle_single_axis(
            scalex, self._autoscaleXon, self._shared_x_axes, 'intervalx',
            'minposx', self.xaxis, self._xmargin, x_stickies, self.set_xbound)
        handle_single_axis(
            scaley, self._autoscaleYon, self._shared_y_axes, 'intervaly',
            'minposy', self.yaxis, self._ymargin, y_stickies, self.set_ybound)

    def _get_axis_list(self):
        return (self.xaxis, self.yaxis)

    # Drawing

    @allow_rasterization
    def draw(self, renderer=None, inframe=False):
        """Draw everything (plot lines, axes, labels)"""
        if renderer is None:
            renderer = self._cachedRenderer

        if renderer is None:
            raise RuntimeError('No renderer defined')
        if not self.get_visible():
            return
        renderer.open_group('axes')
        # prevent triggering call backs during the draw process
        self._stale = True
        locator = self.get_axes_locator()
        if locator:
            pos = locator(self, renderer)
            self.apply_aspect(pos)
        else:
            self.apply_aspect()

        artists = self.get_children()
        artists.remove(self.patch)

        # the frame draws the edges around the axes patch -- we
        # decouple these so the patch can be in the background and the
        # frame in the foreground. Do this before drawing the axis
        # objects so that the spine has the opportunity to update them.
        if not (self.axison and self._frameon):
            for spine in six.itervalues(self.spines):
                artists.remove(spine)

        if self.axison and not inframe:
            if self._axisbelow is True:
                self.xaxis.set_zorder(0.5)
                self.yaxis.set_zorder(0.5)
            elif self._axisbelow is False:
                self.xaxis.set_zorder(2.5)
                self.yaxis.set_zorder(2.5)
            else:
                # 'line': above patches, below lines
                self.xaxis.set_zorder(1.5)
                self.yaxis.set_zorder(1.5)
        else:
            for _axis in self._get_axis_list():
                artists.remove(_axis)

        if inframe:
            artists.remove(self.title)
            artists.remove(self._left_title)
            artists.remove(self._right_title)

        if not self.figure.canvas.is_saving():
            artists = [a for a in artists
                       if not a.get_animated() or a in self.images]
        artists = sorted(artists, key=attrgetter('zorder'))

        # rasterize artists with negative zorder
        # if the minimum zorder is negative, start rasterization
        rasterization_zorder = self._rasterization_zorder
        if (rasterization_zorder is not None and
                artists and artists[0].zorder < rasterization_zorder):
            renderer.start_rasterizing()
            artists_rasterized = [a for a in artists
                                  if a.zorder < rasterization_zorder]
            artists = [a for a in artists
                       if a.zorder >= rasterization_zorder]
        else:
            artists_rasterized = []

        # the patch draws the background rectangle -- the frame below
        # will draw the edges
        if self.axison and self._frameon:
            self.patch.draw(renderer)

        if artists_rasterized:
            for a in artists_rasterized:
                a.draw(renderer)
            renderer.stop_rasterizing()

        mimage._draw_list_compositing_images(renderer, self, artists)

        renderer.close_group('axes')
        self._cachedRenderer = renderer
        self.stale = False

    def draw_artist(self, a):
        """
        This method can only be used after an initial draw which
        caches the renderer.  It is used to efficiently update Axes
        data (axis ticks, labels, etc are not updated)
        """
        if self._cachedRenderer is None:
            raise AttributeError("draw_artist can only be used after an "
                                 "initial draw which caches the renderer")
        a.draw(self._cachedRenderer)

    def redraw_in_frame(self):
        """
        This method can only be used after an initial draw which
        caches the renderer.  It is used to efficiently update Axes
        data (axis ticks, labels, etc are not updated)
        """
        if self._cachedRenderer is None:
            raise AttributeError("redraw_in_frame can only be used after an "
                                 "initial draw which caches the renderer")
        self.draw(self._cachedRenderer, inframe=True)

    def get_renderer_cache(self):
        return self._cachedRenderer

    # Axes rectangle characteristics

    def get_frame_on(self):
        """
        Get whether the axes rectangle patch is drawn.
        """
        return self._frameon

    def set_frame_on(self, b):
        """
        Set whether the axes rectangle patch is drawn.

        .. ACCEPTS: bool

        Parameters
        ----------
        b : bool
        """
        self._frameon = b
        self.stale = True

    def get_axisbelow(self):
        """
        Get whether axis ticks and gridlines are above or below most artists.
        """
        return self._axisbelow

    def set_axisbelow(self, b):
        """
        Set whether axis ticks and gridlines are above or below most artists.

        .. ACCEPTS: [ bool | 'line' ]

        Parameters
        ----------
        b : bool or 'line'
        """
        self._axisbelow = validate_axisbelow(b)
        self.stale = True

    @docstring.dedent_interpd
    def grid(self, b=None, which='major', axis='both', **kwargs):
        """
        Turn the axes grids on or off.

        Set the axes grids on or off; *b* is a boolean.

        If *b* is *None* and ``len(kwargs)==0``, toggle the grid state.  If
        *kwargs* are supplied, it is assumed that you want a grid and *b*
        is thus set to *True*.

        *which* can be 'major' (default), 'minor', or 'both' to control
        whether major tick grids, minor tick grids, or both are affected.

        *axis* can be 'both' (default), 'x', or 'y' to control which
        set of gridlines are drawn.

        *kwargs* are used to set the grid line properties, e.g.,::

           ax.grid(color='r', linestyle='-', linewidth=2)

        Valid :class:`~matplotlib.lines.Line2D` kwargs are

        %(Line2D)s

        """
        if len(kwargs):
            b = True
        elif b is not None:
            b = _string_to_bool(b)

        if axis == 'x' or axis == 'both':
            self.xaxis.grid(b, which=which, **kwargs)
        if axis == 'y' or axis == 'both':
            self.yaxis.grid(b, which=which, **kwargs)

    def ticklabel_format(self, **kwargs):
        """
        Change the `~matplotlib.ticker.ScalarFormatter` used by
        default for linear axes.

        Optional keyword arguments:

          ==============   =========================================
          Keyword          Description
          ==============   =========================================
          *style*          [ 'sci' (or 'scientific') | 'plain' ]
                           plain turns off scientific notation
          *scilimits*      (m, n), pair of integers; if *style*
                           is 'sci', scientific notation will
                           be used for numbers outside the range
                           10`m`:sup: to 10`n`:sup:.
                           Use (0,0) to include all numbers.
          *useOffset*      [ bool | offset ]; if True,
                           the offset will be calculated as needed;
                           if False, no offset will be used; if a
                           numeric offset is specified, it will be
                           used.
          *axis*           [ 'x' | 'y' | 'both' ]
          *useLocale*      If True, format the number according to
                           the current locale.  This affects things
                           such as the character used for the
                           decimal separator.  If False, use
                           C-style (English) formatting.  The
                           default setting is controlled by the
                           axes.formatter.use_locale rcparam.
          *useMathText*    If True, render the offset and scientific
                           notation in mathtext
          ==============   =========================================

        Only the major ticks are affected.
        If the method is called when the
        :class:`~matplotlib.ticker.ScalarFormatter` is not the
        :class:`~matplotlib.ticker.Formatter` being used, an
        :exc:`AttributeError` will be raised.

        """
        style = kwargs.pop('style', '').lower()
        scilimits = kwargs.pop('scilimits', None)
        useOffset = kwargs.pop('useOffset', None)
        useLocale = kwargs.pop('useLocale', None)
        useMathText = kwargs.pop('useMathText', None)
        axis = kwargs.pop('axis', 'both').lower()
        if scilimits is not None:
            try:
                m, n = scilimits
                m + n + 1  # check that both are numbers
            except (ValueError, TypeError):
                raise ValueError("scilimits must be a sequence of 2 integers")
        if style[:3] == 'sci':
            sb = True
        elif style == 'plain':
            sb = False
        elif style == 'comma':
            raise NotImplementedError("comma style remains to be added")
        elif style == '':
            sb = None
        else:
            raise ValueError("%s is not a valid style value")
        try:
            if sb is not None:
                if axis == 'both' or axis == 'x':
                    self.xaxis.major.formatter.set_scientific(sb)
                if axis == 'both' or axis == 'y':
                    self.yaxis.major.formatter.set_scientific(sb)
            if scilimits is not None:
                if axis == 'both' or axis == 'x':
                    self.xaxis.major.formatter.set_powerlimits(scilimits)
                if axis == 'both' or axis == 'y':
                    self.yaxis.major.formatter.set_powerlimits(scilimits)
            if useOffset is not None:
                if axis == 'both' or axis == 'x':
                    self.xaxis.major.formatter.set_useOffset(useOffset)
                if axis == 'both' or axis == 'y':
                    self.yaxis.major.formatter.set_useOffset(useOffset)
            if useLocale is not None:
                if axis == 'both' or axis == 'x':
                    self.xaxis.major.formatter.set_useLocale(useLocale)
                if axis == 'both' or axis == 'y':
                    self.yaxis.major.formatter.set_useLocale(useLocale)
            if useMathText is not None:
                if axis == 'both' or axis == 'x':
                    self.xaxis.major.formatter.set_useMathText(useMathText)
                if axis == 'both' or axis == 'y':
                    self.yaxis.major.formatter.set_useMathText(useMathText)
        except AttributeError:
            raise AttributeError(
                "This method only works with the ScalarFormatter.")

    def locator_params(self, axis='both', tight=None, **kwargs):
        """
        Control behavior of tick locators.

        Parameters
        ----------
        axis : ['both' | 'x' | 'y'], optional
            The axis on which to operate.

        tight : bool or None, optional
            Parameter passed to :meth:`autoscale_view`.
            Default is None, for no change.

        Other Parameters
        ----------------
        **kw :
            Remaining keyword arguments are passed to directly to the
            :meth:`~matplotlib.ticker.MaxNLocator.set_params` method.

        Typically one might want to reduce the maximum number
        of ticks and use tight bounds when plotting small
        subplots, for example::

            ax.locator_params(tight=True, nbins=4)

        Because the locator is involved in autoscaling,
        :meth:`autoscale_view` is called automatically after
        the parameters are changed.

        This presently works only for the
        :class:`~matplotlib.ticker.MaxNLocator` used
        by default on linear axes, but it may be generalized.
        """
        _x = axis in ['x', 'both']
        _y = axis in ['y', 'both']
        if _x:
            self.xaxis.get_major_locator().set_params(**kwargs)
        if _y:
            self.yaxis.get_major_locator().set_params(**kwargs)
        self.autoscale_view(tight=tight, scalex=_x, scaley=_y)

    def tick_params(self, axis='both', **kwargs):
        """Change the appearance of ticks, tick labels, and gridlines.

        Parameters
        ----------
        axis : {'x', 'y', 'both'}, optional
            Which axis to apply the parameters to.

        Other Parameters
        ----------------

        axis : {'x', 'y', 'both'}
            Axis on which to operate; default is 'both'.

        reset : bool
            If *True*, set all parameters to defaults
            before processing other keyword arguments.  Default is
            *False*.

        which : {'major', 'minor', 'both'}
            Default is 'major'; apply arguments to *which* ticks.

        direction : {'in', 'out', 'inout'}
            Puts ticks inside the axes, outside the axes, or both.

        length : float
            Tick length in points.

        width : float
            Tick width in points.

        color : color
            Tick color; accepts any mpl color spec.

        pad : float
            Distance in points between tick and label.

        labelsize : float or str
            Tick label font size in points or as a string (e.g., 'large').

        labelcolor : color
            Tick label color; mpl color spec.

        colors : color
            Changes the tick color and the label color to the same value:
            mpl color spec.

        zorder : float
            Tick and label zorder.

        bottom, top, left, right : bool
            Whether to draw the respective ticks.

        labelbottom, labeltop, labelleft, labelright : bool
            Whether to draw the respective tick labels.

        labelrotation : float
            Tick label rotation

        grid_color : color
            Changes the gridline color to the given mpl color spec.

        grid_alpha : float
            Transparency of gridlines: 0 (transparent) to 1 (opaque).

        grid_linewidth : float
            Width of gridlines in points.

        grid_linestyle : string
            Any valid :class:`~matplotlib.lines.Line2D` line style spec.

        Examples
        --------

        Usage ::

            ax.tick_params(direction='out', length=6, width=2, colors='r',
                           grid_color='r', grid_alpha=0.5)

        This will make all major ticks be red, pointing out of the box,
        and with dimensions 6 points by 2 points.  Tick labels will
        also be red.  Gridlines will be red and translucent.

        """
        if axis in ['x', 'both']:
            xkw = dict(kwargs)
            xkw.pop('left', None)
            xkw.pop('right', None)
            xkw.pop('labelleft', None)
            xkw.pop('labelright', None)
            self.xaxis.set_tick_params(**xkw)
        if axis in ['y', 'both']:
            ykw = dict(kwargs)
            ykw.pop('top', None)
            ykw.pop('bottom', None)
            ykw.pop('labeltop', None)
            ykw.pop('labelbottom', None)
            self.yaxis.set_tick_params(**ykw)

    def set_axis_off(self):
        """Turn off the axis."""
        self.axison = False
        self.stale = True

    def set_axis_on(self):
        """Turn on the axis."""
        self.axison = True
        self.stale = True

    # data limits, ticks, tick labels, and formatting

    def invert_xaxis(self):
        """Invert the x-axis."""
        self.set_xlim(self.get_xlim()[::-1], auto=None)

    def xaxis_inverted(self):
        """Return whether the x-axis is inverted."""
        left, right = self.get_xlim()
        return right < left

    def get_xbound(self):
        """Return the lower and upper x-axis bounds, in increasing order."""
        left, right = self.get_xlim()
        if left < right:
            return left, right
        else:
            return right, left

    def set_xbound(self, lower=None, upper=None):
        """
        Set the lower and upper numerical bounds of the x-axis.

        This method will honor axes inversion regardless of parameter order.
        It will not change the _autoscaleXon attribute.

        .. ACCEPTS: (lower: float, upper: float)
        """
        if upper is None and iterable(lower):
            lower, upper = lower

        old_lower, old_upper = self.get_xbound()

        if lower is None:
            lower = old_lower
        if upper is None:
            upper = old_upper

        if self.xaxis_inverted():
            if lower < upper:
                self.set_xlim(upper, lower, auto=None)
            else:
                self.set_xlim(lower, upper, auto=None)
        else:
            if lower < upper:
                self.set_xlim(lower, upper, auto=None)
            else:
                self.set_xlim(upper, lower, auto=None)

    def get_xlim(self):
        """
        Get the x-axis range

        Returns
        -------
        xlimits : tuple
            Returns the current x-axis limits as the tuple
            (`left`, `right`).

        Notes
        -----
        The x-axis may be inverted, in which case the `left` value will
        be greater than the `right` value.

        """
        return tuple(self.viewLim.intervalx)

    def _validate_converted_limits(self, limit, convert):
        """
        Raise ValueError if converted limits are non-finite.

        Note that this function also accepts None as a limit argument.

        Returns
        -------
        The limit value after call to convert(), or None if limit is None.

        """
        if limit is not None:
            converted_limit = convert(limit)
            if (isinstance(converted_limit, float) and
                    (not np.isreal(converted_limit) or
                        not np.isfinite(converted_limit))):
                raise ValueError("Axis limits cannot be NaN or Inf")
            return converted_limit

    def set_xlim(self, left=None, right=None, emit=True, auto=False, **kw):
        """
        Set the data limits for the x-axis

        .. ACCEPTS: (left: float, right: float)

        Parameters
        ----------
        left : scalar, optional
            The left xlim (default: None, which leaves the left limit
            unchanged).

        right : scalar, optional
            The right xlim (default: None, which leaves the right limit
            unchanged).

        emit : bool, optional
            Whether to notify observers of limit change (default: True).

        auto : bool or None, optional
            Whether to turn on autoscaling of the x-axis. True turns on,
            False turns off (default action), None leaves unchanged.

        xlimits : tuple, optional
            The left and right xlims may be passed as the tuple
            (`left`, `right`) as the first positional argument (or as
            the `left` keyword argument).

        Returns
        -------
        xlimits : tuple
            Returns the new x-axis limits as (`left`, `right`).

        Notes
        -----
        The `left` value may be greater than the `right` value, in which
        case the x-axis values will decrease from left to right.

        Examples
        --------
        >>> set_xlim(left, right)
        >>> set_xlim((left, right))
        >>> left, right = set_xlim(left, right)

        One limit may be left unchanged.

        >>> set_xlim(right=right_lim)

        Limits may be passed in reverse order to flip the direction of
        the x-axis. For example, suppose `x` represents the number of
        years before present. The x-axis limits might be set like the
        following so 5000 years ago is on the left of the plot and the
        present is on the right.

        >>> set_xlim(5000, 0)

        """
        if 'xmin' in kw:
            left = kw.pop('xmin')
        if 'xmax' in kw:
            right = kw.pop('xmax')
        if kw:
            raise ValueError("unrecognized kwargs: %s" % list(kw))

        if right is None and iterable(left):
            left, right = left

        self._process_unit_info(xdata=(left, right))
        left = self._validate_converted_limits(left, self.convert_xunits)
        right = self._validate_converted_limits(right, self.convert_xunits)

        old_left, old_right = self.get_xlim()
        if left is None:
            left = old_left
        if right is None:
            right = old_right

        if left == right:
            warnings.warn(
                ('Attempting to set identical left==right results\n'
                 'in singular transformations; automatically expanding.\n'
                 'left=%s, right=%s') % (left, right))
        left, right = mtransforms.nonsingular(left, right, increasing=False)

        if self.get_xscale() == 'log' and (left <= 0.0 or right <= 0.0):
            warnings.warn(
                'Attempted to set non-positive xlimits for log-scale axis; '
                'invalid limits will be ignored.')
        left, right = self.xaxis.limit_range_for_scale(left, right)

        self.viewLim.intervalx = (left, right)
        if auto is not None:
            self._autoscaleXon = bool(auto)

        if emit:
            self.callbacks.process('xlim_changed', self)
            # Call all of the other x-axes that are shared with this one
            for other in self._shared_x_axes.get_siblings(self):
                if other is not self:
                    other.set_xlim(self.viewLim.intervalx,
                                   emit=False, auto=auto)
                    if (other.figure != self.figure and
                            other.figure.canvas is not None):
                        other.figure.canvas.draw_idle()
        self.stale = True
        return left, right

    def get_xscale(self):
        return self.xaxis.get_scale()
    get_xscale.__doc__ = "Return the xaxis scale string: %s""" % (
        ", ".join(mscale.get_scale_names()))

    def set_xscale(self, value, **kwargs):
        """
        Set the x-axis scale.

        .. ACCEPTS: [ 'linear' | 'log' | 'symlog' | 'logit' | ... ]

        Parameters
        ----------
        value : {"linear", "log", "symlog", "logit"}
            scaling strategy to apply

        Notes
        -----
        Different kwargs are accepted, depending on the scale. See
        the `~matplotlib.scale` module for more information.

        See also
        --------
        matplotlib.scale.LinearScale : linear transform

        matplotlib.scale.LogTransform : log transform

        matplotlib.scale.SymmetricalLogTransform : symlog transform

        matplotlib.scale.LogisticTransform : logit transform
        """
        g = self.get_shared_x_axes()
        for ax in g.get_siblings(self):
            ax.xaxis._set_scale(value, **kwargs)
            ax._update_transScale()
            ax.stale = True

        self.autoscale_view(scaley=False)

    def get_xticks(self, minor=False):
        """Return the x ticks as a list of locations"""
        return self.xaxis.get_ticklocs(minor=minor)

    def set_xticks(self, ticks, minor=False):
        """
        Set the x ticks with list of *ticks*

        .. ACCEPTS: list of tick locations.

        Parameters
        ----------
        ticks : list
            List of x-axis tick locations.

        minor : bool, optional
            If ``False`` sets major ticks, if ``True`` sets minor ticks.
            Default is ``False``.
        """
        ret = self.xaxis.set_ticks(ticks, minor=minor)
        self.stale = True
        return ret

    def get_xmajorticklabels(self):
        """
        Get the major x tick labels.

        Returns
        -------
        labels : list
            List of :class:`~matplotlib.text.Text` instances
        """
        return cbook.silent_list('Text xticklabel',
                                 self.xaxis.get_majorticklabels())

    def get_xminorticklabels(self):
        """
        Get the minor x tick labels.

        Returns
        -------
        labels : list
            List of :class:`~matplotlib.text.Text` instances
        """
        return cbook.silent_list('Text xticklabel',
                                 self.xaxis.get_minorticklabels())

    def get_xticklabels(self, minor=False, which=None):
        """
        Get the x tick labels as a list of :class:`~matplotlib.text.Text`
        instances.

        Parameters
        ----------
        minor : bool, optional
           If True return the minor ticklabels,
           else return the major ticklabels.

        which : None, ('minor', 'major', 'both')
           Overrides `minor`.

           Selects which ticklabels to return

        Returns
        -------
        ret : list
           List of :class:`~matplotlib.text.Text` instances.
        """
        return cbook.silent_list('Text xticklabel',
                                 self.xaxis.get_ticklabels(minor=minor,
                                                           which=which))

    def set_xticklabels(self, labels, fontdict=None, minor=False, **kwargs):
        """
        Set the x-tick labels with list of string labels.

        .. ACCEPTS: list of string labels

        Parameters
        ----------
        labels : list of str
            List of string labels.

        fontdict : dict, optional
            A dictionary controlling the appearance of the ticklabels.
            The default `fontdict` is::

               {'fontsize': rcParams['axes.titlesize'],
                'fontweight': rcParams['axes.titleweight'],
                'verticalalignment': 'baseline',
                'horizontalalignment': loc}

        minor : bool, optional
            Whether to set the minor ticklabels rather than the major ones.

        Returns
        -------
        A list of `~.text.Text` instances.

        Other Parameters
        -----------------
        **kwargs : `~.text.Text` properties.
        """
        if fontdict is not None:
            kwargs.update(fontdict)
        ret = self.xaxis.set_ticklabels(labels,
                                        minor=minor, **kwargs)
        self.stale = True
        return ret

    def invert_yaxis(self):
        """Invert the y-axis."""
        self.set_ylim(self.get_ylim()[::-1], auto=None)

    def yaxis_inverted(self):
        """Return whether the y-axis is inverted."""
        bottom, top = self.get_ylim()
        return top < bottom

    def get_ybound(self):
        """Return the lower and upper y-axis bounds, in increasing order."""
        bottom, top = self.get_ylim()
        if bottom < top:
            return bottom, top
        else:
            return top, bottom

    def set_ybound(self, lower=None, upper=None):
        """
        Set the lower and upper numerical bounds of the y-axis.
        This method will honor axes inversion regardless of parameter order.
        It will not change the _autoscaleYon attribute.

        .. ACCEPTS: (lower: float, upper: float)
        """
        if upper is None and iterable(lower):
            lower, upper = lower

        old_lower, old_upper = self.get_ybound()

        if lower is None:
            lower = old_lower
        if upper is None:
            upper = old_upper

        if self.yaxis_inverted():
            if lower < upper:
                self.set_ylim(upper, lower, auto=None)
            else:
                self.set_ylim(lower, upper, auto=None)
        else:
            if lower < upper:
                self.set_ylim(lower, upper, auto=None)
            else:
                self.set_ylim(upper, lower, auto=None)

    def get_ylim(self):
        """
        Get the y-axis range

        Returns
        -------
        ylimits : tuple
            Returns the current y-axis limits as the tuple
            (`bottom`, `top`).

        Notes
        -----
        The y-axis may be inverted, in which case the `bottom` value
        will be greater than the `top` value.

        """
        return tuple(self.viewLim.intervaly)

    def set_ylim(self, bottom=None, top=None, emit=True, auto=False, **kw):
        """
        Set the data limits for the y-axis

        .. ACCEPTS: (bottom: float, top: float)

        Parameters
        ----------
        bottom : scalar, optional
            The bottom ylim (default: None, which leaves the bottom
            limit unchanged).

        top : scalar, optional
            The top ylim (default: None, which leaves the top limit
            unchanged).

        emit : bool, optional
            Whether to notify observers of limit change (default: True).

        auto : bool or None, optional
            Whether to turn on autoscaling of the y-axis. True turns on,
            False turns off (default action), None leaves unchanged.

        ylimits : tuple, optional
            The bottom and top yxlims may be passed as the tuple
            (`bottom`, `top`) as the first positional argument (or as
            the `bottom` keyword argument).

        Returns
        -------
        ylimits : tuple
            Returns the new y-axis limits as (`bottom`, `top`).

        Notes
        -----
        The `bottom` value may be greater than the `top` value, in which
        case the y-axis values will decrease from bottom to top.

        Examples
        --------
        >>> set_ylim(bottom, top)
        >>> set_ylim((bottom, top))
        >>> bottom, top = set_ylim(bottom, top)

        One limit may be left unchanged.

        >>> set_ylim(top=top_lim)

        Limits may be passed in reverse order to flip the direction of
        the y-axis. For example, suppose `y` represents depth of the
        ocean in m. The y-axis limits might be set like the following
        so 5000 m depth is at the bottom of the plot and the surface,
        0 m, is at the top.

        >>> set_ylim(5000, 0)
        """
        if 'ymin' in kw:
            bottom = kw.pop('ymin')
        if 'ymax' in kw:
            top = kw.pop('ymax')
        if kw:
            raise ValueError("unrecognized kwargs: %s" % list(kw))

        if top is None and iterable(bottom):
            bottom, top = bottom

        bottom = self._validate_converted_limits(bottom, self.convert_yunits)
        top = self._validate_converted_limits(top, self.convert_yunits)

        old_bottom, old_top = self.get_ylim()

        if bottom is None:
            bottom = old_bottom
        if top is None:
            top = old_top

        if bottom == top:
            warnings.warn(
                ('Attempting to set identical bottom==top results\n'
                 'in singular transformations; automatically expanding.\n'
                 'bottom=%s, top=%s') % (bottom, top))

        bottom, top = mtransforms.nonsingular(bottom, top, increasing=False)

        if self.get_yscale() == 'log' and (bottom <= 0.0 or top <= 0.0):
            warnings.warn(
                'Attempted to set non-positive ylimits for log-scale axis; '
                'invalid limits will be ignored.')
        bottom, top = self.yaxis.limit_range_for_scale(bottom, top)

        self.viewLim.intervaly = (bottom, top)
        if auto is not None:
            self._autoscaleYon = bool(auto)

        if emit:
            self.callbacks.process('ylim_changed', self)
            # Call all of the other y-axes that are shared with this one
            for other in self._shared_y_axes.get_siblings(self):
                if other is not self:
                    other.set_ylim(self.viewLim.intervaly,
                                   emit=False, auto=auto)
                    if (other.figure != self.figure and
                            other.figure.canvas is not None):
                        other.figure.canvas.draw_idle()
        self.stale = True
        return bottom, top

    def get_yscale(self):
        return self.yaxis.get_scale()
    get_yscale.__doc__ = "Return the yaxis scale string: %s""" % (
        ", ".join(mscale.get_scale_names()))

    def set_yscale(self, value, **kwargs):
        """
        Set the y-axis scale.

        .. ACCEPTS: [ 'linear' | 'log' | 'symlog' | 'logit' | ... ]

        Parameters
        ----------
        value : {"linear", "log", "symlog", "logit"}
            scaling strategy to apply

        Notes
        -----
        Different kwargs are accepted, depending on the scale. See
        the `~matplotlib.scale` module for more information.

        See also
        --------
        matplotlib.scale.LinearScale : linear transform

        matplotlib.scale.LogTransform : log transform

        matplotlib.scale.SymmetricalLogTransform : symlog transform

        matplotlib.scale.LogisticTransform : logit transform
        """
        g = self.get_shared_y_axes()
        for ax in g.get_siblings(self):
            ax.yaxis._set_scale(value, **kwargs)
            ax._update_transScale()
            ax.stale = True
        self.autoscale_view(scalex=False)

    def get_yticks(self, minor=False):
        """Return the y ticks as a list of locations"""
        return self.yaxis.get_ticklocs(minor=minor)

    def set_yticks(self, ticks, minor=False):
        """
        Set the y ticks with list of *ticks*

        .. ACCEPTS: list of tick locations.

        Parameters
        ----------
        ticks : sequence
            List of y-axis tick locations

        minor : bool, optional
            If ``False`` sets major ticks, if ``True`` sets minor ticks.
            Default is ``False``.
        """
        ret = self.yaxis.set_ticks(ticks, minor=minor)
        return ret

    def get_ymajorticklabels(self):
        """
        Get the major y tick labels.

        Returns
        -------
        labels : list
            List of :class:`~matplotlib.text.Text` instances
        """
        return cbook.silent_list('Text yticklabel',
                                 self.yaxis.get_majorticklabels())

    def get_yminorticklabels(self):
        """
        Get the minor y tick labels.

        Returns
        -------
        labels : list
            List of :class:`~matplotlib.text.Text` instances
        """
        return cbook.silent_list('Text yticklabel',
                                 self.yaxis.get_minorticklabels())

    def get_yticklabels(self, minor=False, which=None):
        """
        Get the x tick labels as a list of :class:`~matplotlib.text.Text`
        instances.

        Parameters
        ----------
        minor : bool
           If True return the minor ticklabels,
           else return the major ticklabels

        which : None, ('minor', 'major', 'both')
           Overrides `minor`.

           Selects which ticklabels to return

        Returns
        -------
        ret : list
           List of :class:`~matplotlib.text.Text` instances.
        """
        return cbook.silent_list('Text yticklabel',
                                 self.yaxis.get_ticklabels(minor=minor,
                                                           which=which))

    def set_yticklabels(self, labels, fontdict=None, minor=False, **kwargs):
        """
        Set the y-tick labels with list of strings labels.

        .. ACCEPTS: list of string labels

        Parameters
        ----------
        labels : list of str
            list of string labels

        fontdict : dict, optional
            A dictionary controlling the appearance of the ticklabels.
            The default `fontdict` is::

               {'fontsize': rcParams['axes.titlesize'],
                'fontweight': rcParams['axes.titleweight'],
                'verticalalignment': 'baseline',
                'horizontalalignment': loc}

        minor : bool, optional
            Whether to set the minor ticklabels rather than the major ones.

        Returns
        -------
        A list of `~.text.Text` instances.

        Other Parameters
        ----------------
        **kwargs : `~.text.Text` properties.
        """
        if fontdict is not None:
            kwargs.update(fontdict)
        return self.yaxis.set_ticklabels(labels,
                                         minor=minor, **kwargs)

    def xaxis_date(self, tz=None):
        """
        Sets up x-axis ticks and labels that treat the x data as dates.

        Parameters
        ----------
        tz : string or :class:`tzinfo` instance, optional
            Timezone string or timezone. Defaults to rc value.
        """
        # should be enough to inform the unit conversion interface
        # dates are coming in
        self.xaxis.axis_date(tz)

    def yaxis_date(self, tz=None):
        """
        Sets up y-axis ticks and labels that treat the y data as dates.

        Parameters
        ----------
        tz : string or :class:`tzinfo` instance, optional
            Timezone string or timezone. Defaults to rc value.
        """
        self.yaxis.axis_date(tz)

    def format_xdata(self, x):
        """
        Return *x* string formatted.  This function will use the attribute
        self.fmt_xdata if it is callable, else will fall back on the xaxis
        major formatter
        """
        try:
            return self.fmt_xdata(x)
        except TypeError:
            func = self.xaxis.get_major_formatter().format_data_short
            val = func(x)
            return val

    def format_ydata(self, y):
        """
        Return y string formatted.  This function will use the
        :attr:`fmt_ydata` attribute if it is callable, else will fall
        back on the yaxis major formatter
        """
        try:
            return self.fmt_ydata(y)
        except TypeError:
            func = self.yaxis.get_major_formatter().format_data_short
            val = func(y)
            return val

    def format_coord(self, x, y):
        """Return a format string formatting the *x*, *y* coord"""
        if x is None:
            xs = '???'
        else:
            xs = self.format_xdata(x)
        if y is None:
            ys = '???'
        else:
            ys = self.format_ydata(y)
        return 'x=%s y=%s' % (xs, ys)

    def minorticks_on(self):
        'Add autoscaling minor ticks to the axes.'
        for ax in (self.xaxis, self.yaxis):
            scale = ax.get_scale()
            if scale == 'log':
                s = ax._scale
                ax.set_minor_locator(mticker.LogLocator(s.base, s.subs))
            elif scale == 'symlog':
                s = ax._scale
                ax.set_minor_locator(
                    mticker.SymmetricalLogLocator(s._transform, s.subs))
            else:
                ax.set_minor_locator(mticker.AutoMinorLocator())

    def minorticks_off(self):
        """Remove minor ticks from the axes."""
        self.xaxis.set_minor_locator(mticker.NullLocator())
        self.yaxis.set_minor_locator(mticker.NullLocator())

    # Interactive manipulation

    def can_zoom(self):
        """
        Return *True* if this axes supports the zoom box button functionality.
        """
        return True

    def can_pan(self):
        """
        Return *True* if this axes supports any pan/zoom button functionality.
        """
        return True

    def get_navigate(self):
        """
        Get whether the axes responds to navigation commands
        """
        return self._navigate

    def set_navigate(self, b):
        """
        Set whether the axes responds to navigation toolbar commands

        .. ACCEPTS: bool

        Parameters
        ----------
        b : bool
        """
        self._navigate = b

    def get_navigate_mode(self):
        """
        Get the navigation toolbar button status: 'PAN', 'ZOOM', or None
        """
        return self._navigate_mode

    def set_navigate_mode(self, b):
        """
        Set the navigation toolbar button status;

        .. warning::
            this is not a user-API function.

        """
        self._navigate_mode = b

    def _get_view(self):
        """
        Save information required to reproduce the current view.

        Called before a view is changed, such as during a pan or zoom
        initiated by the user. You may return any information you deem
        necessary to describe the view.

        .. note::

            Intended to be overridden by new projection types, but if not, the
            default implementation saves the view limits. You *must* implement
            :meth:`_set_view` if you implement this method.
        """
        xmin, xmax = self.get_xlim()
        ymin, ymax = self.get_ylim()
        return (xmin, xmax, ymin, ymax)

    def _set_view(self, view):
        """
        Apply a previously saved view.

        Called when restoring a view, such as with the navigation buttons.

        .. note::

            Intended to be overridden by new projection types, but if not, the
            default implementation restores the view limits. You *must*
            implement :meth:`_get_view` if you implement this method.
        """
        xmin, xmax, ymin, ymax = view
        self.set_xlim((xmin, xmax))
        self.set_ylim((ymin, ymax))

    def _set_view_from_bbox(self, bbox, direction='in',
                            mode=None, twinx=False, twiny=False):
        """
        Update view from a selection bbox.

        .. note::

            Intended to be overridden by new projection types, but if not, the
            default implementation sets the view limits to the bbox directly.

        Parameters
        ----------

        bbox : 4-tuple or 3 tuple
            * If bbox is a 4 tuple, it is the selected bounding box limits,
                in *display* coordinates.
            * If bbox is a 3 tuple, it is an (xp, yp, scl) triple, where
                (xp,yp) is the center of zooming and scl the scale factor to
                zoom by.

        direction : str
            The direction to apply the bounding box.
                * `'in'` - The bounding box describes the view directly, i.e.,
                           it zooms in.
                * `'out'` - The bounding box describes the size to make the
                            existing view, i.e., it zooms out.

        mode : str or None
            The selection mode, whether to apply the bounding box in only the
            `'x'` direction, `'y'` direction or both (`None`).

        twinx : bool
            Whether this axis is twinned in the *x*-direction.

        twiny : bool
            Whether this axis is twinned in the *y*-direction.
        """
        Xmin, Xmax = self.get_xlim()
        Ymin, Ymax = self.get_ylim()

        if len(bbox) == 3:
            # Zooming code
            xp, yp, scl = bbox

            # Should not happen
            if scl == 0:
                scl = 1.

            # direction = 'in'
            if scl > 1:
                direction = 'in'
            else:
                direction = 'out'
                scl = 1/scl

            # get the limits of the axes
            tranD2C = self.transData.transform
            xmin, ymin = tranD2C((Xmin, Ymin))
            xmax, ymax = tranD2C((Xmax, Ymax))

            # set the range
            xwidth = xmax - xmin
            ywidth = ymax - ymin
            xcen = (xmax + xmin)*.5
            ycen = (ymax + ymin)*.5
            xzc = (xp*(scl - 1) + xcen)/scl
            yzc = (yp*(scl - 1) + ycen)/scl

            bbox = [xzc - xwidth/2./scl, yzc - ywidth/2./scl,
                    xzc + xwidth/2./scl, yzc + ywidth/2./scl]
        elif len(bbox) != 4:
            # should be len 3 or 4 but nothing else
            warnings.warn(
                "Warning in _set_view_from_bbox: bounding box is not a tuple "
                "of length 3 or 4. Ignoring the view change.")
            return

        # Just grab bounding box
        lastx, lasty, x, y = bbox

        # zoom to rect
        inverse = self.transData.inverted()
        lastx, lasty = inverse.transform_point((lastx, lasty))
        x, y = inverse.transform_point((x, y))

        if twinx:
            x0, x1 = Xmin, Xmax
        else:
            if Xmin < Xmax:
                if x < lastx:
                    x0, x1 = x, lastx
                else:
                    x0, x1 = lastx, x
                if x0 < Xmin:
                    x0 = Xmin
                if x1 > Xmax:
                    x1 = Xmax
            else:
                if x > lastx:
                    x0, x1 = x, lastx
                else:
                    x0, x1 = lastx, x
                if x0 > Xmin:
                    x0 = Xmin
                if x1 < Xmax:
                    x1 = Xmax

        if twiny:
            y0, y1 = Ymin, Ymax
        else:
            if Ymin < Ymax:
                if y < lasty:
                    y0, y1 = y, lasty
                else:
                    y0, y1 = lasty, y
                if y0 < Ymin:
                    y0 = Ymin
                if y1 > Ymax:
                    y1 = Ymax
            else:
                if y > lasty:
                    y0, y1 = y, lasty
                else:
                    y0, y1 = lasty, y
                if y0 > Ymin:
                    y0 = Ymin
                if y1 < Ymax:
                    y1 = Ymax

        if direction == 'in':
            if mode == 'x':
                self.set_xlim((x0, x1))
            elif mode == 'y':
                self.set_ylim((y0, y1))
            else:
                self.set_xlim((x0, x1))
                self.set_ylim((y0, y1))
        elif direction == 'out':
            if self.get_xscale() == 'log':
                alpha = np.log(Xmax / Xmin) / np.log(x1 / x0)
                rx1 = pow(Xmin / x0, alpha) * Xmin
                rx2 = pow(Xmax / x0, alpha) * Xmin
            else:
                alpha = (Xmax - Xmin) / (x1 - x0)
                rx1 = alpha * (Xmin - x0) + Xmin
                rx2 = alpha * (Xmax - x0) + Xmin
            if self.get_yscale() == 'log':
                alpha = np.log(Ymax / Ymin) / np.log(y1 / y0)
                ry1 = pow(Ymin / y0, alpha) * Ymin
                ry2 = pow(Ymax / y0, alpha) * Ymin
            else:
                alpha = (Ymax - Ymin) / (y1 - y0)
                ry1 = alpha * (Ymin - y0) + Ymin
                ry2 = alpha * (Ymax - y0) + Ymin

            if mode == 'x':
                self.set_xlim((rx1, rx2))
            elif mode == 'y':
                self.set_ylim((ry1, ry2))
            else:
                self.set_xlim((rx1, rx2))
                self.set_ylim((ry1, ry2))

    def start_pan(self, x, y, button):
        """
        Called when a pan operation has started.

        *x*, *y* are the mouse coordinates in display coords.
        button is the mouse button number:

        * 1: LEFT
        * 2: MIDDLE
        * 3: RIGHT

        .. note::

            Intended to be overridden by new projection types.

        """
        self._pan_start = cbook.Bunch(
            lim=self.viewLim.frozen(),
            trans=self.transData.frozen(),
            trans_inverse=self.transData.inverted().frozen(),
            bbox=self.bbox.frozen(),
            x=x,
            y=y)

    def end_pan(self):
        """
        Called when a pan operation completes (when the mouse button
        is up.)

        .. note::

            Intended to be overridden by new projection types.

        """
        del self._pan_start

    def drag_pan(self, button, key, x, y):
        """
        Called when the mouse moves during a pan operation.

        *button* is the mouse button number:

        * 1: LEFT
        * 2: MIDDLE
        * 3: RIGHT

        *key* is a "shift" key

        *x*, *y* are the mouse coordinates in display coords.

        .. note::

            Intended to be overridden by new projection types.

        """
        def format_deltas(key, dx, dy):
            if key == 'control':
                if abs(dx) > abs(dy):
                    dy = dx
                else:
                    dx = dy
            elif key == 'x':
                dy = 0
            elif key == 'y':
                dx = 0
            elif key == 'shift':
                if 2 * abs(dx) < abs(dy):
                    dx = 0
                elif 2 * abs(dy) < abs(dx):
                    dy = 0
                elif abs(dx) > abs(dy):
                    dy = dy / abs(dy) * abs(dx)
                else:
                    dx = dx / abs(dx) * abs(dy)
            return dx, dy

        p = self._pan_start
        dx = x - p.x
        dy = y - p.y
        if dx == 0 and dy == 0:
            return
        if button == 1:
            dx, dy = format_deltas(key, dx, dy)
            result = p.bbox.translated(-dx, -dy).transformed(p.trans_inverse)
        elif button == 3:
            try:
                dx = -dx / self.bbox.width
                dy = -dy / self.bbox.height
                dx, dy = format_deltas(key, dx, dy)
                if self.get_aspect() != 'auto':
                    dx = dy = 0.5 * (dx + dy)
                alpha = np.power(10.0, (dx, dy))
                start = np.array([p.x, p.y])
                oldpoints = p.lim.transformed(p.trans)
                newpoints = start + alpha * (oldpoints - start)
                result = (mtransforms.Bbox(newpoints)
                          .transformed(p.trans_inverse))
            except OverflowError:
                warnings.warn('Overflow while panning')
                return

        valid = np.isfinite(result.transformed(p.trans))
        points = result.get_points().astype(object)
        # Just ignore invalid limits (typically, underflow in log-scale).
        points[~valid] = None
        self.set_xlim(points[:, 0])
        self.set_ylim(points[:, 1])

    @cbook.deprecated("2.1")
    def get_cursor_props(self):
        """
        Return the cursor propertiess as a (*linewidth*, *color*)
        tuple, where *linewidth* is a float and *color* is an RGBA
        tuple
        """
        return self._cursorProps

    @cbook.deprecated("2.1")
    def set_cursor_props(self, *args):
        """Set the cursor property as

        Call signature ::

          ax.set_cursor_props(linewidth, color)

        or::

          ax.set_cursor_props((linewidth, color))

        ACCEPTS: a (*float*, *color*) tuple
        """
        if len(args) == 1:
            lw, c = args[0]
        elif len(args) == 2:
            lw, c = args
        else:
            raise ValueError('args must be a (linewidth, color) tuple')
        c = mcolors.to_rgba(c)
        self._cursorProps = lw, c

    def get_children(self):
        """return a list of child artists"""
        children = []
        children.extend(self.collections)
        children.extend(self.patches)
        children.extend(self.lines)
        children.extend(self.texts)
        children.extend(self.artists)
        children.extend(six.itervalues(self.spines))
        children.append(self.xaxis)
        children.append(self.yaxis)
        children.append(self.title)
        children.append(self._left_title)
        children.append(self._right_title)
        children.extend(self.tables)
        children.extend(self.images)
        if self.legend_ is not None:
            children.append(self.legend_)
        children.append(self.patch)
        return children

    def contains(self, mouseevent):
        """
        Test whether the mouse event occurred in the axes.

        Returns *True* / *False*, {}
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)
        return self.patch.contains(mouseevent)

    def contains_point(self, point):
        """
        Returns *True* if the point (tuple of x,y) is inside the axes
        (the area defined by the its patch). A pixel coordinate is
        required.

        """
        return self.patch.contains_point(point, radius=1.0)

    def pick(self, *args):
        """Trigger pick event

        Call signature::

            pick(mouseevent)

        each child artist will fire a pick event if mouseevent is over
        the artist and the artist has picker set
        """
        martist.Artist.pick(self, args[0])

    def get_default_bbox_extra_artists(self):
        return [artist for artist in self.get_children()
                if artist.get_visible()]

    def get_tightbbox(self, renderer, call_axes_locator=True):
        """
        Return the tight bounding box of the axes.
        The dimension of the Bbox in canvas coordinate.

        If *call_axes_locator* is *False*, it does not call the
        _axes_locator attribute, which is necessary to get the correct
        bounding box. ``call_axes_locator==False`` can be used if the
        caller is only intereted in the relative size of the tightbbox
        compared to the axes bbox.
        """

        bb = []

        if not self.get_visible():
            return None

        locator = self.get_axes_locator()
        if locator and call_axes_locator:
            pos = locator(self, renderer)
            self.apply_aspect(pos)
        else:
            self.apply_aspect()

        bb.append(self.get_window_extent(renderer))

        if self.title.get_visible():
            bb.append(self.title.get_window_extent(renderer))
        if self._left_title.get_visible():
            bb.append(self._left_title.get_window_extent(renderer))
        if self._right_title.get_visible():
            bb.append(self._right_title.get_window_extent(renderer))

        bb_xaxis = self.xaxis.get_tightbbox(renderer)
        if bb_xaxis:
            bb.append(bb_xaxis)

        bb_yaxis = self.yaxis.get_tightbbox(renderer)
        if bb_yaxis:
            bb.append(bb_yaxis)

        for child in self.get_children():
            if isinstance(child, OffsetBox) and child.get_visible():
                bb.append(child.get_window_extent(renderer))
            elif isinstance(child, Legend) and child.get_visible():
                bb.append(child._legend_box.get_window_extent(renderer))

        _bbox = mtransforms.Bbox.union(
            [b for b in bb if b.width != 0 or b.height != 0])

        return _bbox

    def _make_twin_axes(self, *kl, **kwargs):
        """
        Make a twinx axes of self. This is used for twinx and twiny.
        """
        # Typically, SubplotBase._make_twin_axes is called instead of this.
        # There is also an override in axes_grid1/axes_divider.py.
        if 'sharex' in kwargs and 'sharey' in kwargs:
            raise ValueError("Twinned Axes may share only one axis.")
        ax2 = self.figure.add_axes(self.get_position(True), *kl, **kwargs)
        self.set_adjustable('datalim')
        ax2.set_adjustable('datalim')
        self._twinned_axes.join(self, ax2)
        return ax2

    def twinx(self):
        """
        Create a twin Axes sharing the xaxis

        Create a new Axes instance with an invisible x-axis and an independent
        y-axis positioned opposite to the original one (i.e. at right). The
        x-axis autoscale setting will be inherited from the original Axes.
        To ensure that the tick marks of both y-axes align, see
        `~matplotlib.ticker.LinearLocator`

        Returns
        -------
        ax_twin : Axes
            The newly created Axes instance

        Notes
        -----
        For those who are 'picking' artists while using twinx, pick
        events are only called for the artists in the top-most axes.
        """
        ax2 = self._make_twin_axes(sharex=self)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position('right')
        ax2.yaxis.set_offset_position('right')
        ax2.set_autoscalex_on(self.get_autoscalex_on())
        self.yaxis.tick_left()
        ax2.xaxis.set_visible(False)
        ax2.patch.set_visible(False)
        return ax2

    def twiny(self):
        """
        Create a twin Axes sharing the yaxis

        Create a new Axes instance with an invisible y-axis and an independent
        x-axis positioned opposite to the original one (i.e. at top). The
        y-axis autoscale setting will be inherited from the original Axes.
        To ensure that the tick marks of both x-axes align, see
        `~matplotlib.ticker.LinearLocator`

        Returns
        -------
        ax_twin : Axes
            The newly created Axes instance

        Notes
        -----
        For those who are 'picking' artists while using twiny, pick
        events are only called for the artists in the top-most axes.
        """

        ax2 = self._make_twin_axes(sharey=self)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')
        ax2.set_autoscaley_on(self.get_autoscaley_on())
        self.xaxis.tick_bottom()
        ax2.yaxis.set_visible(False)
        ax2.patch.set_visible(False)
        return ax2

    def get_shared_x_axes(self):
        """Return a reference to the shared axes Grouper object for x axes."""
        return self._shared_x_axes

    def get_shared_y_axes(self):
        """Return a reference to the shared axes Grouper object for y axes."""
        return self._shared_y_axes
