from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange, zip, zip_longest

import functools
import itertools
import logging
import math
import warnings

import numpy as np
from numpy import ma

import matplotlib
from matplotlib import _preprocess_data

import matplotlib.cbook as cbook
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.contour as mcontour
import matplotlib.category as _  # <-registers a category unit converter
import matplotlib.dates as _  # <-registers a date unit converter
import matplotlib.docstring as docstring
import matplotlib.image as mimage
import matplotlib.legend as mlegend
import matplotlib.lines as mlines
import matplotlib.markers as mmarkers
import matplotlib.mlab as mlab
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.quiver as mquiver
import matplotlib.stackplot as mstack
import matplotlib.streamplot as mstream
import matplotlib.table as mtable
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import matplotlib.tri as mtri
from matplotlib.cbook import (
    _backports, mplDeprecation, warn_deprecated,
    STEP_LOOKUP_MAP, iterable, safe_first_element)
from matplotlib.container import BarContainer, ErrorbarContainer, StemContainer
from matplotlib.axes._base import _AxesBase, _process_plot_format

_log = logging.getLogger(__name__)

rcParams = matplotlib.rcParams

_alias_map = {'color': ['c'],
              'linewidth': ['lw'],
              'linestyle': ['ls'],
              'facecolor': ['fc'],
              'edgecolor': ['ec'],
              'markerfacecolor': ['mfc'],
              'markeredgecolor': ['mec'],
              'markeredgewidth': ['mew'],
              'markersize': ['ms'],
             }


def _plot_args_replacer(args, data):
    if len(args) == 1:
        return ["y"]
    elif len(args) == 2:
        # this can be two cases: x,y or y,c
        if not args[1] in data:
            # this is not in data, so just assume that it is something which
            # will not get replaced (color spec or array like).
            return ["y", "c"]
        # it's data, but could be a color code like 'ro' or 'b--'
        # -> warn the user in that case...
        try:
            _process_plot_format(args[1])
        except ValueError:
            pass
        else:
            warnings.warn(
                "Second argument {!r} is ambiguous: could be a color spec but "
                "is in data; using as data.  Either rename the entry in data "
                "or use three arguments to plot.".format(args[1]),
                RuntimeWarning, stacklevel=3)
        return ["x", "y"]
    elif len(args) == 3:
        return ["x", "y", "c"]
    else:
        raise ValueError("Using arbitrary long args with data is not "
                         "supported due to ambiguity of arguments.\nUse "
                         "multiple plotting calls instead.")


# The axes module contains all the wrappers to plotting functions.
# All the other methods should go in the _AxesBase class.

class Axes(_AxesBase):
    """
    The :class:`Axes` contains most of the figure elements:
    :class:`~matplotlib.axis.Axis`, :class:`~matplotlib.axis.Tick`,
    :class:`~matplotlib.lines.Line2D`, :class:`~matplotlib.text.Text`,
    :class:`~matplotlib.patches.Polygon`, etc., and sets the
    coordinate system.

    The :class:`Axes` instance supports callbacks through a callbacks
    attribute which is a :class:`~matplotlib.cbook.CallbackRegistry`
    instance.  The events you can connect to are 'xlim_changed' and
    'ylim_changed' and the callback will be called with func(*ax*)
    where *ax* is the :class:`Axes` instance.
    """
    ### Labelling, legend and texts

    aname = 'Axes'

    def get_title(self, loc="center"):
        """
        Get an axes title.

        Get one of the three available axes titles. The available titles
        are positioned above the axes in the center, flush with the left
        edge, and flush with the right edge.

        Parameters
        ----------
        loc : {'center', 'left', 'right'}, str, optional
            Which title to get, defaults to 'center'.

        Returns
        -------
        title : str
            The title text string.

        """
        try:
            title = {'left': self._left_title,
                     'center': self.title,
                     'right': self._right_title}[loc.lower()]
        except KeyError:
            raise ValueError("'%s' is not a valid location" % loc)
        return title.get_text()

    def set_title(self, label, fontdict=None, loc="center", pad=None,
                    **kwargs):
        """
        Set a title for the axes.

        Set one of the three available axes titles. The available titles
        are positioned above the axes in the center, flush with the left
        edge, and flush with the right edge.

        Parameters
        ----------
        label : str
            Text to use for the title

        fontdict : dict
            A dictionary controlling the appearance of the title text,
            the default `fontdict` is::

               {'fontsize': rcParams['axes.titlesize'],
                'fontweight' : rcParams['axes.titleweight'],
                'verticalalignment': 'baseline',
                'horizontalalignment': loc}

        loc : {'center', 'left', 'right'}, str, optional
            Which title to set, defaults to 'center'

        pad : float
            The offset of the title from the top of the axes, in points.
            Default is ``None`` to use rcParams['axes.titlepad'].

        Returns
        -------
        text : :class:`~matplotlib.text.Text`
            The matplotlib text instance representing the title

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.text.Text` properties
            Other keyword arguments are text properties, see
            :class:`~matplotlib.text.Text` for a list of valid text
            properties.
        """
        try:
            title = {'left': self._left_title,
                     'center': self.title,
                     'right': self._right_title}[loc.lower()]
        except KeyError:
            raise ValueError("'%s' is not a valid location" % loc)
        default = {
            'fontsize': rcParams['axes.titlesize'],
            'fontweight': rcParams['axes.titleweight'],
            'verticalalignment': 'baseline',
            'horizontalalignment': loc.lower()}
        if pad is None:
            pad = rcParams['axes.titlepad']
        self._set_title_offset_trans(float(pad))
        title.set_text(label)
        title.update(default)
        if fontdict is not None:
            title.update(fontdict)
        title.update(kwargs)
        return title

    def get_xlabel(self):
        """
        Get the xlabel text string.
        """
        label = self.xaxis.get_label()
        return label.get_text()

    def set_xlabel(self, xlabel, fontdict=None, labelpad=None, **kwargs):
        """
        Set the label for the x-axis.

        Parameters
        ----------
        xlabel : str
            The label text.

        labelpad : scalar, optional, default: None
            Spacing in points between the label and the x-axis.

        Other Parameters
        ----------------
        **kwargs : `.Text` properties
            `.Text` properties control the appearance of the label.

        See also
        --------
        text : for information on how override and the optional args work
        """
        if labelpad is not None:
            self.xaxis.labelpad = labelpad
        return self.xaxis.set_label_text(xlabel, fontdict, **kwargs)

    def get_ylabel(self):
        """
        Get the ylabel text string.
        """
        label = self.yaxis.get_label()
        return label.get_text()

    def set_ylabel(self, ylabel, fontdict=None, labelpad=None, **kwargs):
        """
        Set the label for the y-axis.

        Parameters
        ----------
        ylabel : str
            The label text.

        labelpad : scalar, optional, default: None
            Spacing in points between the label and the y-axis.

        Other Parameters
        ----------------
        **kwargs : `.Text` properties
            `.Text` properties control the appearance of the label.

        See also
        --------
        text : for information on how override and the optional args work

        """
        if labelpad is not None:
            self.yaxis.labelpad = labelpad
        return self.yaxis.set_label_text(ylabel, fontdict, **kwargs)

    def get_legend_handles_labels(self, legend_handler_map=None):
        """
        Return handles and labels for legend

        ``ax.legend()`` is equivalent to ::

          h, l = ax.get_legend_handles_labels()
          ax.legend(h, l)

        """

        # pass through to legend.
        handles, labels = mlegend._get_legend_handles_labels([self],
                legend_handler_map)
        return handles, labels

    @docstring.dedent_interpd
    def legend(self, *args, **kwargs):
        """
        Places a legend on the axes.

        Call signatures::

            legend()
            legend(labels)
            legend(handles, labels)

        The call signatures correspond to three different ways how to use
        this method.

        **1. Automatic detection of elements to be shown in the legend**

        The elements to be added to the legend are automatically determined,
        when you do not pass in any extra arguments.

        In this case, the labels are taken from the artist. You can specify
        them either at artist creation or by calling the
        :meth:`~.Artist.set_label` method on the artist::

            line, = ax.plot([1, 2, 3], label='Inline label')
            ax.legend()

        or::

            line.set_label('Label via method')
            line, = ax.plot([1, 2, 3])
            ax.legend()

        Specific lines can be excluded from the automatic legend element
        selection by defining a label starting with an underscore.
        This is default for all artists, so calling `Axes.legend` without
        any arguments and without setting the labels manually will result in
        no legend being drawn.


        **2. Labeling existing plot elements**

        To make a legend for lines which already exist on the axes
        (via plot for instance), simply call this function with an iterable
        of strings, one for each legend item. For example::

            ax.plot([1, 2, 3])
            ax.legend(['A simple line'])

        Note: This way of using is discouraged, because the relation between
        plot elements and labels is only implicit by their order and can
        easily be mixed up.


        **3. Explicitly defining the elements in the legend**

        For full control of which artists have a legend entry, it is possible
        to pass an iterable of legend artists followed by an iterable of
        legend labels respectively::

            legend((line1, line2, line3), ('label1', 'label2', 'label3'))

        Parameters
        ----------

        handles : sequence of `.Artist`, optional
            A list of Artists (lines, patches) to be added to the legend.
            Use this together with *labels*, if you need full control on what
            is shown in the legend and the automatic mechanism described above
            is not sufficient.

            The length of handles and labels should be the same in this
            case. If they are not, they are truncated to the smaller length.

        labels : sequence of strings, optional
            A list of labels to show next to the artists.
            Use this together with *handles*, if you need full control on what
            is shown in the legend and the automatic mechanism described above
            is not sufficient.

        Other Parameters
        ----------------

        loc : int or string or pair of floats, default: 'upper right'
            The location of the legend. Possible codes are:

                ===============   =============
                Location String   Location Code
                ===============   =============
                'best'            0
                'upper right'     1
                'upper left'      2
                'lower left'      3
                'lower right'     4
                'right'           5
                'center left'     6
                'center right'    7
                'lower center'    8
                'upper center'    9
                'center'          10
                ===============   =============


            Alternatively can be a 2-tuple giving ``x, y`` of the lower-left
            corner of the legend in axes coordinates (in which case
            ``bbox_to_anchor`` will be ignored).

        bbox_to_anchor : `.BboxBase` or pair of floats
            Specify any arbitrary location for the legend in `bbox_transform`
            coordinates (default Axes coordinates).

            For example, to put the legend's upper right hand corner in the
            center of the axes the following keywords can be used::

               loc='upper right', bbox_to_anchor=(0.5, 0.5)

        ncol : integer
            The number of columns that the legend has. Default is 1.

        prop : None or :class:`matplotlib.font_manager.FontProperties` or dict
            The font properties of the legend. If None (default), the current
            :data:`matplotlib.rcParams` will be used.

        fontsize : int or float or {'xx-small', 'x-small', 'small', 'medium', \
'large', 'x-large', 'xx-large'}
            Controls the font size of the legend. If the value is numeric the
            size will be the absolute font size in points. String values are
            relative to the current default font size. This argument is only
            used if `prop` is not specified.

        numpoints : None or int
            The number of marker points in the legend when creating a legend
            entry for a `.Line2D` (line).
            Default is ``None``, which will take the value from
            :rc:`legend.numpoints`.

        scatterpoints : None or int
            The number of marker points in the legend when creating
            a legend entry for a `.PathCollection` (scatter plot).
            Default is ``None``, which will take the value from
            :rc:`legend.scatterpoints`.

        scatteryoffsets : iterable of floats
            The vertical offset (relative to the font size) for the markers
            created for a scatter plot legend entry. 0.0 is at the base the
            legend text, and 1.0 is at the top. To draw all markers at the
            same height, set to ``[0.5]``. Default is ``[0.375, 0.5, 0.3125]``.

        markerscale : None or int or float
            The relative size of legend markers compared with the originally
            drawn ones.
            Default is ``None``, which will take the value from
            :rc:`legend.markerscale`.

        markerfirst : bool
            If *True*, legend marker is placed to the left of the legend label.
            If *False*, legend marker is placed to the right of the legend
            label.
            Default is *True*.

        frameon : None or bool
            Control whether the legend should be drawn on a patch
            (frame).
            Default is ``None``, which will take the value from
            :rc:`legend.frameon`.

        fancybox : None or bool
            Control whether round edges should be enabled around the
            :class:`~matplotlib.patches.FancyBboxPatch` which makes up the
            legend's background.
            Default is ``None``, which will take the value from
            :rc:`legend.fancybox`.

        shadow : None or bool
            Control whether to draw a shadow behind the legend.
            Default is ``None``, which will take the value from
            :rc:`legend.shadow`.

        framealpha : None or float
            Control the alpha transparency of the legend's background.
            Default is ``None``, which will take the value from
            :rc:`legend.framealpha`.  If shadow is activated and
            *framealpha* is ``None``, the default value is ignored.

        facecolor : None or "inherit" or a color spec
            Control the legend's background color.
            Default is ``None``, which will take the value from
            :rc:`legend.facecolor`.  If ``"inherit"``, it will take
            :rc:`axes.facecolor`.

        edgecolor : None or "inherit" or a color spec
            Control the legend's background patch edge color.
            Default is ``None``, which will take the value from
            :rc:`legend.edgecolor` If ``"inherit"``, it will take
            :rc:`axes.edgecolor`.

        mode : {"expand", None}
            If `mode` is set to ``"expand"`` the legend will be horizontally
            expanded to fill the axes area (or `bbox_to_anchor` if defines
            the legend's size).

        bbox_transform : None or :class:`matplotlib.transforms.Transform`
            The transform for the bounding box (`bbox_to_anchor`). For a value
            of ``None`` (default) the Axes'
            :data:`~matplotlib.axes.Axes.transAxes` transform will be used.

        title : str or None
            The legend's title. Default is no title (``None``).

        borderpad : float or None
            The fractional whitespace inside the legend border.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.borderpad`.

        labelspacing : float or None
            The vertical space between the legend entries.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.labelspacing`.

        handlelength : float or None
            The length of the legend handles.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.handlelength`.

        handletextpad : float or None
            The pad between the legend handle and text.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.handletextpad`.

        borderaxespad : float or None
            The pad between the axes and legend border.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.borderaxespad`.

        columnspacing : float or None
            The spacing between columns.
            Measured in font-size units.
            Default is ``None``, which will take the value from
            :rc:`legend.columnspacing`.

        handler_map : dict or None
            The custom dictionary mapping instances or types to a legend
            handler. This `handler_map` updates the default handler map
            found at :func:`matplotlib.legend.Legend.get_legend_handler_map`.

        Returns
        -------

        :class:`matplotlib.legend.Legend` instance

        Notes
        -----

        Not all kinds of artist are supported by the legend command. See
        :ref:`sphx_glr_tutorials_intermediate_legend_guide.py` for details.

        Examples
        --------

        .. plot:: gallery/api/legend.py

        """
        handles, labels, extra_args, kwargs = mlegend._parse_legend_args(
                [self],
                *args,
                **kwargs)
        if len(extra_args):
            raise TypeError('legend only accepts two non-keyword arguments')
        self.legend_ = mlegend.Legend(self, handles, labels, **kwargs)
        self.legend_._remove_method = lambda h: setattr(self, 'legend_', None)
        return self.legend_

    def text(self, x, y, s, fontdict=None, withdash=False, **kwargs):
        """
        Add text to the axes.

        Add the text *s* to the axes at location *x*, *y* in data coordinates.

        Parameters
        ----------
        x, y : scalars
            The position to place the text. By default, this is in data
            coordinates. The coordinate system can be changed using the
            *transform* parameter.

        s : str
            The text.

        fontdict : dictionary, optional, default: None
            A dictionary to override the default text properties. If fontdict
            is None, the defaults are determined by your rc parameters.

        withdash : boolean, optional, default: False
            Creates a `~matplotlib.text.TextWithDash` instance instead of a
            `~matplotlib.text.Text` instance.

        Returns
        -------
        text : `.Text`
            The created `.Text` instance.

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.text.Text` properties.
            Other miscellaneous text parameters.

        Examples
        --------
        Individual keyword arguments can be used to override any given
        parameter::

            >>> text(x, y, s, fontsize=12)

        The default transform specifies that text is in data coords,
        alternatively, you can specify text in axis coords (0,0 is
        lower-left and 1,1 is upper-right).  The example below places
        text in the center of the axes::

            >>> text(0.5, 0.5, 'matplotlib', horizontalalignment='center',
            ...      verticalalignment='center', transform=ax.transAxes)

        You can put a rectangular box around the text instance (e.g., to
        set a background color) by using the keyword `bbox`.  `bbox` is
        a dictionary of `~matplotlib.patches.Rectangle`
        properties.  For example::

            >>> text(x, y, s, bbox=dict(facecolor='red', alpha=0.5))
        """
        default = {
            'verticalalignment': 'baseline',
            'horizontalalignment': 'left',
            'transform': self.transData,
            'clip_on': False}

        # At some point if we feel confident that TextWithDash
        # is robust as a drop-in replacement for Text and that
        # the performance impact of the heavier-weight class
        # isn't too significant, it may make sense to eliminate
        # the withdash kwarg and simply delegate whether there's
        # a dash to TextWithDash and dashlength.
        if withdash:
            t = mtext.TextWithDash(
                x=x, y=y, text=s)
        else:
            t = mtext.Text(
                x=x, y=y, text=s)

        t.update(default)
        if fontdict is not None:
            t.update(fontdict)
        t.update(kwargs)

        t.set_clip_path(self.patch)
        self._add_text(t)
        return t

    @docstring.dedent_interpd
    def annotate(self, *args, **kwargs):
        a = mtext.Annotation(*args, **kwargs)
        a.set_transform(mtransforms.IdentityTransform())
        if 'clip_on' in kwargs:
            a.set_clip_path(self.patch)
        self._add_text(a)
        return a
    annotate.__doc__ = mtext.Annotation.__init__.__doc__
    #### Lines and spans

    @docstring.dedent_interpd
    def axhline(self, y=0, xmin=0, xmax=1, **kwargs):
        """
        Add a horizontal line across the axis.

        Parameters
        ----------
        y : scalar, optional, default: 0
            y position in data coordinates of the horizontal line.

        xmin : scalar, optional, default: 0
            Should be between 0 and 1, 0 being the far left of the plot, 1 the
            far right of the plot.

        xmax : scalar, optional, default: 1
            Should be between 0 and 1, 0 being the far left of the plot, 1 the
            far right of the plot.

        Returns
        -------
        :class:`~matplotlib.lines.Line2D`

        Other Parameters
        ----------------
        **kwargs :
            Valid kwargs are :class:`~matplotlib.lines.Line2D` properties,
            with the exception of 'transform':

            %(Line2D)s

        See also
        --------
        hlines : Add horizontal lines in data coordinates.
        axhspan : Add a horizontal span (rectangle) across the axis.

        Examples
        --------

        * draw a thick red hline at 'y' = 0 that spans the xrange::

            >>> axhline(linewidth=4, color='r')

        * draw a default hline at 'y' = 1 that spans the xrange::

            >>> axhline(y=1)

        * draw a default hline at 'y' = .5 that spans the middle half of
          the xrange::

            >>> axhline(y=.5, xmin=0.25, xmax=0.75)

        """
        if "transform" in kwargs:
            raise ValueError(
                "'transform' is not allowed as a kwarg;"
                + "axhline generates its own transform.")
        ymin, ymax = self.get_ybound()

        # We need to strip away the units for comparison with
        # non-unitized bounds
        self._process_unit_info(ydata=y, kwargs=kwargs)
        yy = self.convert_yunits(y)
        scaley = (yy < ymin) or (yy > ymax)

        trans = self.get_yaxis_transform(which='grid')
        l = mlines.Line2D([xmin, xmax], [y, y], transform=trans, **kwargs)
        self.add_line(l)
        self.autoscale_view(scalex=False, scaley=scaley)
        return l

    @docstring.dedent_interpd
    def axvline(self, x=0, ymin=0, ymax=1, **kwargs):
        """
        Add a vertical line across the axes.

        Parameters
        ----------
        x : scalar, optional, default: 0
            x position in data coordinates of the vertical line.

        ymin : scalar, optional, default: 0
            Should be between 0 and 1, 0 being the bottom of the plot, 1 the
            top of the plot.

        ymax : scalar, optional, default: 1
            Should be between 0 and 1, 0 being the bottom of the plot, 1 the
            top of the plot.

        Returns
        -------
        :class:`~matplotlib.lines.Line2D`

        Other Parameters
        ----------------
        **kwargs :
            Valid kwargs are :class:`~matplotlib.lines.Line2D` properties,
            with the exception of 'transform':

            %(Line2D)s

        Examples
        --------
        * draw a thick red vline at *x* = 0 that spans the yrange::

            >>> axvline(linewidth=4, color='r')

        * draw a default vline at *x* = 1 that spans the yrange::

            >>> axvline(x=1)

        * draw a default vline at *x* = .5 that spans the middle half of
          the yrange::

            >>> axvline(x=.5, ymin=0.25, ymax=0.75)

        See also
        --------
        vlines : Add vertical lines in data coordinates.
        axvspan : Add a vertical span (rectangle) across the axis.
        """

        if "transform" in kwargs:
            raise ValueError(
                "'transform' is not allowed as a kwarg;"
                + "axvline generates its own transform.")
        xmin, xmax = self.get_xbound()

        # We need to strip away the units for comparison with
        # non-unitized bounds
        self._process_unit_info(xdata=x, kwargs=kwargs)
        xx = self.convert_xunits(x)
        scalex = (xx < xmin) or (xx > xmax)

        trans = self.get_xaxis_transform(which='grid')
        l = mlines.Line2D([x, x], [ymin, ymax], transform=trans, **kwargs)
        self.add_line(l)
        self.autoscale_view(scalex=scalex, scaley=False)
        return l

    @docstring.dedent_interpd
    def axhspan(self, ymin, ymax, xmin=0, xmax=1, **kwargs):
        """
        Add a horizontal span (rectangle) across the axis.

        Draw a horizontal span (rectangle) from *ymin* to *ymax*.
        With the default values of *xmin* = 0 and *xmax* = 1, this
        always spans the xrange, regardless of the xlim settings, even
        if you change them, e.g., with the :meth:`set_xlim` command.
        That is, the horizontal extent is in axes coords: 0=left,
        0.5=middle, 1.0=right but the *y* location is in data
        coordinates.

        Parameters
        ----------
        ymin : float
               Lower limit of the horizontal span in data units.
        ymax : float
               Upper limit of the horizontal span in data units.
        xmin : float, optional, default: 0
               Lower limit of the vertical span in axes (relative
               0-1) units.
        xmax : float, optional, default: 1
               Upper limit of the vertical span in axes (relative
               0-1) units.

        Returns
        -------
        Polygon : `~matplotlib.patches.Polygon`

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.patches.Polygon` properties.

        %(Polygon)s

        See Also
        --------
        axvspan : Add a vertical span across the axes.
        """
        trans = self.get_yaxis_transform(which='grid')

        # process the unit information
        self._process_unit_info([xmin, xmax], [ymin, ymax], kwargs=kwargs)

        # first we need to strip away the units
        xmin, xmax = self.convert_xunits([xmin, xmax])
        ymin, ymax = self.convert_yunits([ymin, ymax])

        verts = (xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)
        p = mpatches.Polygon(verts, **kwargs)
        p.set_transform(trans)
        self.add_patch(p)
        self.autoscale_view(scalex=False)
        return p

    def axvspan(self, xmin, xmax, ymin=0, ymax=1, **kwargs):
        """
        Add a vertical span (rectangle) across the axes.

        Draw a vertical span (rectangle) from `xmin` to `xmax`.  With
        the default values of `ymin` = 0 and `ymax` = 1. This always
        spans the yrange, regardless of the ylim settings, even if you
        change them, e.g., with the :meth:`set_ylim` command.  That is,
        the vertical extent is in axes coords: 0=bottom, 0.5=middle,
        1.0=top but the y location is in data coordinates.

        Parameters
        ----------
        xmin : scalar
            Number indicating the first X-axis coordinate of the vertical
            span rectangle in data units.
        xmax : scalar
            Number indicating the second X-axis coordinate of the vertical
            span rectangle in data units.
        ymin : scalar, optional
            Number indicating the first Y-axis coordinate of the vertical
            span rectangle in relative Y-axis units (0-1). Default to 0.
        ymax : scalar, optional
            Number indicating the second Y-axis coordinate of the vertical
            span rectangle in relative Y-axis units (0-1). Default to 1.

        Returns
        -------
        rectangle : matplotlib.patches.Polygon
            Vertical span (rectangle) from (xmin, ymin) to (xmax, ymax).

        Other Parameters
        ----------------
        **kwargs
            Optional parameters are properties of the class
            matplotlib.patches.Polygon.

        See Also
        --------
        axhspan : Add a horizontal span across the axes.

        Examples
        --------
        Draw a vertical, green, translucent rectangle from x = 1.25 to
        x = 1.55 that spans the yrange of the axes.

        >>> axvspan(1.25, 1.55, facecolor='g', alpha=0.5)

        """
        trans = self.get_xaxis_transform(which='grid')

        # process the unit information
        self._process_unit_info([xmin, xmax], [ymin, ymax], kwargs=kwargs)

        # first we need to strip away the units
        xmin, xmax = self.convert_xunits([xmin, xmax])
        ymin, ymax = self.convert_yunits([ymin, ymax])

        verts = [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)]
        p = mpatches.Polygon(verts, **kwargs)
        p.set_transform(trans)
        self.add_patch(p)
        self.autoscale_view(scaley=False)
        return p

    @_preprocess_data(replace_names=["y", "xmin", "xmax", "colors"],
                      label_namer="y")
    def hlines(self, y, xmin, xmax, colors='k', linestyles='solid',
               label='', **kwargs):
        """
        Plot horizontal lines at each *y* from *xmin* to *xmax*.

        Parameters
        ----------
        y : scalar or sequence of scalar
            y-indexes where to plot the lines.

        xmin, xmax : scalar or 1D array_like
            Respective beginning and end of each line. If scalars are
            provided, all lines will have same length.

        colors : array_like of colors, optional, default: 'k'

        linestyles : ['solid' | 'dashed' | 'dashdot' | 'dotted'], optional

        label : string, optional, default: ''

        Returns
        -------
        lines : `~matplotlib.collections.LineCollection`

        Other Parameters
        ----------------
        **kwargs :  `~matplotlib.collections.LineCollection` properties.

        See also
        --------
        vlines : vertical lines
        axhline: horizontal line across the axes
        """

        # We do the conversion first since not all unitized data is uniform
        # process the unit information
        self._process_unit_info([xmin, xmax], y, kwargs=kwargs)
        y = self.convert_yunits(y)
        xmin = self.convert_xunits(xmin)
        xmax = self.convert_xunits(xmax)

        if not iterable(y):
            y = [y]
        if not iterable(xmin):
            xmin = [xmin]
        if not iterable(xmax):
            xmax = [xmax]

        y, xmin, xmax = cbook.delete_masked_points(y, xmin, xmax)

        y = np.ravel(y)
        xmin = np.resize(xmin, y.shape)
        xmax = np.resize(xmax, y.shape)

        verts = [((thisxmin, thisy), (thisxmax, thisy))
                 for thisxmin, thisxmax, thisy in zip(xmin, xmax, y)]
        lines = mcoll.LineCollection(verts, colors=colors,
                                     linestyles=linestyles, label=label)
        self.add_collection(lines, autolim=False)
        lines.update(kwargs)

        if len(y) > 0:
            minx = min(xmin.min(), xmax.min())
            maxx = max(xmin.max(), xmax.max())
            miny = y.min()
            maxy = y.max()

            corners = (minx, miny), (maxx, maxy)

            self.update_datalim(corners)
            self.autoscale_view()

        return lines

    @_preprocess_data(replace_names=["x", "ymin", "ymax", "colors"],
                      label_namer="x")
    def vlines(self, x, ymin, ymax, colors='k', linestyles='solid',
               label='', **kwargs):
        """
        Plot vertical lines.

        Plot vertical lines at each *x* from *ymin* to *ymax*.

        Parameters
        ----------
        x : scalar or 1D array_like
            x-indexes where to plot the lines.

        ymin, ymax : scalar or 1D array_like
            Respective beginning and end of each line. If scalars are
            provided, all lines will have same length.

        colors : array_like of colors, optional, default: 'k'

        linestyles : ['solid' | 'dashed' | 'dashdot' | 'dotted'], optional

        label : string, optional, default: ''

        Returns
        -------
        lines : `~matplotlib.collections.LineCollection`

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.collections.LineCollection` properties.

        See also
        --------
        hlines : horizontal lines
        axvline: vertical line across the axes
        """

        self._process_unit_info(xdata=x, ydata=[ymin, ymax], kwargs=kwargs)

        # We do the conversion first since not all unitized data is uniform
        x = self.convert_xunits(x)
        ymin = self.convert_yunits(ymin)
        ymax = self.convert_yunits(ymax)

        if not iterable(x):
            x = [x]
        if not iterable(ymin):
            ymin = [ymin]
        if not iterable(ymax):
            ymax = [ymax]

        x, ymin, ymax = cbook.delete_masked_points(x, ymin, ymax)

        x = np.ravel(x)
        ymin = np.resize(ymin, x.shape)
        ymax = np.resize(ymax, x.shape)

        verts = [((thisx, thisymin), (thisx, thisymax))
                 for thisx, thisymin, thisymax in zip(x, ymin, ymax)]
        lines = mcoll.LineCollection(verts, colors=colors,
                                     linestyles=linestyles, label=label)
        self.add_collection(lines, autolim=False)
        lines.update(kwargs)

        if len(x) > 0:
            minx = x.min()
            maxx = x.max()
            miny = min(ymin.min(), ymax.min())
            maxy = max(ymin.max(), ymax.max())

            corners = (minx, miny), (maxx, maxy)
            self.update_datalim(corners)
            self.autoscale_view()

        return lines

    @_preprocess_data(replace_names=["positions", "lineoffsets",
                                     "linelengths", "linewidths",
                                     "colors", "linestyles"],
                      label_namer=None)
    @docstring.dedent_interpd
    def eventplot(self, positions, orientation='horizontal', lineoffsets=1,
                  linelengths=1, linewidths=None, colors=None,
                  linestyles='solid', **kwargs):
        """
        Plot identical parallel lines at the given positions.

        *positions* should be a 1D or 2D array-like object, with each row
        corresponding to a row or column of lines.

        This type of plot is commonly used in neuroscience for representing
        neural events, where it is usually called a spike raster, dot raster,
        or raster plot.

        However, it is useful in any situation where you wish to show the
        timing or position of multiple sets of discrete events, such as the
        arrival times of people to a business on each day of the month or the
        date of hurricanes each year of the last century.

        Parameters
        ----------
        positions : 1D or 2D array-like object
            Each value is an event. If *positions* is a 2D array-like, each
            row corresponds to a row or a column of lines (depending on the
            *orientation* parameter).

        orientation : {'horizontal', 'vertical'}, optional
            Controls the direction of the event collections:

                - 'horizontal' : the lines are arranged horizontally in rows,
                  and are vertical.
                - 'vertical' : the lines are arranged vertically in columns,
                  and are horizontal.

        lineoffsets : scalar or sequence of scalars, optional, default: 1
            The offset of the center of the lines from the origin, in the
            direction orthogonal to *orientation*.

        linelengths : scalar or sequence of scalars, optional, default: 1
            The total height of the lines (i.e. the lines stretches from
            ``lineoffset - linelength/2`` to ``lineoffset + linelength/2``).

        linewidths : scalar, scalar sequence or None, optional, default: None
            The line width(s) of the event lines, in points. If it is None,
            defaults to its rcParams setting.

        colors : color, sequence of colors or None, optional, default: None
            The color(s) of the event lines. If it is None, defaults to its
            rcParams setting.

        linestyles : str or tuple or a sequence of such values, optional
            Default is 'solid'. Valid strings are ['solid', 'dashed',
            'dashdot', 'dotted', '-', '--', '-.', ':']. Dash tuples
            should be of the form::

                (offset, onoffseq),

            where *onoffseq* is an even length tuple of on and off ink
            in points.

        **kwargs : optional
            Other keyword arguments are line collection properties.  See
            :class:`~matplotlib.collections.LineCollection` for a list of
            the valid properties.

        Returns
        -------

        A list of :class:`matplotlib.collections.EventCollection` objects that
        were added.

        Notes
        -----

        For *linelengths*, *linewidths*, *colors*, and *linestyles*, if only
        a single value is given, that value is applied to all lines.  If an
        array-like is given, it must have the same length as *positions*, and
        each value will be applied to the corresponding row of the array.

        Examples
        --------

        .. plot:: gallery/lines_bars_and_markers/eventplot_demo.py
        """
        self._process_unit_info(xdata=positions,
                                ydata=[lineoffsets, linelengths],
                                kwargs=kwargs)

        # We do the conversion first since not all unitized data is uniform
        positions = self.convert_xunits(positions)
        lineoffsets = self.convert_yunits(lineoffsets)
        linelengths = self.convert_yunits(linelengths)

        if not iterable(positions):
            positions = [positions]
        elif any(iterable(position) for position in positions):
            positions = [np.asanyarray(position) for position in positions]
        else:
            positions = [np.asanyarray(positions)]

        if len(positions) == 0:
            return []

        # prevent 'singular' keys from **kwargs dict from overriding the effect
        # of 'plural' keyword arguments (e.g. 'color' overriding 'colors')
        colors = cbook.local_over_kwdict(colors, kwargs, 'color')
        linewidths = cbook.local_over_kwdict(linewidths, kwargs, 'linewidth')
        linestyles = cbook.local_over_kwdict(linestyles, kwargs, 'linestyle')

        if not iterable(lineoffsets):
            lineoffsets = [lineoffsets]
        if not iterable(linelengths):
            linelengths = [linelengths]
        if not iterable(linewidths):
            linewidths = [linewidths]
        if not iterable(colors):
            colors = [colors]
        if hasattr(linestyles, 'lower') or not iterable(linestyles):
            linestyles = [linestyles]

        lineoffsets = np.asarray(lineoffsets)
        linelengths = np.asarray(linelengths)
        linewidths = np.asarray(linewidths)

        if len(lineoffsets) == 0:
            lineoffsets = [None]
        if len(linelengths) == 0:
            linelengths = [None]
        if len(linewidths) == 0:
            lineoffsets = [None]
        if len(linewidths) == 0:
            lineoffsets = [None]
        if len(colors) == 0:
            colors = [None]
        try:
            # Early conversion of the colors into RGBA values to take care
            # of cases like colors='0.5' or colors='C1'.  (Issue #8193)
            colors = mcolors.to_rgba_array(colors)
        except ValueError:
            # Will fail if any element of *colors* is None. But as long
            # as len(colors) == 1 or len(positions), the rest of the
            # code should process *colors* properly.
            pass

        if len(lineoffsets) == 1 and len(positions) != 1:
            lineoffsets = np.tile(lineoffsets, len(positions))
            lineoffsets[0] = 0
            lineoffsets = np.cumsum(lineoffsets)
        if len(linelengths) == 1:
            linelengths = np.tile(linelengths, len(positions))
        if len(linewidths) == 1:
            linewidths = np.tile(linewidths, len(positions))
        if len(colors) == 1:
            colors = list(colors)
            colors = colors * len(positions)
        if len(linestyles) == 1:
            linestyles = [linestyles] * len(positions)

        if len(lineoffsets) != len(positions):
            raise ValueError('lineoffsets and positions are unequal sized '
                             'sequences')
        if len(linelengths) != len(positions):
            raise ValueError('linelengths and positions are unequal sized '
                             'sequences')
        if len(linewidths) != len(positions):
            raise ValueError('linewidths and positions are unequal sized '
                             'sequences')
        if len(colors) != len(positions):
            raise ValueError('colors and positions are unequal sized '
                             'sequences')
        if len(linestyles) != len(positions):
            raise ValueError('linestyles and positions are unequal sized '
                             'sequences')

        colls = []
        for position, lineoffset, linelength, linewidth, color, linestyle in \
            zip(positions, lineoffsets, linelengths, linewidths,
                           colors, linestyles):
            coll = mcoll.EventCollection(position,
                                         orientation=orientation,
                                         lineoffset=lineoffset,
                                         linelength=linelength,
                                         linewidth=linewidth,
                                         color=color,
                                         linestyle=linestyle)
            self.add_collection(coll, autolim=False)
            coll.update(kwargs)
            colls.append(coll)

        if len(positions) > 0:
            # try to get min/max
            min_max = [(np.min(_p), np.max(_p)) for _p in positions
                       if len(_p) > 0]
            # if we have any non-empty positions, try to autoscale
            if len(min_max) > 0:
                mins, maxes = zip(*min_max)
                minpos = np.min(mins)
                maxpos = np.max(maxes)

                minline = (lineoffsets - linelengths).min()
                maxline = (lineoffsets + linelengths).max()

                if (orientation is not None and
                        orientation.lower() == "vertical"):
                    corners = (minline, minpos), (maxline, maxpos)
                else:  # "horizontal", None or "none" (see EventCollection)
                    corners = (minpos, minline), (maxpos, maxline)
                self.update_datalim(corners)
                self.autoscale_view()

        return colls

    # ### Basic plotting
    # The label_naming happens in `matplotlib.axes._base._plot_args`
    @_preprocess_data(replace_names=["x", "y"],
                      positional_parameter_names=_plot_args_replacer,
                      label_namer=None)
    @docstring.dedent_interpd
    def plot(self, *args, **kwargs):
        """
        Plot y versus x as lines and/or markers.

        Call signatures::

            plot([x], y, [fmt], data=None, **kwargs)
            plot([x], y, [fmt], [x2], y2, [fmt2], ..., **kwargs)

        The coordinates of the points or line nodes are given by *x*, *y*.

        The optional parameter *fmt* is a convenient way for defining basic
        formatting like color, marker and linestyle. It's a shortcut string
        notation described in the *Notes* section below.

        >>> plot(x, y)        # plot x and y using default line style and color
        >>> plot(x, y, 'bo')  # plot x and y using blue circle markers
        >>> plot(y)           # plot y using x as index array 0..N-1
        >>> plot(y, 'r+')     # ditto, but with red plusses

        You can use `.Line2D` properties as keyword arguments for more
        control on the  appearance. Line properties and *fmt* can be mixed.
        The following two calls yield identical results:

        >>> plot(x, y, 'go--', linewidth=2, markersize=12)
        >>> plot(x, y, color='green', marker='o', linestyle='dashed',
                linewidth=2, markersize=12)

        When conflicting with *fmt*, keyword arguments take precedence.

        **Plotting labelled data**

        There's a convenient way for plotting objects with labelled data (i.e.
        data that can be accessed by index ``obj['y']``). Instead of giving
        the data in *x* and *y*, you can provide the object in the *data*
        parameter and just give the labels for *x* and *y*::

        >>> plot('xlabel', 'ylabel', data=obj)

        All indexable objects are supported. This could e.g. be a `dict`, a
        `pandas.DataFame` or a structured numpy array.


        **Plotting multiple sets of data**

        There are various ways to plot multiple sets of data.

        - The most straight forward way is just to call `plot` multiple times.
          Example:

          >>> plot(x1, y1, 'bo')
          >>> plot(x2, y2, 'go')

        - Alternatively, if your data is already a 2d array, you can pass it
          directly to *x*, *y*. A separate data set will be drawn for every
          column.

          Example: an array ``a`` where the first column represents the *x*
          values and the other columns are the *y* columns::

          >>> plot(a[0], a[1:])

        - The third way is to specify multiple sets of *[x]*, *y*, *[fmt]*
          groups::

          >>> plot(x1, y1, 'g^', x2, y2, 'g-')

          In this case, any additional keyword argument applies to all
          datasets. Also this syntax cannot be combined with the *data*
          parameter.

        By default, each line is assigned a different style specified by a
        'style cycle'. The *fmt* and line property parameters are only
        necessary if you want explicit deviations from these defaults.
        Alternatively, you can also change the style cycle using the
        'axes.prop_cycle' rcParam.

        Parameters
        ----------
        x, y : array-like or scalar
            The horizontal / vertical coordinates of the data points.
            *x* values are optional. If not given, they default to
            ``[0, ..., N-1]``.

            Commonly, these parameters are arrays of length N. However,
            scalars are supported as well (equivalent to an array with
            constant value).

            The parameters can also be 2-dimensional. Then, the columns
            represent separate data sets.

        fmt : str, optional
            A format string, e.g. 'ro' for red circles. See the *Notes*
            section for a full description of the format strings.

            Format strings are just an abbreviation for quickly setting
            basic line properties. All of these and more can also be
            controlled by keyword arguments.

        data : indexable object, optional
            An object with labelled data. If given, provide the label names to
            plot in *x* and *y*.

            .. note::
                Technically there's a slight ambiguity in calls where the
                second label is a valid *fmt*. `plot('n', 'o', data=obj)`
                could be `plt(x, y)` or `plt(y, fmt)`. In such cases,
                the former interpretation is chosen, but a warning is issued.
                You may suppress the warning by adding an empty format string
                `plot('n', 'o', '', data=obj)`.


        Other Parameters
        ----------------
        scalex, scaley : bool, optional, default: True
            These parameters determined if the view limits are adapted to
            the data limits. The values are passed on to `autoscale_view`.

        **kwargs : `.Line2D` properties, optional
            *kwargs* are used to specify properties like a line label (for
            auto legends), linewidth, antialiasing, marker face color.
            Example::

            >>> plot([1,2,3], [1,2,3], 'go-', label='line 1', linewidth=2)
            >>> plot([1,2,3], [1,4,9], 'rs',  label='line 2')

            If you make multiple lines with one plot command, the kwargs
            apply to all those lines.

            Here is a list of available `.Line2D` properties:

            %(Line2D)s

        Returns
        -------
        lines
            A list of `.Line2D` objects representing the plotted data.


        See Also
        --------
        scatter : XY scatter plot with markers of variing size and/or color (
            sometimes also called bubble chart).


        Notes
        -----
        **Format Strings**

        A format string consists of a part for color, marker and line::

            fmt = '[color][marker][line]'

        Each of them is optional. If not provided, the value from the style
        cycle is used. Exception: If ``line`` is given, but no ``marker``,
        the data will be a line without markers.

        **Colors**

        The following color abbreviations are supported:

        =============    ===============================
        character        color
        =============    ===============================
        ``'b'``          blue
        ``'g'``          green
        ``'r'``          red
        ``'c'``          cyan
        ``'m'``          magenta
        ``'y'``          yellow
        ``'k'``          black
        ``'w'``          white
        =============    ===============================

        If the color is the only part of the format string, you can
        additionally use any  `matplotlib.colors` spec, e.g. full names
        (``'green'``) or hex strings (``'#008000'``).

        **Markers**

        =============    ===============================
        character        description
        =============    ===============================
        ``'.'``          point marker
        ``','``          pixel marker
        ``'o'``          circle marker
        ``'v'``          triangle_down marker
        ``'^'``          triangle_up marker
        ``'<'``          triangle_left marker
        ``'>'``          triangle_right marker
        ``'1'``          tri_down marker
        ``'2'``          tri_up marker
        ``'3'``          tri_left marker
        ``'4'``          tri_right marker
        ``'s'``          square marker
        ``'p'``          pentagon marker
        ``'*'``          star marker
        ``'h'``          hexagon1 marker
        ``'H'``          hexagon2 marker
        ``'+'``          plus marker
        ``'x'``          x marker
        ``'D'``          diamond marker
        ``'d'``          thin_diamond marker
        ``'|'``          vline marker
        ``'_'``          hline marker
        =============    ===============================

        **Line Styles**

        =============    ===============================
        character        description
        =============    ===============================
        ``'-'``          solid line style
        ``'--'``         dashed line style
        ``'-.'``         dash-dot line style
        ``':'``          dotted line style
        =============    ===============================

        Example format strings::

            'b'    # blue markers with default shape
            'ro'   # red circles
            'g-'   # green solid line
            '--'   # dashed line with default color
            'k^:'  # black triangle_up markers connected by a dotted line

        """
        scalex = kwargs.pop('scalex', True)
        scaley = kwargs.pop('scaley', True)

        if not self._hold:
            self.cla()
        lines = []

        kwargs = cbook.normalize_kwargs(kwargs, _alias_map)

        for line in self._get_lines(*args, **kwargs):
            self.add_line(line)
            lines.append(line)

        self.autoscale_view(scalex=scalex, scaley=scaley)
        return lines

    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    @docstring.dedent_interpd
    def plot_date(self, x, y, fmt='o', tz=None, xdate=True, ydate=False,
                  **kwargs):
        """
        Plot data that contains dates.

        Similar to `.plot`, this plots *y* vs. *x* as lines or markers.
        However, the axis labels are formatted as dates depending on *xdate*
        and *ydate*.

        Parameters
        ----------
        x, y : array-like
            The coordinates of the data points. If *xdate* or *ydate* is
            *True*, the respective values *x* or *y* are interpreted as
            :ref:`Matplotlib dates <date-format>`.

        fmt : str, optional
            The plot format string. For details, see the corresponding
            parameter in `.plot`.

        tz : [ *None* | timezone string | :class:`tzinfo` instance]
            The time zone to use in labeling dates. If *None*, defaults to
            rcParam ``timezone``.

        xdate : bool, optional, default: True
            If *True*, the *x*-axis will be interpreted as Matplotlib dates.

        ydate : bool, optional, default: False
            If *True*, the *y*-axis will be interpreted as Matplotlib dates.


        Returns
        -------
        lines
            A list of `~.Line2D` objects representing the plotted data.


        Other Parameters
        ----------------
        **kwargs
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s


        See Also
        --------
        matplotlib.dates : Helper functions on dates.
        matplotlib.dates.date2num : Convert dates to num.
        matplotlib.dates.num2date : Convert num to dates.
        matplotlib.dates.drange : Create an equally spaced sequence of dates.


        Notes
        -----
        If you are using custom date tickers and formatters, it may be
        necessary to set the formatters/locators after the call to
        `.plot_date`. `.plot_date` will set the default tick locator to
        `.AutoDateLocator` (if the tick locator is not already set to a
        `.DateLocator` instance) and the default tick formatter to
        `.AutoDateFormatter` (if the tick formatter is not already set to a
        `.DateFormatter` instance).
        """

        if not self._hold:
            self.cla()

        if xdate:
            self.xaxis_date(tz)
        if ydate:
            self.yaxis_date(tz)

        ret = self.plot(x, y, fmt, **kwargs)

        self.autoscale_view()

        return ret

    # @_preprocess_data() # let 'plot' do the unpacking..
    @docstring.dedent_interpd
    def loglog(self, *args, **kwargs):
        """
        Make a plot with log scaling on both the x and y axis.

        Call signatures::

            loglog([x], y, [fmt], data=None, **kwargs)
            loglog([x], y, [fmt], [x2], y2, [fmt2], ..., **kwargs)

        This is just a thin wrapper around `.plot` which additionally changes
        both the x-axis and the y-axis to log scaling. All of the concepts and
        parameters of plot can be used here as well.

        The additional parameters *basex/y*, *subsx/y* and *nonposx/y* control
        the x/y-axis properties. They are just forwarded to `.Axes.set_xscale`
        and `.Axes.set_yscale`.

        Parameters
        ----------
        basex, basey : scalar, optional, default 10
            Base of the x/y logarithm.

        subsx, subsy : sequence, optional
            The location of the minor x/y ticks. If *None*, reasonable
            locations are automatically chosen depending on the number of
            decades in the plot.
            See `.Axes.set_xscale` / `.Axes.set_yscale` for details.

        nonposx, nonposy : {'mask', 'clip'}, optional, default 'mask'
            Non-positive values in x or y can be masked as invalid, or clipped
            to a very small positive number.

        Returns
        -------
        lines
            A list of `~.Line2D` objects representing the plotted data.

        Other Parameters
        ----------------
        **kwargs
            All parameters supported by `.plot`.
        """
        if not self._hold:
            self.cla()

        dx = {k: kwargs.pop(k) for k in ['basex', 'subsx', 'nonposx']
                if k in kwargs}
        dy = {k: kwargs.pop(k) for k in ['basey', 'subsy', 'nonposy']
                if k in kwargs}

        self.set_xscale('log', **dx)
        self.set_yscale('log', **dy)

        b = self._hold
        self._hold = True  # we've already processed the hold
        l = self.plot(*args, **kwargs)
        self._hold = b  # restore the hold

        return l

    # @_preprocess_data() # let 'plot' do the unpacking..
    @docstring.dedent_interpd
    def semilogx(self, *args, **kwargs):
        """
        Make a plot with log scaling on the x axis.

        Call signatures::

            semilogx([x], y, [fmt], data=None, **kwargs)
            semilogx([x], y, [fmt], [x2], y2, [fmt2], ..., **kwargs)

        This is just a thin wrapper around `.plot` which additionally changes
        the x-axis to log scaling. All of the concepts and parameters of plot
        can be used here as well.

        The additional parameters *basex*, *subsx* and *nonposx* control the
        x-axis properties. They are just forwarded to `.Axes.set_xscale`.

        Parameters
        ----------
        basex : scalar, optional, default 10
            Base of the x logarithm.

        subsx : array_like, optional
            The location of the minor xticks. If *None*, reasonable locations
            are automatically chosen depending on the number of decades in the
            plot. See `.Axes.set_xscale` for details.

        nonposx : {'mask', 'clip'}, optional, default 'mask'
            Non-positive values in x can be masked as invalid, or clipped to a
            very small positive number.

        Returns
        -------
        lines
            A list of `~.Line2D` objects representing the plotted data.

        Other Parameters
        ----------------
        **kwargs
            All parameters supported by `.plot`.
        """
        if not self._hold:
            self.cla()
        d = {k: kwargs.pop(k) for k in ['basex', 'subsx', 'nonposx']
                if k in kwargs}

        self.set_xscale('log', **d)
        b = self._hold
        self._hold = True  # we've already processed the hold
        l = self.plot(*args, **kwargs)
        self._hold = b  # restore the hold
        return l

    # @_preprocess_data() # let 'plot' do the unpacking..
    @docstring.dedent_interpd
    def semilogy(self, *args, **kwargs):
        """
        Make a plot with log scaling on the y axis.

        Call signatures::

            semilogy([x], y, [fmt], data=None, **kwargs)
            semilogy([x], y, [fmt], [x2], y2, [fmt2], ..., **kwargs)

        This is just a thin wrapper around `.plot` which additionally changes
        the y-axis to log scaling. All of the concepts and parameters of plot
        can be used here as well.

        The additional parameters *basey*, *subsy* and *nonposy* control the
        y-axis properties. They are just forwarded to `.Axes.set_yscale`.

        Parameters
        ----------
        basey : scalar, optional, default 10
            Base of the y logarithm.

        subsy : array_like, optional
            The location of the minor yticks. If *None*, reasonable locations
            are automatically chosen depending on the number of decades in the
            plot. See `.Axes.set_yscale` for details.

        nonposy : {'mask', 'clip'}, optional, default 'mask'
            Non-positive values in y can be masked as invalid, or clipped to a
            very small positive number.

        Returns
        -------
        lines
            A list of `~.Line2D` objects representing the plotted data.

        Other Parameters
        ----------------
        **kwargs
            All parameters supported by `.plot`.
        """
        if not self._hold:
            self.cla()
        d = {k: kwargs.pop(k) for k in ['basey', 'subsy', 'nonposy']
                if k in kwargs}
        self.set_yscale('log', **d)
        b = self._hold
        self._hold = True  # we've already processed the hold
        l = self.plot(*args, **kwargs)
        self._hold = b  # restore the hold

        return l

    @_preprocess_data(replace_names=["x"], label_namer="x")
    def acorr(self, x, **kwargs):
        """
        Plot the autocorrelation of *x*.

        Parameters
        ----------

        x : sequence of scalar

        hold : bool, optional, *deprecated*, default: True

        detrend : callable, optional, default: `mlab.detrend_none`
            *x* is detrended by the *detrend* callable. Default is no
            normalization.

        normed : bool, optional, default: True
            If ``True``, input vectors are normalised to unit length.

        usevlines : bool, optional, default: True
            If ``True``, `Axes.vlines` is used to plot the vertical lines from
            the origin to the acorr. Otherwise, `Axes.plot` is used.

        maxlags : integer, optional, default: 10
            Number of lags to show. If ``None``, will return all
            ``2 * len(x) - 1`` lags.

        Returns
        -------
        lags : array (lenth ``2*maxlags+1``)
            lag vector.
        c : array  (length ``2*maxlags+1``)
            auto correlation vector.
        line : `.LineCollection` or `.Line2D`
            `.Artist` added to the axes of the correlation.

             `.LineCollection` if *usevlines* is True
             `.Line2D` if *usevlines* is False
        b : `.Line2D` or None
            Horizontal line at 0 if *usevlines* is True
            None *usevlines* is False

        Other Parameters
        ----------------
        linestyle : `~matplotlib.lines.Line2D` prop, optional, default: None
            Only used if usevlines is ``False``.

        marker : string, optional, default: 'o'

        Notes
        -----
        The cross correlation is performed with :func:`numpy.correlate` with
        ``mode = 2``.
        """
        if "hold" in kwargs:
            warnings.warn("the 'hold' kwarg is deprecated", mplDeprecation)
        return self.xcorr(x, x, **kwargs)

    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    def xcorr(self, x, y, normed=True, detrend=mlab.detrend_none,
              usevlines=True, maxlags=10, **kwargs):
        """
        Plot the cross correlation between *x* and *y*.

        The correlation with lag k is defined as sum_n x[n+k] * conj(y[n]).

        Parameters
        ----------
        x : sequence of scalars of length n

        y : sequence of scalars of length n

        hold : bool, optional, *deprecated*, default: True

        detrend : callable, optional, default: `mlab.detrend_none`
            *x* is detrended by the *detrend* callable. Default is no
            normalization.

        normed : bool, optional, default: True
            If ``True``, input vectors are normalised to unit length.

        usevlines : bool, optional, default: True
            If ``True``, `Axes.vlines` is used to plot the vertical lines from
            the origin to the acorr. Otherwise, `Axes.plot` is used.

        maxlags : int, optional
            Number of lags to show. If None, will return all ``2 * len(x) - 1``
            lags. Default is 10.

        Returns
        -------
        lags : array (lenth ``2*maxlags+1``)
            lag vector.
        c : array  (length ``2*maxlags+1``)
            auto correlation vector.
        line : `.LineCollection` or `.Line2D`
            `.Artist` added to the axes of the correlation

             `.LineCollection` if *usevlines* is True
             `.Line2D` if *usevlines* is False
        b : `.Line2D` or None
            Horizontal line at 0 if *usevlines* is True
            None *usevlines* is False

        Other Parameters
        ----------------
        linestyle : `~matplotlib.lines.Line2D` property, optional
            Only used if usevlines is ``False``.

        marker : string, optional
            Default is 'o'.

        Notes
        -----
        The cross correlation is performed with :func:`numpy.correlate` with
        ``mode = 2``.
        """
        if "hold" in kwargs:
            warnings.warn("the 'hold' kwarg is deprecated", mplDeprecation)

        Nx = len(x)
        if Nx != len(y):
            raise ValueError('x and y must be equal length')

        x = detrend(np.asarray(x))
        y = detrend(np.asarray(y))

        correls = np.correlate(x, y, mode=2)

        if normed:
            correls /= np.sqrt(np.dot(x, x) * np.dot(y, y))

        if maxlags is None:
            maxlags = Nx - 1

        if maxlags >= Nx or maxlags < 1:
            raise ValueError('maxlags must be None or strictly '
                             'positive < %d' % Nx)

        lags = np.arange(-maxlags, maxlags + 1)
        correls = correls[Nx - 1 - maxlags:Nx + maxlags]

        if usevlines:
            a = self.vlines(lags, [0], correls, **kwargs)
            # Make label empty so only vertical lines get a legend entry
            kwargs.pop('label', '')
            b = self.axhline(**kwargs)
        else:
            kwargs.setdefault('marker', 'o')
            kwargs.setdefault('linestyle', 'None')
            a, = self.plot(lags, correls, **kwargs)
            b = None
        return lags, correls, a, b

    #### Specialized plotting

    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    def step(self, x, y, *args, **kwargs):
        """
        Make a step plot.

        Call signatures::

            step(x, y, [fmt], *, data=None, where='pre', **kwargs)
            step(x, y, [fmt], x2, y2, [fmt2], ..., *, where='pre', **kwargs)

        This is just a thin wrapper around `.plot` which changes some
        formatting options. Most of the concepts and parameters of plot can be
        used here as well.

        Parameters
        ----------
        x : array_like
            1-D sequence of x positions. It is assumed, but not checked, that
            it is uniformly increasing.

        y : array_like
            1-D sequence of y levels.

        fmt : str, optional
            A format string, e.g. 'g' for a green line. See `.plot` for a more
            detailed description.

            Note: While full format strings are accepted, it is recommended to
            only specify the color. Line styles are currently ignored (use
            the keyword argument *linestyle* instead). Markers are accepted
            and plotted on the given positions, however, this is a rarely
            needed feature for step plots.

        data : indexable object, optional
            An object with labelled data. If given, provide the label names to
            plot in *x* and *y*.

        where : {'pre', 'post', 'mid'}, optional, default 'pre'
            Define where the steps should be placed:

            - 'pre': The y value is continued constantly to the left from
              every *x* position, i.e. the interval ``(x[i-1], x[i]]`` has the
              value ``y[i]``.
            - 'post': The y value is continued constantly to the right from
              every *x* position, i.e. the interval ``[x[i], x[i+1])`` has the
              value ``y[i]``.
            - 'mid': Steps occur half-way between the *x* positions.

        Returns
        -------
        lines
            A list of `.Line2D` objects representing the plotted data.

        Other Parameters
        ----------------
        **kwargs
            Additional parameters are the same as those for `.plot`.

        Notes
        -----
        .. [notes section required to get data note injection right]
        """
        where = kwargs.pop('where', 'pre')
        if where not in ('pre', 'post', 'mid'):
            raise ValueError("'where' argument to step must be "
                             "'pre', 'post' or 'mid'")
        usr_linestyle = kwargs.pop('linestyle', '')
        kwargs['linestyle'] = 'steps-' + where + usr_linestyle

        return self.plot(x, y, *args, **kwargs)

    @_preprocess_data(replace_names=["x", "left",
                                     "height", "width",
                                     "y", "bottom",
                                     "color", "edgecolor", "linewidth",
                                     "tick_label", "xerr", "yerr",
                                     "ecolor"],
                      label_namer=None,
                      replace_all_args=True
                      )
    @docstring.dedent_interpd
    def bar(self, *args, **kwargs):
        r"""
        Make a bar plot.

        Call signatures::

           bar(x, height, *, align='center', **kwargs)
           bar(x, height, width, *, align='center', **kwargs)
           bar(x, height, width, bottom, *, align='center', **kwargs)

        The bars are positioned at *x* with the given *align* ment. Their
        dimensions are given by *width* and *height*. The vertical baseline
        is *bottom* (default 0).

        Each of *x*, *height*, *width*, and *bottom* may either be a scalar
        applying to all bars, or it may be a sequence of length N providing a
        separate value for each bar.


        Parameters
        ----------
        x : sequence of scalars
            The x coordinates of the bars. See also *align* for the
            alignment of the bars to the coordinates.

        height : scalar or sequence of scalars
            The height(s) of the bars.

        width : scalar or array-like, optional
            The width(s) of the bars (default: 0.8).

        bottom : scalar or array-like, optional
            The y coordinate(s) of the bars bases (default: 0).

        align : {'center', 'edge'}, optional, default: 'center'
            Alignment of the bars to the *x* coordinates:

            - 'center': Center the base on the *x* positions.
            - 'edge': Align the left edges of the bars with the *x* positions.

            To align the bars on the right edge pass a negative *width* and
            ``align='edge'``.

        Returns
        -------
        `.BarContainer`
            Container with all the bars and optionally errorbars.

        Other Parameters
        ----------------
        color : scalar or array-like, optional
            The colors of the bar faces.

        edgecolor : scalar or array-like, optional
            The colors of the bar edges.

        linewidth : scalar or array-like, optional
            Width of the bar edge(s). If 0, don't draw edges.

        tick_label : string or array-like, optional
            The tick labels of the bars.
            Default: None (Use default numeric labels.)

        xerr, yerr : scalar or array-like of shape(N,) or shape(2,N), optional
            If not *None*, add horizontal / vertical errorbars to the bar tips.
            The values are +/- sizes relative to the data:

            - scalar: symmetric +/- values for all bars
            - shape(N,): symmetric +/- values for each bar
            - shape(2,N): separate + and - values for each bar

            Default: None

        ecolor : scalar or array-like, optional, default: 'black'
            The line color of the errorbars.

        capsize : scalar, optional
           The length of the error bar caps in points.
           Default: None, which will take the value from
           :rc:`errorbar.capsize`.

        error_kw : dict, optional
            Dictionary of kwargs to be passed to the `~.Axes.errorbar`
            method. Values of *ecolor* or *capsize* defined here take
            precedence over the independent kwargs.

        log : bool, optional, default: False
            If *True*, set the y-axis to be log scale.

        orientation : {'vertical',  'horizontal'}, optional
            *This is for internal use only.* Please use `barh` for
            horizontal bar plots. Default: 'vertical'.

        See also
        --------
        barh: Plot a horizontal bar plot.

        Notes
        -----
        The optional arguments *color*, *edgecolor*, *linewidth*,
        *xerr*, and *yerr* can be either scalars or sequences of
        length equal to the number of bars.  This enables you to use
        bar as the basis for stacked bar charts, or candlestick plots.
        Detail: *xerr* and *yerr* are passed directly to
        :meth:`errorbar`, so they can also have shape 2xN for
        independent specification of lower and upper errors.

        Other optional kwargs:

        %(Rectangle)s

        """
        kwargs = cbook.normalize_kwargs(kwargs, mpatches._patch_alias_map)
        # this is using the lambdas to do the arg/kwarg unpacking rather
        # than trying to re-implement all of that logic our selves.
        matchers = [
            (lambda x, height, width=0.8, bottom=None, **kwargs:
             (False, x, height, width, bottom, kwargs)),
            (lambda left, height, width=0.8, bottom=None, **kwargs:
             (True, left, height, width, bottom, kwargs)),
        ]
        exps = []
        for matcher in matchers:
            try:
                dp, x, height, width, y, kwargs = matcher(*args, **kwargs)
            except TypeError as e:
                # This can only come from a no-match as there is
                # no other logic in the matchers.
                exps.append(e)
            else:
                break
        else:
            raise exps[0]
        # if we matched the second-case, then the user passed in
        # left=val as a kwarg which we want to deprecate
        if dp:
            warnings.warn(
                "The *left* kwarg to `bar` is deprecated use *x* instead. "
                "Support for *left* will be removed in Matplotlib 3.0",
                mplDeprecation, stacklevel=2)
        if not self._hold:
            self.cla()
        color = kwargs.pop('color', None)
        if color is None:
            color = self._get_patches_for_fill.get_next_color()
        edgecolor = kwargs.pop('edgecolor', None)
        linewidth = kwargs.pop('linewidth', None)

        # Because xerr and yerr will be passed to errorbar,
        # most dimension checking and processing will be left
        # to the errorbar method.
        xerr = kwargs.pop('xerr', None)
        yerr = kwargs.pop('yerr', None)
        error_kw = kwargs.pop('error_kw', dict())
        ecolor = kwargs.pop('ecolor', 'k')
        capsize = kwargs.pop('capsize', rcParams["errorbar.capsize"])
        error_kw.setdefault('ecolor', ecolor)
        error_kw.setdefault('capsize', capsize)

        if rcParams['_internal.classic_mode']:
            align = kwargs.pop('align', 'edge')
        else:
            align = kwargs.pop('align', 'center')

        orientation = kwargs.pop('orientation', 'vertical')
        log = kwargs.pop('log', False)
        label = kwargs.pop('label', '')
        tick_labels = kwargs.pop('tick_label', None)

        adjust_ylim = False
        adjust_xlim = False

        if orientation == 'vertical':
            if y is None:
                if self.get_yscale() == 'log':
                    adjust_ylim = True
                y = 0

        elif orientation == 'horizontal':
            if x is None:
                if self.get_xscale() == 'log':
                    adjust_xlim = True
                x = 0

        if orientation == 'vertical':
            self._process_unit_info(xdata=x, ydata=height, kwargs=kwargs)
            if log:
                self.set_yscale('log', nonposy='clip')
        elif orientation == 'horizontal':
            self._process_unit_info(xdata=width, ydata=y, kwargs=kwargs)
            if log:
                self.set_xscale('log', nonposx='clip')
        else:
            raise ValueError('invalid orientation: %s' % orientation)

        # lets do some conversions now since some types cannot be
        # subtracted uniformly
        if self.xaxis is not None:
            x = self.convert_xunits(x)
            width = self.convert_xunits(width)
            if xerr is not None:
                xerr = self.convert_xunits(xerr)

        if self.yaxis is not None:
            y = self.convert_yunits(y)
            height = self.convert_yunits(height)
            if yerr is not None:
                yerr = self.convert_yunits(yerr)

        x, height, width, y, linewidth = np.broadcast_arrays(
            # Make args iterable too.
            np.atleast_1d(x), height, width, y, linewidth)

        # Now that units have been converted, set the tick locations.
        if orientation == 'vertical':
            tick_label_axis = self.xaxis
            tick_label_position = x
        elif orientation == 'horizontal':
            tick_label_axis = self.yaxis
            tick_label_position = y

        linewidth = itertools.cycle(np.atleast_1d(linewidth))
        color = itertools.chain(itertools.cycle(mcolors.to_rgba_array(color)),
                                # Fallback if color == "none".
                                itertools.repeat([0, 0, 0, 0]))
        if edgecolor is None:
            edgecolor = itertools.repeat(None)
        else:
            edgecolor = itertools.chain(
                itertools.cycle(mcolors.to_rgba_array(edgecolor)),
                # Fallback if edgecolor == "none".
                itertools.repeat([0, 0, 0, 0]))

        # We will now resolve the alignment and really have
        # left, bottom, width, height vectors
        if align == 'center':
            if orientation == 'vertical':
                left = x - width / 2
                bottom = y
            elif orientation == 'horizontal':
                bottom = y - height / 2
                left = x
        elif align == 'edge':
            left = x
            bottom = y
        else:
            raise ValueError('invalid alignment: %s' % align)

        patches = []
        args = zip(left, bottom, width, height, color, edgecolor, linewidth)
        for l, b, w, h, c, e, lw in args:
            r = mpatches.Rectangle(
                xy=(l, b), width=w, height=h,
                facecolor=c,
                edgecolor=e,
                linewidth=lw,
                label='_nolegend_',
                )
            r.update(kwargs)
            r.get_path()._interpolation_steps = 100
            if orientation == 'vertical':
                r.sticky_edges.y.append(b)
            elif orientation == 'horizontal':
                r.sticky_edges.x.append(l)
            self.add_patch(r)
            patches.append(r)

        holdstate = self._hold
        self._hold = True  # ensure hold is on before plotting errorbars

        if xerr is not None or yerr is not None:
            if orientation == 'vertical':
                # using list comps rather than arrays to preserve unit info
                ex = [l + 0.5 * w for l, w in zip(left, width)]
                ey = [b + h for b, h in zip(bottom, height)]

            elif orientation == 'horizontal':
                # using list comps rather than arrays to preserve unit info
                ex = [l + w for l, w in zip(left, width)]
                ey = [b + 0.5 * h for b, h in zip(bottom, height)]

            error_kw.setdefault("label", '_nolegend_')

            errorbar = self.errorbar(ex, ey,
                                     yerr=yerr, xerr=xerr,
                                     fmt='none', **error_kw)
        else:
            errorbar = None

        self._hold = holdstate  # restore previous hold state

        if adjust_xlim:
            xmin, xmax = self.dataLim.intervalx
            xmin = min(w for w in width if w > 0)
            if xerr is not None:
                xmin = xmin - np.max(xerr)
            xmin = max(xmin * 0.9, 1e-100)
            self.dataLim.intervalx = (xmin, xmax)

        if adjust_ylim:
            ymin, ymax = self.dataLim.intervaly
            ymin = min(h for h in height if h > 0)
            if yerr is not None:
                ymin = ymin - np.max(yerr)
            ymin = max(ymin * 0.9, 1e-100)
            self.dataLim.intervaly = (ymin, ymax)
        self.autoscale_view()

        bar_container = BarContainer(patches, errorbar, label=label)
        self.add_container(bar_container)

        if tick_labels is not None:
            tick_labels = _backports.broadcast_to(tick_labels, len(patches))
            tick_label_axis.set_ticks(tick_label_position)
            tick_label_axis.set_ticklabels(tick_labels)

        return bar_container

    @docstring.dedent_interpd
    def barh(self, *args, **kwargs):
        r"""
        Make a horizontal bar plot.

        Call signatures::

           bar(y, width, *, align='center', **kwargs)
           bar(y, width, height, *, align='center', **kwargs)
           bar(y, width, height, left, *, align='center', **kwargs)

        The bars are positioned at *y* with the given *align*. Their
        dimensions are given by *width* and *height*. The horizontal baseline
        is *left* (default 0).

        Each of *y*, *width*, *height*, and *left* may either be a scalar
        applying to all bars, or it may be a sequence of length N providing a
        separate value for each bar.


        Parameters
        ----------
        y : scalar or array-like
            The y coordinates of the bars. See also *align* for the
            alignment of the bars to the coordinates.

        width : scalar or array-like
            The width(s) of the bars.

        height : sequence of scalars, optional, default: 0.8
            The heights of the bars.

        left : sequence of scalars
            The x coordinates of the left sides of the bars (default: 0).

        align : {'center', 'edge'}, optional, default: 'center'
            Alignment of the base to the *y* coordinates*:

            - 'center': Center the bars on the *y* positions.
            - 'edge': Align the bottom edges of the bars with the *y*
              positions.

            To align the bars on the top edge pass a negative *height* and
            ``align='edge'``.

        Returns
        -------
        `.BarContainer`
            Container with all the bars and optionally errorbars.

        Other Parameters
        ----------------
        color : scalar or array-like, optional
            The colors of the bar faces.

        edgecolor : scalar or array-like, optional
            The colors of the bar edges.

        linewidth : scalar or array-like, optional
            Width of the bar edge(s). If 0, don't draw edges.

        tick_label : string or array-like, optional
            The tick labels of the bars.
            Default: None (Use default numeric labels.)

        xerr, yerr : scalar or array-like of shape(N,) or shape(2,N), optional
            If not ``None``, add horizontal / vertical errorbars to the
            bar tips. The values are +/- sizes relative to the data:

            - scalar: symmetric +/- values for all bars
            - shape(N,): symmetric +/- values for each bar
            - shape(2,N): separate + and - values for each bar

            Default: None

        ecolor : scalar or array-like, optional, default: 'black'
            The line color of the errorbars.

        capsize : scalar, optional
           The length of the error bar caps in points.
           Default: None, which will take the value from
           :rc:`errorbar.capsize`.

        error_kw : dict, optional
            Dictionary of kwargs to be passed to the `~.Axes.errorbar`
            method. Values of *ecolor* or *capsize* defined here take
            precedence over the independent kwargs.

        log : bool, optional, default: False
            If ``True``, set the x-axis to be log scale.

        See also
        --------
        bar: Plot a vertical bar plot.

        Notes
        -----
        The optional arguments *color*, *edgecolor*, *linewidth*,
        *xerr*, and *yerr* can be either scalars or sequences of
        length equal to the number of bars.  This enables you to use
        bar as the basis for stacked bar charts, or candlestick plots.
        Detail: *xerr* and *yerr* are passed directly to
        :meth:`errorbar`, so they can also have shape 2xN for
        independent specification of lower and upper errors.

        Other optional kwargs:

        %(Rectangle)s

        """
        # this is using the lambdas to do the arg/kwarg unpacking rather
        # than trying to re-implement all of that logic our selves.
        matchers = [
            (lambda y, width, height=0.8, left=None, **kwargs:
             (False, y, width, height, left, kwargs)),
            (lambda bottom, width, height=0.8, left=None, **kwargs:
             (True, bottom, width, height, left, kwargs)),
        ]
        excs = []
        for matcher in matchers:
            try:
                dp, y, width, height, left, kwargs = matcher(*args, **kwargs)
            except TypeError as e:
                # This can only come from a no-match as there is
                # no other logic in the matchers.
                excs.append(e)
            else:
                break
        else:
            raise excs[0]

        if dp:
            warnings.warn(
                "The *bottom* kwarg to `barh` is deprecated use *y* instead. "
                "Support for *bottom* will be removed in Matplotlib 3.0",
                mplDeprecation, stacklevel=2)
        kwargs.setdefault('orientation', 'horizontal')
        patches = self.bar(x=left, height=height, width=width,
                           bottom=y, **kwargs)
        return patches

    @_preprocess_data(label_namer=None)
    @docstring.dedent_interpd
    def broken_barh(self, xranges, yrange, **kwargs):
        """
        Plot a horizontal sequence of rectangles.

        A rectangle is drawn for each element of *xranges*. All rectangles
        have the same vertical position and size defined by *yrange*.

        This is a convenience function for instantiating a
        `.BrokenBarHCollection`, adding it to the axes and autoscaling the
        view.

        Parameters
        ----------
        xranges : sequence of tuples (*xmin*, *xwidth*)
            The x-positions and extends of the rectangles. For each tuple
            (*xmin*, *xwidth*) a rectangle is drawn from *xmin* to *xmin* +
            *xwidth*.
        yranges : (*ymin*, *ymax*)
            The y-position and extend for all the rectangles.

        Other Parameters
        ----------------
        **kwargs : :class:`.BrokenBarHCollection` properties

            Each *kwarg* can be either a single argument applying to all
            rectangles, e.g.::

                facecolors='black'

            or a sequence of arguments over which is cycled, e.g.::

                facecolors=('black', 'blue')

            would create interleaving black and blue rectangles.

            Supported keywords:

            %(BrokenBarHCollection)s

        Returns
        -------
        :class:`matplotlib.collections.BrokenBarHCollection`

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        # process the unit information
        if len(xranges):
            xdata = cbook.safe_first_element(xranges)
        else:
            xdata = None
        if len(yrange):
            ydata = cbook.safe_first_element(yrange)
        else:
            ydata = None
        self._process_unit_info(xdata=xdata,
                                ydata=ydata,
                                kwargs=kwargs)
        xranges = self.convert_xunits(xranges)
        yrange = self.convert_yunits(yrange)

        col = mcoll.BrokenBarHCollection(xranges, yrange, **kwargs)
        self.add_collection(col, autolim=True)
        self.autoscale_view()

        return col

    @_preprocess_data(replace_all_args=True, label_namer=None)
    def stem(self, *args, **kwargs):
        """
        Create a stem plot.

        A stem plot plots vertical lines at each *x* location from the baseline
        to *y*, and places a marker there.

        Call signature::

          stem([x,] y, linefmt=None, markerfmt=None, basefmt=None)

        The x-positions are optional. The formats may be provided either as
        positional or as keyword-arguments.

        Parameters
        ----------
        x : array-like, optional
            The x-positions of the stems. Default: (0, 1, ..., len(y) - 1).

        y : array-like
            The y-values of the stem heads.

        linefmt : str, optional
            A string defining the properties of the vertical lines. Usually,
            this will be a color or a color and a linestyle:

            =========  =============
            Character  Line Style
            =========  =============
            ``'-'``    solid line
            ``'--'``   dashed line
            ``'-.'``   dash-dot line
            ``':'``    dotted line
            =========  =============

            Default: 'C0-', i.e. solid line with the first color of the color
            cycle.

            Note: While it is technically possible to specify valid formats
            other than color or color and linestyle (e.g. 'rx' or '-.'), this
            is beyond the intention of the method and will most likely not
            result in a reasonable reasonable plot.

        markerfmt : str, optional
            A string defining the properties of the markers at the stem heads.
            Default: 'C0o', i.e. filled circles with the first color of the
            color cycle.

        basefmt : str, optional
            A format string defining the properties of the baseline.

            Default: 'C3-' ('C2-' in classic mode).

        bottom : float, optional, default: 0
            The y-position of the baseline.

        label : str, optional, default: None
            The label to use for the stems in legends.


        Other Parameters
        ----------------
        **kwargs
            No other parameters are supported. They are currently ignored
            silently for backward compatibility. This behavior is deprecated.
            Future versions will not accept any other parameters and will
            raise a TypeError instead.


        Returns
        -------
        :class:`~matplotlib.container.StemContainer`
            The stemcontainer may be treated like a tuple
            (*markerline*, *stemlines*, *baseline*)


        Notes
        -----

        .. seealso::
            The MATLAB function
            `stem <http://www.mathworks.com/help/techdoc/ref/stem.html>`_
            which inspired this method.

        """

        # kwargs handling
        # We would like to have a signature with explicit kewords:
        # stem(*args, linefmt=None, markerfmt=None, basefmt=None,
        #      bottom=0, label=None)
        # Unfortunately,  this is not supported in Python 2.x. There, *args
        # can only exist after keyword arguments.
        linefmt = kwargs.pop('linefmt', None)
        markerfmt = kwargs.pop('markerfmt', None)
        basefmt = kwargs.pop('basefmt', None)
        bottom = kwargs.pop('bottom', None)
        if bottom is None:
            bottom = 0
        label = kwargs.pop('label', None)
        if kwargs:
            warn_deprecated(since='2.2',
                            message="stem() got an unexpected keyword "
                                    "argument '%s'. This will raise a "
                                    "TypeError in future versions." % (
                                next(k for k in kwargs), )
                            )

        remember_hold = self._hold
        if not self._hold:
            self.cla()
        self._hold = True

        # Assume there's at least one data array
        y = np.asarray(args[0])
        args = args[1:]

        # Try a second one
        try:
            second = np.asarray(args[0], dtype=float)
            x, y = y, second
            args = args[1:]
        except (IndexError, ValueError):
            # The second array doesn't make sense, or it doesn't exist
            second = np.arange(len(y))
            x = second

        # defaults for formats
        if linefmt is None:
            try:
                # fallback to positional argument
                linefmt = args[0]
            except IndexError:
                linecolor = 'C0'
                linemarker = 'None'
                linestyle = '-'
            else:
                linestyle, linemarker, linecolor = \
                    _process_plot_format(linefmt)
        else:
            linestyle, linemarker, linecolor = _process_plot_format(linefmt)

        if markerfmt is None:
            try:
                # fallback to positional argument
                markerfmt = args[1]
            except IndexError:
                markercolor = 'C0'
                markermarker = 'o'
                markerstyle = 'None'
            else:
                markerstyle, markermarker, markercolor = \
                    _process_plot_format(markerfmt)
        else:
            markerstyle, markermarker, markercolor = \
                _process_plot_format(markerfmt)

        if basefmt is None:
            try:
                # fallback to positional argument
                basefmt = args[2]
            except IndexError:
                if rcParams['_internal.classic_mode']:
                    basecolor = 'C2'
                else:
                    basecolor = 'C3'
                basemarker = 'None'
                basestyle = '-'
            else:
                basestyle, basemarker, basecolor = \
                    _process_plot_format(basefmt)
        else:
            basestyle, basemarker, basecolor = _process_plot_format(basefmt)

        markerline, = self.plot(x, y, color=markercolor, linestyle=markerstyle,
                                marker=markermarker, label="_nolegend_")

        stemlines = []
        for thisx, thisy in zip(x, y):
            l, = self.plot([thisx, thisx], [bottom, thisy],
                           color=linecolor, linestyle=linestyle,
                           marker=linemarker, label="_nolegend_")
            stemlines.append(l)

        baseline, = self.plot([np.min(x), np.max(x)], [bottom, bottom],
                              color=basecolor, linestyle=basestyle,
                              marker=basemarker, label="_nolegend_")

        self._hold = remember_hold

        stem_container = StemContainer((markerline, stemlines, baseline),
                                       label=label)
        self.add_container(stem_container)

        return stem_container

    @_preprocess_data(replace_names=["x", "explode", "labels", "colors"],
                      label_namer=None)
    def pie(self, x, explode=None, labels=None, colors=None,
            autopct=None, pctdistance=0.6, shadow=False, labeldistance=1.1,
            startangle=None, radius=None, counterclock=True,
            wedgeprops=None, textprops=None, center=(0, 0),
            frame=False, rotatelabels=False):
        """
        Plot a pie chart.

        Make a pie chart of array *x*.  The fractional area of each wedge is
        given by ``x/sum(x)``.  If ``sum(x) < 1``, then the values of *x* give
        the fractional area directly and the array will not be normalized. The
        resulting pie will have an empty wedge of size ``1 - sum(x)``.

        The wedges are plotted counterclockwise, by default starting from the
        x-axis.

        Parameters
        ----------
        x : array-like
            The wedge sizes.

        explode : array-like, optional, default: None
            If not *None*, is a ``len(x)`` array which specifies the fraction
            of the radius with which to offset each wedge.

        labels : list, optional, default: None
            A sequence of strings providing the labels for each wedge

        colors : array-like, optional, default: None
            A sequence of matplotlib color args through which the pie chart
            will cycle.  If *None*, will use the colors in the currently
            active cycle.

        autopct : None (default), string, or function, optional
            If not *None*, is a string or function used to label the wedges
            with their numeric value.  The label will be placed inside the
            wedge.  If it is a format string, the label will be ``fmt%pct``.
            If it is a function, it will be called.

        pctdistance : float, optional, default: 0.6
            The ratio between the center of each pie slice and the start of
            the text generated by *autopct*.  Ignored if *autopct* is *None*.

        shadow : bool, optional, default: False
            Draw a shadow beneath the pie.

        labeldistance : float, optional, default: 1.1
            The radial distance at which the pie labels are drawn

        startangle : float, optional, default: None
            If not *None*, rotates the start of the pie chart by *angle*
            degrees counterclockwise from the x-axis.

        radius : float, optional, default: None
            The radius of the pie, if *radius* is *None* it will be set to 1.

        counterclock : bool, optional, default: True
            Specify fractions direction, clockwise or counterclockwise.

        wedgeprops : dict, optional, default: None
            Dict of arguments passed to the wedge objects making the pie.
            For example, you can pass in ``wedgeprops = {'linewidth': 3}``
            to set the width of the wedge border lines equal to 3.
            For more details, look at the doc/arguments of the wedge object.
            By default ``clip_on=False``.

        textprops : dict, optional, default: None
            Dict of arguments to pass to the text objects.

        center :  list of float, optional, default: (0, 0)
            Center position of the chart. Takes value (0, 0) or is a sequence
            of 2 scalars.

        frame : bool, optional, default: False
            Plot axes frame with the chart if true.

        rotatelabels : bool, optional, default: False
            Rotate each label to the angle of the corresponding slice if true.

        Returns
        -------
        patches : list
            A sequence of :class:`matplotlib.patches.Wedge` instances

        texts : list
            A list of the label :class:`matplotlib.text.Text` instances.

        autotexts : list
            A list of :class:`~matplotlib.text.Text` instances for the numeric
            labels. This will only be returned if the parameter *autopct* is
            not *None*.

        Notes
        -----
        The pie chart will probably look best if the figure and axes are
        square, or the Axes aspect is equal.
        """
        x = np.array(x, np.float32)

        sx = x.sum()
        if sx > 1:
            x /= sx

        if labels is None:
            labels = [''] * len(x)
        if explode is None:
            explode = [0] * len(x)
        if len(x) != len(labels):
            raise ValueError("'label' must be of length 'x'")
        if len(x) != len(explode):
            raise ValueError("'explode' must be of length 'x'")
        if colors is None:
            get_next_color = self._get_patches_for_fill.get_next_color
        else:
            color_cycle = itertools.cycle(colors)

            def get_next_color():
                return next(color_cycle)

        if radius is None:
            radius = 1

        # Starting theta1 is the start fraction of the circle
        if startangle is None:
            theta1 = 0
        else:
            theta1 = startangle / 360.0

        # set default values in wedge_prop
        if wedgeprops is None:
            wedgeprops = {}
        wedgeprops.setdefault('clip_on', False)

        if textprops is None:
            textprops = {}
        textprops.setdefault('clip_on', False)

        texts = []
        slices = []
        autotexts = []

        i = 0
        for frac, label, expl in zip(x, labels, explode):
            x, y = center
            theta2 = (theta1 + frac) if counterclock else (theta1 - frac)
            thetam = 2 * np.pi * 0.5 * (theta1 + theta2)
            x += expl * math.cos(thetam)
            y += expl * math.sin(thetam)

            w = mpatches.Wedge((x, y), radius, 360. * min(theta1, theta2),
                               360. * max(theta1, theta2),
                               facecolor=get_next_color(),
                               **wedgeprops)
            slices.append(w)
            self.add_patch(w)
            w.set_label(label)

            if shadow:
                # make sure to add a shadow after the call to
                # add_patch so the figure and transform props will be
                # set
                shad = mpatches.Shadow(w, -0.02, -0.02)
                shad.set_zorder(0.9 * w.get_zorder())
                shad.set_label('_nolegend_')
                self.add_patch(shad)

            xt = x + labeldistance * radius * math.cos(thetam)
            yt = y + labeldistance * radius * math.sin(thetam)
            label_alignment_h = xt > 0 and 'left' or 'right'
            label_alignment_v = 'center'
            label_rotation = 'horizontal'
            if rotatelabels:
                label_alignment_v = yt > 0 and 'bottom' or 'top'
                label_rotation = np.rad2deg(thetam) + (0 if xt > 0 else 180)

            t = self.text(xt, yt, label,
                          size=rcParams['xtick.labelsize'],
                          horizontalalignment=label_alignment_h,
                          verticalalignment=label_alignment_v,
                          rotation=label_rotation,
                          **textprops)

            texts.append(t)

            if autopct is not None:
                xt = x + pctdistance * radius * math.cos(thetam)
                yt = y + pctdistance * radius * math.sin(thetam)
                if isinstance(autopct, six.string_types):
                    s = autopct % (100. * frac)
                elif callable(autopct):
                    s = autopct(100. * frac)
                else:
                    raise TypeError(
                        'autopct must be callable or a format string')

                t = self.text(xt, yt, s,
                              horizontalalignment='center',
                              verticalalignment='center',
                              **textprops)

                autotexts.append(t)

            theta1 = theta2
            i += 1

        if not frame:
            self.set_frame_on(False)

            self.set_xlim((-1.25 + center[0],
                           1.25 + center[0]))
            self.set_ylim((-1.25 + center[1],
                           1.25 + center[1]))
            self.set_xticks([])
            self.set_yticks([])

        if autopct is None:
            return slices, texts
        else:
            return slices, texts, autotexts

    @_preprocess_data(replace_names=["x", "y", "xerr", "yerr"],
                      label_namer="y")
    @docstring.dedent_interpd
    def errorbar(self, x, y, yerr=None, xerr=None,
                 fmt='', ecolor=None, elinewidth=None, capsize=None,
                 barsabove=False, lolims=False, uplims=False,
                 xlolims=False, xuplims=False, errorevery=1, capthick=None,
                 **kwargs):
        """
        Plot y versus x as lines and/or markers with attached errorbars.

        *x*, *y* define the data locations, *xerr*, *yerr* define the errorbar
        sizes. By default, this draws the data markers/lines as well the
        errorbars. Use fmt='none' to draw errorbars without any data markers.

        Parameters
        ----------
        x, y : scalar or array-like
            The data positions.

        xerr, yerr : scalar or array-like, shape(N,) or shape(2,N), optional
            The errorbar sizes:

            - scalar: Symmetric +/- values for all data points.
            - shape(N,): Symmetric +/-values for each data point.
            - shape(2,N): Separate + and - values for each data point.
            - *None*: No errorbar.

        fmt : plot format string, optional, default: ''
            The format for the data points / data lines. See `.plot` for
            details.

            Use 'none' (case insensitive) to plot errorbars without any data
            markers.

        ecolor : mpl color, optional, default: None
            A matplotlib color arg which gives the color the errorbar lines.
            If None, use the color of the line connecting the markers.

        elinewidth : scalar, optional, default: None
            The linewidth of the errorbar lines. If None, the linewidth of
            the current style is used.

        capsize : scalar, optional, default: None
            The length of the error bar caps in points. If None, it will take
            the value from :rc:`errorbar.capsize`.

        capthick : scalar, optional, default: None
            An alias to the keyword argument *markeredgewidth* (a.k.a. *mew*).
            This setting is a more sensible name for the property that
            controls the thickness of the error bar cap in points. For
            backwards compatibility, if *mew* or *markeredgewidth* are given,
            then they will over-ride *capthick*. This may change in future
            releases.

        barsabove : bool, optional, default: False
            If True, will plot the errorbars above the plot
            symbols. Default is below.

        lolims, uplims, xlolims, xuplims : bool, optional, default: None
            These arguments can be used to indicate that a value gives only
            upper/lower limits. In that case a caret symbol is used to
            indicate this. *lims*-arguments may be of the same type as *xerr*
            and *yerr*.  To use limits with inverted axes, :meth:`set_xlim`
            or :meth:`set_ylim` must be called before :meth:`errorbar`.

        errorevery : positive integer, optional, default: 1
            Subsamples the errorbars. e.g., if errorevery=5, errorbars for
            every 5-th datapoint will be plotted. The data plot itself still
            shows all data points.

        Returns
        -------
        :class:`~.container.ErrorbarContainer`
            The container contains:

            - plotline: :class:`~matplotlib.lines.Line2D` instance of
              x, y plot markers and/or line.
            - caplines: A tuple of :class:`~matplotlib.lines.Line2D` instances
              of the error bar caps.
            - barlinecols: A tuple of
              :class:`~matplotlib.collections.LineCollection` with the
              horizontal and vertical error ranges.

        Other Parameters
        ----------------
        **kwargs :
            All other keyword arguments are passed on to the plot
            command for the markers. For example, this code makes big red
            squares with thick green edges::

                x,y,yerr = rand(3,10)
                errorbar(x, y, yerr, marker='s', mfc='red',
                         mec='green', ms=20, mew=4)

            where *mfc*, *mec*, *ms* and *mew* are aliases for the longer
            property names, *markerfacecolor*, *markeredgecolor*, *markersize*
            and *markeredgewidth*.

            Valid kwargs for the marker properties are `.Lines2D` properties:

            %(Line2D)s

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        kwargs = cbook.normalize_kwargs(kwargs, _alias_map)
        # anything that comes in as 'None', drop so the default thing
        # happens down stream
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        kwargs.setdefault('zorder', 2)

        if errorevery < 1:
            raise ValueError(
                'errorevery has to be a strictly positive integer')

        self._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)
        if not self._hold:
            self.cla()
        holdstate = self._hold
        self._hold = True

        plot_line = (fmt.lower() != 'none')
        label = kwargs.pop("label", None)

        if fmt == '':
            fmt_style_kwargs = {}
        else:
            fmt_style_kwargs = {k: v for k, v in
                            zip(('linestyle', 'marker', 'color'),
                                _process_plot_format(fmt)) if v is not None}
        if fmt == 'none':
            # Remove alpha=0 color that _process_plot_format returns
            fmt_style_kwargs.pop('color')

        if ('color' in kwargs or 'color' in fmt_style_kwargs or
                ecolor is not None):
            base_style = {}
            if 'color' in kwargs:
                base_style['color'] = kwargs.pop('color')
        else:
            base_style = next(self._get_lines.prop_cycler)

        base_style['label'] = '_nolegend_'
        base_style.update(fmt_style_kwargs)
        if 'color' not in base_style:
            base_style['color'] = 'C0'
        if ecolor is None:
            ecolor = base_style['color']
        # make sure all the args are iterable; use lists not arrays to
        # preserve units
        if not iterable(x):
            x = [x]

        if not iterable(y):
            y = [y]

        if xerr is not None:
            if not iterable(xerr):
                xerr = [xerr] * len(x)

        if yerr is not None:
            if not iterable(yerr):
                yerr = [yerr] * len(y)

        # make the style dict for the 'normal' plot line
        plot_line_style = dict(base_style)
        plot_line_style.update(**kwargs)
        if barsabove:
            plot_line_style['zorder'] = kwargs['zorder'] - .1
        else:
            plot_line_style['zorder'] = kwargs['zorder'] + .1

        # make the style dict for the line collections (the bars)
        eb_lines_style = dict(base_style)
        eb_lines_style.pop('marker', None)
        eb_lines_style.pop('linestyle', None)
        eb_lines_style['color'] = ecolor

        if elinewidth:
            eb_lines_style['linewidth'] = elinewidth
        elif 'linewidth' in kwargs:
            eb_lines_style['linewidth'] = kwargs['linewidth']

        for key in ('transform', 'alpha', 'zorder', 'rasterized'):
            if key in kwargs:
                eb_lines_style[key] = kwargs[key]

        # set up cap style dictionary
        eb_cap_style = dict(base_style)
        # eject any marker information from format string
        eb_cap_style.pop('marker', None)
        eb_lines_style.pop('markerfacecolor', None)
        eb_lines_style.pop('markeredgewidth', None)
        eb_lines_style.pop('markeredgecolor', None)
        eb_cap_style.pop('ls', None)
        eb_cap_style['linestyle'] = 'none'
        if capsize is None:
            capsize = rcParams["errorbar.capsize"]
        if capsize > 0:
            eb_cap_style['markersize'] = 2. * capsize
        if capthick is not None:
            eb_cap_style['markeredgewidth'] = capthick

        # For backwards-compat, allow explicit setting of
        # 'markeredgewidth' to over-ride capthick.
        for key in ('markeredgewidth', 'transform', 'alpha',
                    'zorder', 'rasterized'):
            if key in kwargs:
                eb_cap_style[key] = kwargs[key]
        eb_cap_style['color'] = ecolor

        data_line = None
        if plot_line:
            data_line = mlines.Line2D(x, y, **plot_line_style)
            self.add_line(data_line)

        barcols = []
        caplines = []

        # arrays fine here, they are booleans and hence not units
        def _bool_asarray_helper(d, expected):
            if not iterable(d):
                return np.asarray([d] * expected, bool)
            else:
                return np.asarray(d, bool)

        lolims = _bool_asarray_helper(lolims, len(x))
        uplims = _bool_asarray_helper(uplims, len(x))
        xlolims = _bool_asarray_helper(xlolims, len(x))
        xuplims = _bool_asarray_helper(xuplims, len(x))

        everymask = np.arange(len(x)) % errorevery == 0

        def xywhere(xs, ys, mask):
            """
            return xs[mask], ys[mask] where mask is True but xs and
            ys are not arrays
            """
            assert len(xs) == len(ys)
            assert len(xs) == len(mask)
            xs = [thisx for thisx, b in zip(xs, mask) if b]
            ys = [thisy for thisy, b in zip(ys, mask) if b]
            return xs, ys

        def extract_err(err, data):
            '''private function to compute error bars

            Parameters
            ----------
            err : iterable
                xerr or yerr from errorbar
            data : iterable
                x or y from errorbar
            '''
            try:
                a, b = err
            except (TypeError, ValueError):
                pass
            else:
                if iterable(a) and iterable(b):
                    # using list comps rather than arrays to preserve units
                    low = [thisx - thiserr for (thisx, thiserr)
                           in cbook.safezip(data, a)]
                    high = [thisx + thiserr for (thisx, thiserr)
                            in cbook.safezip(data, b)]
                    return low, high
            # Check if xerr is scalar or symmetric. Asymmetric is handled
            # above. This prevents Nx2 arrays from accidentally
            # being accepted, when the user meant the 2xN transpose.
            # special case for empty lists
            if len(err) > 1:
                fe = safe_first_element(err)
                if (len(err) != len(data) or np.size(fe) > 1):
                    raise ValueError("err must be [ scalar | N, Nx1 "
                                     "or 2xN array-like ]")
            # using list comps rather than arrays to preserve units
            low = [thisx - thiserr for (thisx, thiserr)
                   in cbook.safezip(data, err)]
            high = [thisx + thiserr for (thisx, thiserr)
                    in cbook.safezip(data, err)]
            return low, high

        if xerr is not None:
            left, right = extract_err(xerr, x)
            # select points without upper/lower limits in x and
            # draw normal errorbars for these points
            noxlims = ~(xlolims | xuplims)
            if noxlims.any() or len(noxlims) == 0:
                yo, _ = xywhere(y, right, noxlims & everymask)
                lo, ro = xywhere(left, right, noxlims & everymask)
                barcols.append(self.hlines(yo, lo, ro, **eb_lines_style))
                if capsize > 0:
                    caplines.append(mlines.Line2D(lo, yo, marker='|',
                                                  **eb_cap_style))
                    caplines.append(mlines.Line2D(ro, yo, marker='|',
                                                  **eb_cap_style))

            if xlolims.any():
                yo, _ = xywhere(y, right, xlolims & everymask)
                lo, ro = xywhere(x, right, xlolims & everymask)
                barcols.append(self.hlines(yo, lo, ro, **eb_lines_style))
                rightup, yup = xywhere(right, y, xlolims & everymask)
                if self.xaxis_inverted():
                    marker = mlines.CARETLEFTBASE
                else:
                    marker = mlines.CARETRIGHTBASE
                caplines.append(
                    mlines.Line2D(rightup, yup, ls='None', marker=marker,
                                  **eb_cap_style))
                if capsize > 0:
                    xlo, ylo = xywhere(x, y, xlolims & everymask)
                    caplines.append(mlines.Line2D(xlo, ylo, marker='|',
                                                  **eb_cap_style))

            if xuplims.any():
                yo, _ = xywhere(y, right, xuplims & everymask)
                lo, ro = xywhere(left, x, xuplims & everymask)
                barcols.append(self.hlines(yo, lo, ro, **eb_lines_style))
                leftlo, ylo = xywhere(left, y, xuplims & everymask)
                if self.xaxis_inverted():
                    marker = mlines.CARETRIGHTBASE
                else:
                    marker = mlines.CARETLEFTBASE
                caplines.append(
                    mlines.Line2D(leftlo, ylo, ls='None', marker=marker,
                                  **eb_cap_style))
                if capsize > 0:
                    xup, yup = xywhere(x, y, xuplims & everymask)
                    caplines.append(mlines.Line2D(xup, yup, marker='|',
                                                  **eb_cap_style))

        if yerr is not None:
            lower, upper = extract_err(yerr, y)
            # select points without upper/lower limits in y and
            # draw normal errorbars for these points
            noylims = ~(lolims | uplims)
            if noylims.any() or len(noylims) == 0:
                xo, _ = xywhere(x, lower, noylims & everymask)
                lo, uo = xywhere(lower, upper, noylims & everymask)
                barcols.append(self.vlines(xo, lo, uo, **eb_lines_style))
                if capsize > 0:
                    caplines.append(mlines.Line2D(xo, lo, marker='_',
                                                  **eb_cap_style))
                    caplines.append(mlines.Line2D(xo, uo, marker='_',
                                                  **eb_cap_style))

            if lolims.any():
                xo, _ = xywhere(x, lower, lolims & everymask)
                lo, uo = xywhere(y, upper, lolims & everymask)
                barcols.append(self.vlines(xo, lo, uo, **eb_lines_style))
                xup, upperup = xywhere(x, upper, lolims & everymask)
                if self.yaxis_inverted():
                    marker = mlines.CARETDOWNBASE
                else:
                    marker = mlines.CARETUPBASE
                caplines.append(
                    mlines.Line2D(xup, upperup, ls='None', marker=marker,
                                  **eb_cap_style))
                if capsize > 0:
                    xlo, ylo = xywhere(x, y, lolims & everymask)
                    caplines.append(mlines.Line2D(xlo, ylo, marker='_',
                                                  **eb_cap_style))

            if uplims.any():
                xo, _ = xywhere(x, lower, uplims & everymask)
                lo, uo = xywhere(lower, y, uplims & everymask)
                barcols.append(self.vlines(xo, lo, uo, **eb_lines_style))
                xlo, lowerlo = xywhere(x, lower, uplims & everymask)
                if self.yaxis_inverted():
                    marker = mlines.CARETUPBASE
                else:
                    marker = mlines.CARETDOWNBASE
                caplines.append(
                    mlines.Line2D(xlo, lowerlo, ls='None', marker=marker,
                                  **eb_cap_style))
                if capsize > 0:
                    xup, yup = xywhere(x, y, uplims & everymask)
                    caplines.append(mlines.Line2D(xup, yup, marker='_',
                                                  **eb_cap_style))
        for l in caplines:
            self.add_line(l)

        self.autoscale_view()
        self._hold = holdstate

        errorbar_container = ErrorbarContainer((data_line, tuple(caplines),
                                                tuple(barcols)),
                                               has_xerr=(xerr is not None),
                                               has_yerr=(yerr is not None),
                                               label=label)
        self.containers.append(errorbar_container)

        return errorbar_container  # (l0, caplines, barcols)

    @_preprocess_data(label_namer=None)
    def boxplot(self, x, notch=None, sym=None, vert=None, whis=None,
                positions=None, widths=None, patch_artist=None,
                bootstrap=None, usermedians=None, conf_intervals=None,
                meanline=None, showmeans=None, showcaps=None,
                showbox=None, showfliers=None, boxprops=None,
                labels=None, flierprops=None, medianprops=None,
                meanprops=None, capprops=None, whiskerprops=None,
                manage_xticks=True, autorange=False, zorder=None):
        """
        Make a box and whisker plot.

        Make a box and whisker plot for each column of ``x`` or each
        vector in sequence ``x``.  The box extends from the lower to
        upper quartile values of the data, with a line at the median.
        The whiskers extend from the box to show the range of the
        data.  Flier points are those past the end of the whiskers.

        Parameters
        ----------
        x : Array or a sequence of vectors.
            The input data.

        notch : bool, optional (False)
            If `True`, will produce a notched box plot. Otherwise, a
            rectangular boxplot is produced. The notches represent the
            confidence interval (CI) around the median. See the entry
            for the ``bootstrap`` parameter for information regarding
            how the locations of the notches are computed.

            .. note::

                In cases where the values of the CI are less than the
                lower quartile or greater than the upper quartile, the
                notches will extend beyond the box, giving it a
                distinctive "flipped" appearance. This is expected
                behavior and consistent with other statistical
                visualization packages.

        sym : str, optional
            The default symbol for flier points. Enter an empty string
            ('') if you don't want to show fliers. If `None`, then the
            fliers default to 'b+'  If you want more control use the
            flierprops kwarg.

        vert : bool, optional (True)
            If `True` (default), makes the boxes vertical. If `False`,
            everything is drawn horizontally.

        whis : float, sequence, or string (default = 1.5)
            As a float, determines the reach of the whiskers to the beyond the
            first and third quartiles. In other words, where IQR is the
            interquartile range (`Q3-Q1`), the upper whisker will extend to
            last datum less than `Q3 + whis*IQR`). Similarly, the lower whisker
            will extend to the first datum greater than `Q1 - whis*IQR`.
            Beyond the whiskers, data
            are considered outliers and are plotted as individual
            points. Set this to an unreasonably high value to force the
            whiskers to show the min and max values. Alternatively, set
            this to an ascending sequence of percentile (e.g., [5, 95])
            to set the whiskers at specific percentiles of the data.
            Finally, ``whis`` can be the string ``'range'`` to force the
            whiskers to the min and max of the data.

        bootstrap : int, optional
            Specifies whether to bootstrap the confidence intervals
            around the median for notched boxplots. If ``bootstrap`` is
            None, no bootstrapping is performed, and notches are
            calculated using a Gaussian-based asymptotic approximation
            (see McGill, R., Tukey, J.W., and Larsen, W.A., 1978, and
            Kendall and Stuart, 1967). Otherwise, bootstrap specifies
            the number of times to bootstrap the median to determine its
            95% confidence intervals. Values between 1000 and 10000 are
            recommended.

        usermedians : array-like, optional
            An array or sequence whose first dimension (or length) is
            compatible with ``x``. This overrides the medians computed
            by matplotlib for each element of ``usermedians`` that is not
            `None`. When an element of ``usermedians`` is None, the median
            will be computed by matplotlib as normal.

        conf_intervals : array-like, optional
            Array or sequence whose first dimension (or length) is
            compatible with ``x`` and whose second dimension is 2. When
            the an element of ``conf_intervals`` is not None, the
            notch locations computed by matplotlib are overridden
            (provided ``notch`` is `True`). When an element of
            ``conf_intervals`` is `None`, the notches are computed by the
            method specified by the other kwargs (e.g., ``bootstrap``).

        positions : array-like, optional
            Sets the positions of the boxes. The ticks and limits are
            automatically set to match the positions. Defaults to
            `range(1, N+1)` where N is the number of boxes to be drawn.

        widths : scalar or array-like
            Sets the width of each box either with a scalar or a
            sequence. The default is 0.5, or ``0.15*(distance between
            extreme positions)``, if that is smaller.

        patch_artist : bool, optional (False)
            If `False` produces boxes with the Line2D artist. Otherwise,
            boxes and drawn with Patch artists.

        labels : sequence, optional
            Labels for each dataset. Length must be compatible with
            dimensions  of ``x``.

        manage_xticks : bool, optional (True)
            If the function should adjust the xlim and xtick locations.

        autorange : bool, optional (False)
            When `True` and the data are distributed such that the  25th and
            75th percentiles are equal, ``whis`` is set to ``'range'`` such
            that the whisker ends are at the minimum and maximum of the
            data.

        meanline : bool, optional (False)
            If `True` (and ``showmeans`` is `True`), will try to render
            the mean as a line spanning the full width of the box
            according to ``meanprops`` (see below). Not recommended if
            ``shownotches`` is also True. Otherwise, means will be shown
            as points.

        zorder : scalar, optional (None)
            Sets the zorder of the boxplot.

        Other Parameters
        ----------------
        showcaps : bool, optional (True)
            Show the caps on the ends of whiskers.
        showbox : bool, optional (True)
            Show the central box.
        showfliers : bool, optional (True)
            Show the outliers beyond the caps.
        showmeans : bool, optional (False)
            Show the arithmetic means.
        capprops : dict, optional (None)
            Specifies the style of the caps.
        boxprops : dict, optional (None)
            Specifies the style of the box.
        whiskerprops : dict, optional (None)
            Specifies the style of the whiskers.
        flierprops : dict, optional (None)
            Specifies the style of the fliers.
        medianprops : dict, optional (None)
            Specifies the style of the median.
        meanprops : dict, optional (None)
            Specifies the style of the mean.

        Returns
        -------
        result : dict
          A dictionary mapping each component of the boxplot to a list
          of the :class:`matplotlib.lines.Line2D` instances
          created. That dictionary has the following keys (assuming
          vertical boxplots):

          - ``boxes``: the main body of the boxplot showing the
            quartiles and the median's confidence intervals if
            enabled.

          - ``medians``: horizontal lines at the median of each box.

          - ``whiskers``: the vertical lines extending to the most
            extreme, non-outlier data points.

          - ``caps``: the horizontal lines at the ends of the
            whiskers.

          - ``fliers``: points representing data that extend beyond
            the whiskers (fliers).

          - ``means``: points or lines representing the means.

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """

        # If defined in matplotlibrc, apply the value from rc file
        # Overridden if argument is passed
        if whis is None:
            whis = rcParams['boxplot.whiskers']
        if bootstrap is None:
            bootstrap = rcParams['boxplot.bootstrap']

        bxpstats = cbook.boxplot_stats(x, whis=whis, bootstrap=bootstrap,
                                       labels=labels, autorange=autorange)
        if notch is None:
            notch = rcParams['boxplot.notch']
        if vert is None:
            vert = rcParams['boxplot.vertical']
        if patch_artist is None:
            patch_artist = rcParams['boxplot.patchartist']
        if meanline is None:
            meanline = rcParams['boxplot.meanline']
        if showmeans is None:
            showmeans = rcParams['boxplot.showmeans']
        if showcaps is None:
            showcaps = rcParams['boxplot.showcaps']
        if showbox is None:
            showbox = rcParams['boxplot.showbox']
        if showfliers is None:
            showfliers = rcParams['boxplot.showfliers']

        def _update_dict(dictionary, rc_name, properties):
            """ Loads properties in the dictionary from rc file if not already
            in the dictionary"""
            rc_str = 'boxplot.{0}.{1}'
            if dictionary is None:
                dictionary = dict()
            for prop_dict in properties:
                dictionary.setdefault(prop_dict,
                                rcParams[rc_str.format(rc_name, prop_dict)])
            return dictionary

        # Common property dictionnaries loading from rc
        flier_props = ['color', 'marker', 'markerfacecolor', 'markeredgecolor',
                       'markersize', 'linestyle', 'linewidth']
        default_props = ['color', 'linewidth', 'linestyle']

        boxprops = _update_dict(boxprops, 'boxprops', default_props)
        whiskerprops = _update_dict(whiskerprops, 'whiskerprops',
                                                            default_props)
        capprops = _update_dict(capprops, 'capprops', default_props)
        medianprops = _update_dict(medianprops, 'medianprops', default_props)
        meanprops = _update_dict(meanprops, 'meanprops', default_props)
        flierprops = _update_dict(flierprops, 'flierprops', flier_props)

        if patch_artist:
            boxprops['linestyle'] = 'solid'
            boxprops['edgecolor'] = boxprops.pop('color')

        # if non-default sym value, put it into the flier dictionary
        # the logic for providing the default symbol ('b+') now lives
        # in bxp in the initial value of final_flierprops
        # handle all of the `sym` related logic here so we only have to pass
        # on the flierprops dict.
        if sym is not None:
            # no-flier case, which should really be done with
            # 'showfliers=False' but none-the-less deal with it to keep back
            # compatibility
            if sym == '':
                # blow away existing dict and make one for invisible markers
                flierprops = dict(linestyle='none', marker='', color='none')
                # turn the fliers off just to be safe
                showfliers = False
            # now process the symbol string
            else:
                # process the symbol string
                # discarded linestyle
                _, marker, color = _process_plot_format(sym)
                # if we have a marker, use it
                if marker is not None:
                    flierprops['marker'] = marker
                # if we have a color, use it
                if color is not None:
                    # assume that if color is passed in the user want
                    # filled symbol, if the users want more control use
                    # flierprops
                    flierprops['color'] = color
                    flierprops['markerfacecolor'] = color
                    flierprops['markeredgecolor'] = color

        # replace medians if necessary:
        if usermedians is not None:
            if (len(np.ravel(usermedians)) != len(bxpstats) or
                    np.shape(usermedians)[0] != len(bxpstats)):
                raise ValueError('usermedians length not compatible with x')
            else:
                # reassign medians as necessary
                for stats, med in zip(bxpstats, usermedians):
                    if med is not None:
                        stats['med'] = med

        if conf_intervals is not None:
            if np.shape(conf_intervals)[0] != len(bxpstats):
                err_mess = 'conf_intervals length not compatible with x'
                raise ValueError(err_mess)
            else:
                for stats, ci in zip(bxpstats, conf_intervals):
                    if ci is not None:
                        if len(ci) != 2:
                            raise ValueError('each confidence interval must '
                                             'have two values')
                        else:
                            if ci[0] is not None:
                                stats['cilo'] = ci[0]
                            if ci[1] is not None:
                                stats['cihi'] = ci[1]

        artists = self.bxp(bxpstats, positions=positions, widths=widths,
                           vert=vert, patch_artist=patch_artist,
                           shownotches=notch, showmeans=showmeans,
                           showcaps=showcaps, showbox=showbox,
                           boxprops=boxprops, flierprops=flierprops,
                           medianprops=medianprops, meanprops=meanprops,
                           meanline=meanline, showfliers=showfliers,
                           capprops=capprops, whiskerprops=whiskerprops,
                           manage_xticks=manage_xticks, zorder=zorder)
        return artists

    def bxp(self, bxpstats, positions=None, widths=None, vert=True,
            patch_artist=False, shownotches=False, showmeans=False,
            showcaps=True, showbox=True, showfliers=True,
            boxprops=None, whiskerprops=None, flierprops=None,
            medianprops=None, capprops=None, meanprops=None,
            meanline=False, manage_xticks=True, zorder=None):
        """
        Drawing function for box and whisker plots.

        Make a box and whisker plot for each column of *x* or each
        vector in sequence *x*.  The box extends from the lower to
        upper quartile values of the data, with a line at the median.
        The whiskers extend from the box to show the range of the
        data.  Flier points are those past the end of the whiskers.

        Parameters
        ----------

        bxpstats : list of dicts
          A list of dictionaries containing stats for each boxplot.
          Required keys are:

          - ``med``: The median (scalar float).

          - ``q1``: The first quartile (25th percentile) (scalar
            float).

          - ``q3``: The third quartile (75th percentile) (scalar
            float).

          - ``whislo``: Lower bound of the lower whisker (scalar
            float).

          - ``whishi``: Upper bound of the upper whisker (scalar
            float).

          Optional keys are:

          - ``mean``: The mean (scalar float). Needed if
            ``showmeans=True``.

          - ``fliers``: Data beyond the whiskers (sequence of floats).
            Needed if ``showfliers=True``.

          - ``cilo`` & ``cihi``: Lower and upper confidence intervals
            about the median. Needed if ``shownotches=True``.

          - ``label``: Name of the dataset (string). If available,
            this will be used a tick label for the boxplot

        positions : array-like, default = [1, 2, ..., n]
          Sets the positions of the boxes. The ticks and limits
          are automatically set to match the positions.

        widths : array-like, default = None
          Either a scalar or a vector and sets the width of each
          box. The default is ``0.15*(distance between extreme
          positions)``, clipped to no less than 0.15 and no more than
          0.5.

        vert : bool, default = False
          If `True` (default), makes the boxes vertical.  If `False`,
          makes horizontal boxes.

        patch_artist : bool, default = False
          If `False` produces boxes with the
          `~matplotlib.lines.Line2D` artist.  If `True` produces boxes
          with the `~matplotlib.patches.Patch` artist.

        shownotches : bool, default = False
          If `False` (default), produces a rectangular box plot.
          If `True`, will produce a notched box plot

        showmeans : bool, default = False
          If `True`, will toggle on the rendering of the means

        showcaps  : bool, default = True
          If `True`, will toggle on the rendering of the caps

        showbox  : bool, default = True
          If `True`, will toggle on the rendering of the box

        showfliers : bool, default = True
          If `True`, will toggle on the rendering of the fliers

        boxprops : dict or None (default)
          If provided, will set the plotting style of the boxes

        whiskerprops : dict or None (default)
          If provided, will set the plotting style of the whiskers

        capprops : dict or None (default)
          If provided, will set the plotting style of the caps

        flierprops : dict or None (default)
          If provided will set the plotting style of the fliers

        medianprops : dict or None (default)
          If provided, will set the plotting style of the medians

        meanprops : dict or None (default)
          If provided, will set the plotting style of the means

        meanline : bool, default = False
          If `True` (and *showmeans* is `True`), will try to render the mean
          as a line spanning the full width of the box according to
          *meanprops*. Not recommended if *shownotches* is also True.
          Otherwise, means will be shown as points.

        manage_xticks : bool, default = True
          If the function should adjust the xlim and xtick locations.

        zorder : scalar,  default = None
          The zorder of the resulting boxplot

        Returns
        -------
        result : dict
          A dictionary mapping each component of the boxplot to a list
          of the :class:`matplotlib.lines.Line2D` instances
          created. That dictionary has the following keys (assuming
          vertical boxplots):

          - ``boxes``: the main body of the boxplot showing the
            quartiles and the median's confidence intervals if
            enabled.

          - ``medians``: horizontal lines at the median of each box.

          - ``whiskers``: the vertical lines extending to the most
            extreme, non-outlier data points.

          - ``caps``: the horizontal lines at the ends of the
            whiskers.

          - ``fliers``: points representing data that extend beyond
            the whiskers (fliers).

          - ``means``: points or lines representing the means.

        Examples
        --------

        .. plot:: gallery/statistics/bxp.py

        """
        # lists of artists to be output
        whiskers = []
        caps = []
        boxes = []
        medians = []
        means = []
        fliers = []

        # empty list of xticklabels
        datalabels = []

        # Use default zorder if none specified
        if zorder is None:
            zorder = mlines.Line2D.zorder

        zdelta = 0.1
        # box properties
        if patch_artist:
            final_boxprops = dict(
                linestyle=rcParams['boxplot.boxprops.linestyle'],
                edgecolor=rcParams['boxplot.boxprops.color'],
                facecolor=rcParams['patch.facecolor'],
                linewidth=rcParams['boxplot.boxprops.linewidth']
            )
            if rcParams['_internal.classic_mode']:
                final_boxprops['facecolor'] = 'white'
        else:
            final_boxprops = dict(
                linestyle=rcParams['boxplot.boxprops.linestyle'],
                color=rcParams['boxplot.boxprops.color'],
            )

        final_boxprops['zorder'] = zorder
        if boxprops is not None:
            final_boxprops.update(boxprops)

        # other (cap, whisker) properties
        final_whiskerprops = dict(
            linestyle=rcParams['boxplot.whiskerprops.linestyle'],
            linewidth=rcParams['boxplot.whiskerprops.linewidth'],
            color=rcParams['boxplot.whiskerprops.color'],
        )

        final_capprops = dict(
            linestyle=rcParams['boxplot.capprops.linestyle'],
            linewidth=rcParams['boxplot.capprops.linewidth'],
            color=rcParams['boxplot.capprops.color'],
        )

        final_capprops['zorder'] = zorder
        if capprops is not None:
            final_capprops.update(capprops)

        final_whiskerprops['zorder'] = zorder
        if whiskerprops is not None:
            final_whiskerprops.update(whiskerprops)

        # set up the default flier properties
        final_flierprops = dict(
            linestyle=rcParams['boxplot.flierprops.linestyle'],
            linewidth=rcParams['boxplot.flierprops.linewidth'],
            color=rcParams['boxplot.flierprops.color'],
            marker=rcParams['boxplot.flierprops.marker'],
            markerfacecolor=rcParams['boxplot.flierprops.markerfacecolor'],
            markeredgecolor=rcParams['boxplot.flierprops.markeredgecolor'],
            markersize=rcParams['boxplot.flierprops.markersize'],
        )

        final_flierprops['zorder'] = zorder
        # flier (outlier) properties
        if flierprops is not None:
            final_flierprops.update(flierprops)

        # median line properties
        final_medianprops = dict(
            linestyle=rcParams['boxplot.medianprops.linestyle'],
            linewidth=rcParams['boxplot.medianprops.linewidth'],
            color=rcParams['boxplot.medianprops.color'],
        )
        final_medianprops['zorder'] = zorder + zdelta
        if medianprops is not None:
            final_medianprops.update(medianprops)

        # mean (line or point) properties
        if meanline:
            final_meanprops = dict(
                linestyle=rcParams['boxplot.meanprops.linestyle'],
                linewidth=rcParams['boxplot.meanprops.linewidth'],
                color=rcParams['boxplot.meanprops.color'],
            )
        else:
            final_meanprops = dict(
                linestyle='',
                marker=rcParams['boxplot.meanprops.marker'],
                markerfacecolor=rcParams['boxplot.meanprops.markerfacecolor'],
                markeredgecolor=rcParams['boxplot.meanprops.markeredgecolor'],
                markersize=rcParams['boxplot.meanprops.markersize'],
            )
        final_meanprops['zorder'] = zorder + zdelta
        if meanprops is not None:
            final_meanprops.update(meanprops)

        def to_vc(xs, ys):
            # convert arguments to verts and codes, append (0, 0) (ignored).
            verts = np.append(np.column_stack([xs, ys]), [(0, 0)], 0)
            codes = ([mpath.Path.MOVETO]
                     + [mpath.Path.LINETO] * (len(verts) - 2)
                     + [mpath.Path.CLOSEPOLY])
            return verts, codes

        def patch_list(xs, ys, **kwargs):
            verts, codes = to_vc(xs, ys)
            path = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(path, **kwargs)
            self.add_artist(patch)
            return [patch]

        # vertical or horizontal plot?
        if vert:
            def doplot(*args, **kwargs):
                return self.plot(*args, **kwargs)

            def dopatch(xs, ys, **kwargs):
                return patch_list(xs, ys, **kwargs)

        else:
            def doplot(*args, **kwargs):
                shuffled = []
                for i in xrange(0, len(args), 2):
                    shuffled.extend([args[i + 1], args[i]])
                return self.plot(*shuffled, **kwargs)

            def dopatch(xs, ys, **kwargs):
                xs, ys = ys, xs  # flip X, Y
                return patch_list(xs, ys, **kwargs)

        # input validation
        N = len(bxpstats)
        datashape_message = ("List of boxplot statistics and `{0}` "
                             "values must have same the length")
        # check position
        if positions is None:
            positions = list(xrange(1, N + 1))
        elif len(positions) != N:
            raise ValueError(datashape_message.format("positions"))

        # width
        if widths is None:
            widths = [np.clip(0.15 * np.ptp(positions), 0.15, 0.5)] * N
        elif np.isscalar(widths):
            widths = [widths] * N
        elif len(widths) != N:
            raise ValueError(datashape_message.format("widths"))

        # check and save the `hold` state of the current axes
        if not self._hold:
            self.cla()
        holdStatus = self._hold
        for pos, width, stats in zip(positions, widths, bxpstats):
            # try to find a new label
            datalabels.append(stats.get('label', pos))

            # whisker coords
            whisker_x = np.ones(2) * pos
            whiskerlo_y = np.array([stats['q1'], stats['whislo']])
            whiskerhi_y = np.array([stats['q3'], stats['whishi']])

            # cap coords
            cap_left = pos - width * 0.25
            cap_right = pos + width * 0.25
            cap_x = np.array([cap_left, cap_right])
            cap_lo = np.ones(2) * stats['whislo']
            cap_hi = np.ones(2) * stats['whishi']

            # box and median coords
            box_left = pos - width * 0.5
            box_right = pos + width * 0.5
            med_y = [stats['med'], stats['med']]

            # notched boxes
            if shownotches:
                box_x = [box_left, box_right, box_right, cap_right, box_right,
                         box_right, box_left, box_left, cap_left, box_left,
                         box_left]
                box_y = [stats['q1'], stats['q1'], stats['cilo'],
                         stats['med'], stats['cihi'], stats['q3'],
                         stats['q3'], stats['cihi'], stats['med'],
                         stats['cilo'], stats['q1']]
                med_x = cap_x

            # plain boxes
            else:
                box_x = [box_left, box_right, box_right, box_left, box_left]
                box_y = [stats['q1'], stats['q1'], stats['q3'], stats['q3'],
                         stats['q1']]
                med_x = [box_left, box_right]

            # maybe draw the box:
            if showbox:
                if patch_artist:
                    boxes.extend(dopatch(box_x, box_y, **final_boxprops))
                else:
                    boxes.extend(doplot(box_x, box_y, **final_boxprops))

            # draw the whiskers
            whiskers.extend(doplot(
                whisker_x, whiskerlo_y, **final_whiskerprops
            ))
            whiskers.extend(doplot(
                whisker_x, whiskerhi_y, **final_whiskerprops
            ))

            # maybe draw the caps:
            if showcaps:
                caps.extend(doplot(cap_x, cap_lo, **final_capprops))
                caps.extend(doplot(cap_x, cap_hi, **final_capprops))

            # draw the medians
            medians.extend(doplot(med_x, med_y, **final_medianprops))

            # maybe draw the means
            if showmeans:
                if meanline:
                    means.extend(doplot(
                        [box_left, box_right], [stats['mean'], stats['mean']],
                        **final_meanprops
                    ))
                else:
                    means.extend(doplot(
                        [pos], [stats['mean']], **final_meanprops
                    ))

            # maybe draw the fliers
            if showfliers:
                # fliers coords
                flier_x = np.ones(len(stats['fliers'])) * pos
                flier_y = stats['fliers']

                fliers.extend(doplot(
                    flier_x, flier_y, **final_flierprops
                ))

        # fix our axes/ticks up a little
        if vert:
            setticks = self.set_xticks
            setlim = self.set_xlim
            setlabels = self.set_xticklabels
        else:
            setticks = self.set_yticks
            setlim = self.set_ylim
            setlabels = self.set_yticklabels

        if manage_xticks:
            newlimits = min(positions) - 0.5, max(positions) + 0.5
            setlim(newlimits)
            setticks(positions)
            setlabels(datalabels)

        # reset hold status
        self._hold = holdStatus

        return dict(whiskers=whiskers, caps=caps, boxes=boxes,
                    medians=medians, fliers=fliers, means=means)

    @_preprocess_data(replace_names=["x", "y", "s", "linewidths",
                                     "edgecolors", "c", "facecolor",
                                     "facecolors", "color"],
                      label_namer="y")
    def scatter(self, x, y, s=None, c=None, marker=None, cmap=None, norm=None,
                vmin=None, vmax=None, alpha=None, linewidths=None,
                verts=None, edgecolors=None,
                **kwargs):
        """
        A scatter plot of *y* vs *x* with varying marker size and/or color.

        Parameters
        ----------
        x, y : array_like, shape (n, )
            The data positions.

        s : scalar or array_like, shape (n, ), optional
            The marker size in points**2.
            Default is ``rcParams['lines.markersize'] ** 2``.

        c : color, sequence, or sequence of color, optional, default: 'b'
            The marker color. Possible values:

            - A single color format string.
            - A sequence of color specifications of length n.
            - A sequence of n numbers to be mapped to colors using *cmap* and
              *norm*.
            - A 2-D array in which the rows are RGB or RGBA.

            Note that *c* should not be a single numeric RGB or RGBA sequence
            because that is indistinguishable from an array of values to be
            colormapped. If you want to specify the same RGB or RGBA value for
            all points, use a 2-D array with a single row.

        marker : `~matplotlib.markers.MarkerStyle`, optional, default: 'o'
            The marker style. *marker* can be either an instance of the class
            or the text shorthand for a particular marker.
            See `~matplotlib.markers` for more information marker styles.

        cmap : `~matplotlib.colors.Colormap`, optional, default: None
            A `.Colormap` instance or registered colormap name. *cmap* is only
            used if *c* is an array of floats. If ``None``, defaults to rc
            ``image.cmap``.

        norm : `~matplotlib.colors.Normalize`, optional, default: None
            A `.Normalize` instance is used to scale luminance data to 0, 1.
            *norm* is only used if *c* is an array of floats. If *None*, use
            the default `.colors.Normalize`.

        vmin, vmax : scalar, optional, default: None
            *vmin* and *vmax* are used in conjunction with *norm* to normalize
            luminance data. If None, the respective min and max of the color
            array is used. *vmin* and *vmax* are ignored if you pass a *norm*
            instance.

        alpha : scalar, optional, default: None
            The alpha blending value, between 0 (transparent) and 1 (opaque).

        linewidths : scalar or array_like, optional, default: None
            The linewidth of the marker edges. Note: The default *edgecolors*
            is 'face'. You may want to change this as well.
            If *None*, defaults to rcParams ``lines.linewidth``.

        verts : sequence of (x, y), optional
            If *marker* is *None*, these vertices will be used to construct
            the marker.  The center of the marker is located at (0, 0) in
            normalized units.  The overall marker is rescaled by *s*.

        edgecolors : color or sequence of color, optional, default: 'face'
            The edge color of the marker. Possible values:

            - 'face': The edge color will always be the same as the face color.
            - 'none': No patch boundary will be drawn.
            - A matplotib color.

            For non-filled markers, the *edgecolors* kwarg is ignored and
            forced to 'face' internally.

        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.collections.Collection` properties

        See Also
        --------
        plot : To plot scatter plots when markers are identical in size and
            color.

        Notes
        -----

        * The `.plot` function will be faster for scatterplots where markers
          don't vary in size or color.

        * Any or all of *x*, *y*, *s*, and *c* may be masked arrays, in which
          case all masks will be combined and only unmasked points will be
          plotted.

        * Fundamentally, scatter works with 1-D arrays; *x*, *y*, *s*, and *c*
          may be input as 2-D arrays, but within scatter they will be
          flattened. The exception is *c*, which will be flattened only if its
          size matches the size of *x* and *y*.

        """

        if not self._hold:
            self.cla()

        # Process **kwargs to handle aliases, conflicts with explicit kwargs:

        facecolors = None
        edgecolors = kwargs.pop('edgecolor', edgecolors)
        fc = kwargs.pop('facecolors', None)
        fc = kwargs.pop('facecolor', fc)
        if fc is not None:
            facecolors = fc
        co = kwargs.pop('color', None)
        if co is not None:
            try:
                mcolors.to_rgba_array(co)
            except ValueError:
                raise ValueError("'color' kwarg must be an mpl color"
                                 " spec or sequence of color specs.\n"
                                 "For a sequence of values to be"
                                 " color-mapped, use the 'c' kwarg instead.")
            if edgecolors is None:
                edgecolors = co
            if facecolors is None:
                facecolors = co
            if c is not None:
                raise ValueError("Supply a 'c' kwarg or a 'color' kwarg"
                                 " but not both; they differ but"
                                 " their functionalities overlap.")
        if c is None:
            if facecolors is not None:
                c = facecolors
            else:
                if rcParams['_internal.classic_mode']:
                    c = 'b'  # The original default
                else:
                    c = self._get_patches_for_fill.get_next_color()
            c_none = True
        else:
            c_none = False

        if edgecolors is None and not rcParams['_internal.classic_mode']:
            edgecolors = 'face'

        self._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)

        # np.ma.ravel yields an ndarray, not a masked array,
        # unless its argument is a masked array.
        xy_shape = (np.shape(x), np.shape(y))
        x = np.ma.ravel(x)
        y = np.ma.ravel(y)
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        if s is None:
            if rcParams['_internal.classic_mode']:
                s = 20
            else:
                s = rcParams['lines.markersize'] ** 2.0

        s = np.ma.ravel(s)  # This doesn't have to match x, y in size.

        # After this block, c_array will be None unless
        # c is an array for mapping.  The potential ambiguity
        # with a sequence of 3 or 4 numbers is resolved in
        # favor of mapping, not rgb or rgba.
        if c_none or co is not None:
            c_array = None
        else:
            try:
                c_array = np.asanyarray(c, dtype=float)
                if c_array.shape in xy_shape:
                    c = np.ma.ravel(c_array)
                else:
                    # Wrong size; it must not be intended for mapping.
                    c_array = None
            except ValueError:
                # Failed to make a floating-point array; c must be color specs.
                c_array = None

        if c_array is None:
            try:
                # must be acceptable as PathCollection facecolors
                colors = mcolors.to_rgba_array(c)
            except ValueError:
                # c not acceptable as PathCollection facecolor
                raise ValueError("c of shape {} not acceptable as a color "
                                 "sequence for x with size {}, y with size {}"
                                 .format(c.shape, x.size, y.size))
        else:
            colors = None  # use cmap, norm after collection is created

        # `delete_masked_points` only modifies arguments of the same length as
        # `x`.
        x, y, s, c, colors, edgecolors, linewidths =\
            cbook.delete_masked_points(
                x, y, s, c, colors, edgecolors, linewidths)

        scales = s   # Renamed for readability below.

        # to be API compatible
        if marker is None and verts is not None:
            marker = (verts, 0)
            verts = None

        # load default marker from rcParams
        if marker is None:
            marker = rcParams['scatter.marker']

        if isinstance(marker, mmarkers.MarkerStyle):
            marker_obj = marker
        else:
            marker_obj = mmarkers.MarkerStyle(marker)

        path = marker_obj.get_path().transformed(
            marker_obj.get_transform())
        if not marker_obj.is_filled():
            edgecolors = 'face'
            linewidths = rcParams['lines.linewidth']

        offsets = np.column_stack([x, y])

        collection = mcoll.PathCollection(
                (path,), scales,
                facecolors=colors,
                edgecolors=edgecolors,
                linewidths=linewidths,
                offsets=offsets,
                transOffset=kwargs.pop('transform', self.transData),
                alpha=alpha
                )
        collection.set_transform(mtransforms.IdentityTransform())
        collection.update(kwargs)

        if colors is None:
            if norm is not None and not isinstance(norm, mcolors.Normalize):
                raise ValueError(
                    "'norm' must be an instance of 'mcolors.Normalize'")
            collection.set_array(np.asarray(c))
            collection.set_cmap(cmap)
            collection.set_norm(norm)

            if vmin is not None or vmax is not None:
                collection.set_clim(vmin, vmax)
            else:
                collection.autoscale_None()

        # Classic mode only:
        # ensure there are margins to allow for the
        # finite size of the symbols.  In v2.x, margins
        # are present by default, so we disable this
        # scatter-specific override.
        if rcParams['_internal.classic_mode']:
            if self._xmargin < 0.05 and x.size > 0:
                self.set_xmargin(0.05)
            if self._ymargin < 0.05 and x.size > 0:
                self.set_ymargin(0.05)

        self.add_collection(collection)
        self.autoscale_view()

        return collection

    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    @docstring.dedent_interpd
    def hexbin(self, x, y, C=None, gridsize=100, bins=None,
               xscale='linear', yscale='linear', extent=None,
               cmap=None, norm=None, vmin=None, vmax=None,
               alpha=None, linewidths=None, edgecolors='face',
               reduce_C_function=np.mean, mincnt=None, marginals=False,
               **kwargs):
        """
        Make a hexagonal binning plot.

        Make a hexagonal binning plot of *x* versus *y*, where *x*,
        *y* are 1-D sequences of the same length, *N*. If *C* is *None*
        (the default), this is a histogram of the number of occurrences
        of the observations at (x[i],y[i]).

        If *C* is specified, it specifies values at the coordinate
        (x[i],y[i]). These values are accumulated for each hexagonal
        bin and then reduced according to *reduce_C_function*, which
        defaults to numpy's mean function (np.mean). (If *C* is
        specified, it must also be a 1-D sequence of the same length
        as *x* and *y*.)

        Parameters
        ----------
        x, y : array or masked array

        C : array or masked array, optional, default is *None*

        gridsize : int or (int, int), optional, default is 100
            The number of hexagons in the *x*-direction, default is
            100. The corresponding number of hexagons in the
            *y*-direction is chosen such that the hexagons are
            approximately regular. Alternatively, gridsize can be a
            tuple with two elements specifying the number of hexagons
            in the *x*-direction and the *y*-direction.

        bins : {'log'} or int or sequence, optional, default is *None*
            If *None*, no binning is applied; the color of each hexagon
            directly corresponds to its count value.

            If 'log', use a logarithmic scale for the color
            map. Internally, :math:`log_{10}(i+1)` is used to
            determine the hexagon color.

            If an integer, divide the counts in the specified number
            of bins, and color the hexagons accordingly.

            If a sequence of values, the values of the lower bound of
            the bins to be used.

        xscale : {'linear', 'log'}, optional, default is 'linear'
            Use a linear or log10 scale on the horizontal axis.

        yscale : {'linear', 'log'}, optional, default is 'linear'
            Use a linear or log10 scale on the vertical axis.

        mincnt : int > 0, optional, default is *None*
            If not *None*, only display cells with more than *mincnt*
            number of points in the cell

        marginals : bool, optional, default is *False*
            if marginals is *True*, plot the marginal density as
            colormapped rectagles along the bottom of the x-axis and
            left of the y-axis

        extent : scalar, optional, default is *None*
            The limits of the bins. The default assigns the limits
            based on *gridsize*, *x*, *y*, *xscale* and *yscale*.

            If *xscale* or *yscale* is set to 'log', the limits are
            expected to be the exponent for a power of 10. E.g. for
            x-limits of 1 and 50 in 'linear' scale and y-limits
            of 10 and 1000 in 'log' scale, enter (1, 50, 1, 3).

            Order of scalars is (left, right, bottom, top).

        Other Parameters
        ----------------
        cmap : object, optional, default is *None*
            a :class:`matplotlib.colors.Colormap` instance. If *None*,
            defaults to rc ``image.cmap``.

        norm : object, optional, default is *None*
            :class:`matplotlib.colors.Normalize` instance is used to
            scale luminance data to 0,1.

        vmin, vmax : scalar, optional, default is *None*
            *vmin* and *vmax* are used in conjunction with *norm* to
            normalize luminance data. If *None*, the min and max of the
            color array *C* are used.  Note if you pass a norm instance
            your settings for *vmin* and *vmax* will be ignored.

        alpha : scalar between 0 and 1, optional, default is *None*
            the alpha value for the patches

        linewidths : scalar, optional, default is *None*
            If *None*, defaults to 1.0.

        edgecolors : {'face', 'none', *None*} or color, optional

            If 'face' (the default), draws the edges in the same color as the
            fill color.

            If 'none', no edge is drawn; this can sometimes lead to unsightly
            unpainted pixels between the hexagons.

            If *None*, draws outlines in the default color.

            If a matplotlib color arg, draws outlines in the specified color.

        Returns
        -------
        object
            a :class:`~matplotlib.collections.PolyCollection` instance; use
            :meth:`~matplotlib.collections.PolyCollection.get_array` on
            this :class:`~matplotlib.collections.PolyCollection` to get
            the counts in each hexagon.

            If *marginals* is *True*, horizontal
            bar and vertical bar (both PolyCollections) will be attached
            to the return collection as attributes *hbar* and *vbar*.

        Notes
        -----
        The standard descriptions of all the
        :class:`~matplotlib.collections.Collection` parameters:

            %(Collection)s

        """

        if not self._hold:
            self.cla()

        self._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)

        x, y, C = cbook.delete_masked_points(x, y, C)

        # Set the size of the hexagon grid
        if iterable(gridsize):
            nx, ny = gridsize
        else:
            nx = gridsize
            ny = int(nx / math.sqrt(3))
        # Count the number of data in each hexagon
        x = np.array(x, float)
        y = np.array(y, float)
        if xscale == 'log':
            if np.any(x <= 0.0):
                raise ValueError("x contains non-positive values, so can not"
                                 " be log-scaled")
            x = np.log10(x)
        if yscale == 'log':
            if np.any(y <= 0.0):
                raise ValueError("y contains non-positive values, so can not"
                                 " be log-scaled")
            y = np.log10(y)
        if extent is not None:
            xmin, xmax, ymin, ymax = extent
        else:
            xmin, xmax = (np.min(x), np.max(x)) if len(x) else (0, 1)
            ymin, ymax = (np.min(y), np.max(y)) if len(y) else (0, 1)

            # to avoid issues with singular data, expand the min/max pairs
            xmin, xmax = mtransforms.nonsingular(xmin, xmax, expander=0.1)
            ymin, ymax = mtransforms.nonsingular(ymin, ymax, expander=0.1)

        # In the x-direction, the hexagons exactly cover the region from
        # xmin to xmax. Need some padding to avoid roundoff errors.
        padding = 1.e-9 * (xmax - xmin)
        xmin -= padding
        xmax += padding
        sx = (xmax - xmin) / nx
        sy = (ymax - ymin) / ny

        if marginals:
            xorig = x.copy()
            yorig = y.copy()

        x = (x - xmin) / sx
        y = (y - ymin) / sy
        ix1 = np.round(x).astype(int)
        iy1 = np.round(y).astype(int)
        ix2 = np.floor(x).astype(int)
        iy2 = np.floor(y).astype(int)

        nx1 = nx + 1
        ny1 = ny + 1
        nx2 = nx
        ny2 = ny
        n = nx1 * ny1 + nx2 * ny2

        d1 = (x - ix1) ** 2 + 3.0 * (y - iy1) ** 2
        d2 = (x - ix2 - 0.5) ** 2 + 3.0 * (y - iy2 - 0.5) ** 2
        bdist = (d1 < d2)
        if C is None:
            lattice1 = np.zeros((nx1, ny1))
            lattice2 = np.zeros((nx2, ny2))

            cond1 = (0 <= ix1) * (ix1 < nx1) * (0 <= iy1) * (iy1 < ny1)
            cond2 = (0 <= ix2) * (ix2 < nx2) * (0 <= iy2) * (iy2 < ny2)

            cond1 *= bdist
            cond2 *= np.logical_not(bdist)
            ix1, iy1 = ix1[cond1], iy1[cond1]
            ix2, iy2 = ix2[cond2], iy2[cond2]

            for ix, iy in zip(ix1, iy1):
                lattice1[ix, iy] += 1
            for ix, iy in zip(ix2, iy2):
                lattice2[ix, iy] += 1

            # threshold
            if mincnt is not None:
                lattice1[lattice1 < mincnt] = np.nan
                lattice2[lattice2 < mincnt] = np.nan
            accum = np.hstack((lattice1.ravel(),
                               lattice2.ravel()))
            good_idxs = ~np.isnan(accum)

        else:
            if mincnt is None:
                mincnt = 0

            # create accumulation arrays
            lattice1 = np.empty((nx1, ny1), dtype=object)
            for i in xrange(nx1):
                for j in xrange(ny1):
                    lattice1[i, j] = []
            lattice2 = np.empty((nx2, ny2), dtype=object)
            for i in xrange(nx2):
                for j in xrange(ny2):
                    lattice2[i, j] = []

            for i in xrange(len(x)):
                if bdist[i]:
                    if 0 <= ix1[i] < nx1 and 0 <= iy1[i] < ny1:
                        lattice1[ix1[i], iy1[i]].append(C[i])
                else:
                    if 0 <= ix2[i] < nx2 and 0 <= iy2[i] < ny2:
                        lattice2[ix2[i], iy2[i]].append(C[i])

            for i in xrange(nx1):
                for j in xrange(ny1):
                    vals = lattice1[i, j]
                    if len(vals) > mincnt:
                        lattice1[i, j] = reduce_C_function(vals)
                    else:
                        lattice1[i, j] = np.nan
            for i in xrange(nx2):
                for j in xrange(ny2):
                    vals = lattice2[i, j]
                    if len(vals) > mincnt:
                        lattice2[i, j] = reduce_C_function(vals)
                    else:
                        lattice2[i, j] = np.nan

            accum = np.hstack((lattice1.astype(float).ravel(),
                               lattice2.astype(float).ravel()))
            good_idxs = ~np.isnan(accum)

        offsets = np.zeros((n, 2), float)
        offsets[:nx1 * ny1, 0] = np.repeat(np.arange(nx1), ny1)
        offsets[:nx1 * ny1, 1] = np.tile(np.arange(ny1), nx1)
        offsets[nx1 * ny1:, 0] = np.repeat(np.arange(nx2) + 0.5, ny2)
        offsets[nx1 * ny1:, 1] = np.tile(np.arange(ny2), nx2) + 0.5
        offsets[:, 0] *= sx
        offsets[:, 1] *= sy
        offsets[:, 0] += xmin
        offsets[:, 1] += ymin
        # remove accumulation bins with no data
        offsets = offsets[good_idxs, :]
        accum = accum[good_idxs]

        polygon = np.zeros((6, 2), float)
        polygon[:, 0] = sx * np.array([0.5, 0.5, 0.0, -0.5, -0.5, 0.0])
        polygon[:, 1] = sy * np.array([-0.5, 0.5, 1.0, 0.5, -0.5, -1.0]) / 3.0

        if linewidths is None:
            linewidths = [1.0]

        if xscale == 'log' or yscale == 'log':
            polygons = np.expand_dims(polygon, 0) + np.expand_dims(offsets, 1)
            if xscale == 'log':
                polygons[:, :, 0] = 10.0 ** polygons[:, :, 0]
                xmin = 10.0 ** xmin
                xmax = 10.0 ** xmax
                self.set_xscale(xscale)
            if yscale == 'log':
                polygons[:, :, 1] = 10.0 ** polygons[:, :, 1]
                ymin = 10.0 ** ymin
                ymax = 10.0 ** ymax
                self.set_yscale(yscale)
            collection = mcoll.PolyCollection(
                polygons,
                edgecolors=edgecolors,
                linewidths=linewidths,
                )
        else:
            collection = mcoll.PolyCollection(
                [polygon],
                edgecolors=edgecolors,
                linewidths=linewidths,
                offsets=offsets,
                transOffset=mtransforms.IdentityTransform(),
                offset_position="data"
                )

        if isinstance(norm, mcolors.LogNorm):
            if (accum == 0).any():
                # make sure we have not zeros
                accum += 1

        # autoscale the norm with curren accum values if it hasn't
        # been set
        if norm is not None:
            if norm.vmin is None and norm.vmax is None:
                norm.autoscale(accum)

        # Transform accum if needed
        if bins == 'log':
            accum = np.log10(accum + 1)
        elif bins is not None:
            if not iterable(bins):
                minimum, maximum = min(accum), max(accum)
                bins -= 1  # one less edge than bins
                bins = minimum + (maximum - minimum) * np.arange(bins) / bins
            bins = np.sort(bins)
            accum = bins.searchsorted(accum)

        if norm is not None and not isinstance(norm, mcolors.Normalize):
            raise ValueError(
                "'norm' must be an instance of 'mcolors.Normalize'")
        collection.set_array(accum)
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_alpha(alpha)
        collection.update(kwargs)

        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)
        else:
            collection.autoscale_None()

        corners = ((xmin, ymin), (xmax, ymax))
        self.update_datalim(corners)
        collection.sticky_edges.x[:] = [xmin, xmax]
        collection.sticky_edges.y[:] = [ymin, ymax]
        self.autoscale_view(tight=True)

        # add the collection last
        self.add_collection(collection, autolim=False)
        if not marginals:
            return collection

        if C is None:
            C = np.ones(len(x))

        def coarse_bin(x, y, coarse):
            ind = coarse.searchsorted(x).clip(0, len(coarse) - 1)
            mus = np.zeros(len(coarse))
            for i in range(len(coarse)):
                yi = y[ind == i]
                if len(yi) > 0:
                    mu = reduce_C_function(yi)
                else:
                    mu = np.nan
                mus[i] = mu
            return mus

        coarse = np.linspace(xmin, xmax, gridsize)

        xcoarse = coarse_bin(xorig, C, coarse)
        valid = ~np.isnan(xcoarse)
        verts, values = [], []
        for i, val in enumerate(xcoarse):
            thismin = coarse[i]
            if i < len(coarse) - 1:
                thismax = coarse[i + 1]
            else:
                thismax = thismin + np.diff(coarse)[-1]

            if not valid[i]:
                continue

            verts.append([(thismin, 0),
                          (thismin, 0.05),
                          (thismax, 0.05),
                          (thismax, 0)])
            values.append(val)

        values = np.array(values)
        trans = self.get_xaxis_transform(which='grid')

        hbar = mcoll.PolyCollection(verts, transform=trans, edgecolors='face')

        hbar.set_array(values)
        hbar.set_cmap(cmap)
        hbar.set_norm(norm)
        hbar.set_alpha(alpha)
        hbar.update(kwargs)
        self.add_collection(hbar, autolim=False)

        coarse = np.linspace(ymin, ymax, gridsize)
        ycoarse = coarse_bin(yorig, C, coarse)
        valid = ~np.isnan(ycoarse)
        verts, values = [], []
        for i, val in enumerate(ycoarse):
            thismin = coarse[i]
            if i < len(coarse) - 1:
                thismax = coarse[i + 1]
            else:
                thismax = thismin + np.diff(coarse)[-1]
            if not valid[i]:
                continue
            verts.append([(0, thismin), (0.0, thismax),
                          (0.05, thismax), (0.05, thismin)])
            values.append(val)

        values = np.array(values)

        trans = self.get_yaxis_transform(which='grid')

        vbar = mcoll.PolyCollection(verts, transform=trans, edgecolors='face')
        vbar.set_array(values)
        vbar.set_cmap(cmap)
        vbar.set_norm(norm)
        vbar.set_alpha(alpha)
        vbar.update(kwargs)
        self.add_collection(vbar, autolim=False)

        collection.hbar = hbar
        collection.vbar = vbar

        def on_changed(collection):
            hbar.set_cmap(collection.get_cmap())
            hbar.set_clim(collection.get_clim())
            vbar.set_cmap(collection.get_cmap())
            vbar.set_clim(collection.get_clim())

        collection.callbacksSM.connect('changed', on_changed)

        return collection

    @docstring.dedent_interpd
    def arrow(self, x, y, dx, dy, **kwargs):
        """
        Add an arrow to the axes.

        This draws an arrow from ``(x, y)`` to ``(x+dx, y+dy)``.

        Parameters
        ----------
        x, y : float
            The x/y-coordinate of the arrow base.
        dx, dy : float
            The length of the arrow along x/y-direction.

        Returns
        -------
        arrow : `.FancyArrow`
            The created `.FancyArrow` object.

        Other Parameters
        ----------------
        **kwargs
            Optional kwargs (inherited from `.FancyArrow` patch) control the
            arrow construction and properties:

        %(FancyArrow)s

        Notes
        -----
        The resulting arrow is affected by the axes aspect ratio and limits.
        This may produce an arrow whose head is not square with its stem. To
        create an arrow whose head is square with its stem,
        use :meth:`annotate` for example:

        >>> ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0),
        ...             arrowprops=dict(arrowstyle="->"))

        """
        # Strip away units for the underlying patch since units
        # do not make sense to most patch-like code
        x = self.convert_xunits(x)
        y = self.convert_yunits(y)
        dx = self.convert_xunits(dx)
        dy = self.convert_yunits(dy)

        a = mpatches.FancyArrow(x, y, dx, dy, **kwargs)
        self.add_artist(a)
        return a

    def quiverkey(self, *args, **kw):
        qk = mquiver.QuiverKey(*args, **kw)
        self.add_artist(qk)
        return qk
    quiverkey.__doc__ = mquiver.QuiverKey.quiverkey_doc

    # Handle units for x and y, if they've been passed
    def _quiver_units(self, args, kw):
        if len(args) > 3:
            x, y = args[0:2]
            self._process_unit_info(xdata=x, ydata=y, kwargs=kw)
            x = self.convert_xunits(x)
            y = self.convert_yunits(y)
            return (x, y) + args[2:]
        return args

    # args can by a combination if X, Y, U, V, C and all should be replaced
    @_preprocess_data(replace_all_args=True, label_namer=None)
    def quiver(self, *args, **kw):
        if not self._hold:
            self.cla()

        # Make sure units are handled for x and y values
        args = self._quiver_units(args, kw)

        q = mquiver.Quiver(self, *args, **kw)

        self.add_collection(q, autolim=True)
        self.autoscale_view()
        return q
    quiver.__doc__ = mquiver.Quiver.quiver_doc

    # args can by either Y or y1,y2,... and all should be replaced
    @_preprocess_data(replace_all_args=True, label_namer=None)
    def stackplot(self, x, *args, **kwargs):
        return mstack.stackplot(self, x, *args, **kwargs)
    stackplot.__doc__ = mstack.stackplot.__doc__

    @_preprocess_data(replace_names=["x", "y", "u", "v", "start_points"],
                      label_namer=None)
    def streamplot(self, x, y, u, v, density=1, linewidth=None, color=None,
                   cmap=None, norm=None, arrowsize=1, arrowstyle='-|>',
                   minlength=0.1, transform=None, zorder=None,
                   start_points=None, maxlength=4.0,
                   integration_direction='both'):
        if not self._hold:
            self.cla()
        stream_container = mstream.streamplot(
            self, x, y, u, v,
            density=density,
            linewidth=linewidth,
            color=color,
            cmap=cmap,
            norm=norm,
            arrowsize=arrowsize,
            arrowstyle=arrowstyle,
            minlength=minlength,
            start_points=start_points,
            transform=transform,
            zorder=zorder,
            maxlength=maxlength,
            integration_direction=integration_direction)
        return stream_container
    streamplot.__doc__ = mstream.streamplot.__doc__

    # args can be some combination of X, Y, U, V, C and all should be replaced
    @_preprocess_data(replace_all_args=True, label_namer=None)
    @docstring.dedent_interpd
    def barbs(self, *args, **kw):
        """
        %(barbs_doc)s
        """
        if not self._hold:
            self.cla()

        # Make sure units are handled for x and y values
        args = self._quiver_units(args, kw)

        b = mquiver.Barbs(self, *args, **kw)
        self.add_collection(b, autolim=True)
        self.autoscale_view()
        return b

    @_preprocess_data(replace_names=["x", "y"], label_namer=None,
                      positional_parameter_names=["x", "y", "c"])
    def fill(self, *args, **kwargs):
        """
        Plot filled polygons.

        Parameters
        ----------
        args : sequence of x, y, [color]
            Each polygon is defined by the lists of *x* and *y* positions of
            its nodes, optionally followed by by a *color* specifier. See
            :mod:`matplotlib.colors` for supported color specifiers. The
            standard color cycle is used for polygons without a color
            specifier.

            You can plot multiple polygons by providing multiple *x*, *y*,
            *[color]* groups.

            For example, each of the following is legal::

                ax.fill(x, y)                    # a polygon with default color
                ax.fill(x, y, "b")               # a blue polygon
                ax.fill(x, y, x2, y2)            # two polygons
                ax.fill(x, y, "b", x2, y2, "r")  # a blue and a red polygon

        Returns
        -------
        a list of :class:`~matplotlib.patches.Polygon`

        Other Parameters
        ----------------
        **kwargs : :class:`~matplotlib.patches.Polygon` properties

        Notes
        -----
        Use :meth:`fill_between` if you would like to fill the region between
        two curves.
        """
        if not self._hold:
            self.cla()

        kwargs = cbook.normalize_kwargs(kwargs, _alias_map)

        patches = []
        for poly in self._get_patches_for_fill(*args, **kwargs):
            self.add_patch(poly)
            patches.append(poly)
        self.autoscale_view()
        return patches

    @_preprocess_data(replace_names=["x", "y1", "y2", "where"],
                      label_namer=None)
    @docstring.dedent_interpd
    def fill_between(self, x, y1, y2=0, where=None, interpolate=False,
                     step=None, **kwargs):
        """
        Fill the area between two horizontal curves.

        The curves are defined by the points (*x*, *y1*) and (*x*, *y2*). This
        creates one or multiple polygons describing the filled area.

        You may exclude some horizontal sections from filling using *where*.

        By default, the edges connect the given points directly. Use *step* if
        the filling should be a step function, i.e. constant in between *x*.


        Parameters
        ----------
        x : array (length N)
            The x coordinates of the nodes defining the curves.

        y1 : array (length N) or scalar
            The y coordinates of the nodes defining the first curve.

        y2 : array (length N) or scalar, optional, default: 0
            The y coordinates of the nodes defining the second curve.

        where : array of bool (length N), optional, default: None
            Define *where* to exclude some horizontal regions from being
            filled. The filled regions are defined by the coordinates
            ``x[where]``.  More precisely, fill between ``x[i]`` and ``x[i+1]``
            if ``where[i] and where[i+1]``.  Note that this definition implies
            that an isolated *True* value between two *False* values in
            *where* will not result in filling.  Both sides of the *True*
            position remain unfilled due to the adjacent *False* values.

        interpolate : bool, optional
            This option is only relvant if *where* is used and the two curves
            are crossing each other.

            Semantically, *where* is often used for *y1* > *y2* or similar.
            By default, the nodes of the polygon defining the filled region
            will only be placed at the positions in the *x* array.  Such a
            polygon cannot describe the above semantics close to the
            intersection.  The x-sections containing the intersecion are
            simply clipped.

            Setting *interpolate* to *True* will calculate the actual
            interscection point and extend the filled region up to this point.

        step : {'pre', 'post', 'mid'}, optional
            Define *step* if the filling should be a step function,
            i.e. constant in between *x*. The value determines where the
            step will occur:

            - 'pre': The y value is continued constantly to the left from
              every *x* position, i.e. the interval ``(x[i-1], x[i]]`` has the
              value ``y[i]``.
            - 'post': The y value is continued constantly to the right from
              every *x* position, i.e. the interval ``[x[i], x[i+1])`` has the
              value ``y[i]``.
            - 'mid': Steps occur half-way between the *x* positions.

        Other Parameters
        ----------------
        **kwargs
            All other keyword arguments are passed on to `.PolyCollection`.
            They control the `.Polygon` properties:

            %(PolyCollection)s

        Returns
        -------
        `.PolyCollection`
            A `.PolyCollection` containing the plotted polygons.

        See Also
        --------
        fill_betweenx : Fill between two sets of x-values.

        Notes
        -----
        .. [notes section required to get data note injection right]

        """
        if not rcParams['_internal.classic_mode']:
            color_aliases = mcoll._color_aliases
            kwargs = cbook.normalize_kwargs(kwargs, color_aliases)

            if not any(c in kwargs for c in ('color', 'facecolors')):
                fc = self._get_patches_for_fill.get_next_color()
                kwargs['facecolors'] = fc

        # Handle united data, such as dates
        self._process_unit_info(xdata=x, ydata=y1, kwargs=kwargs)
        self._process_unit_info(ydata=y2)

        # Convert the arrays so we can work with them
        x = ma.masked_invalid(self.convert_xunits(x))
        y1 = ma.masked_invalid(self.convert_yunits(y1))
        y2 = ma.masked_invalid(self.convert_yunits(y2))

        for name, array in [('x', x), ('y1', y1), ('y2', y2)]:
            if array.ndim > 1:
                raise ValueError('Input passed into argument "%r"' % name +
                                 'is not 1-dimensional.')

        if where is None:
            where = True
        where = where & ~functools.reduce(np.logical_or,
                                          map(np.ma.getmask, [x, y1, y2]))

        x, y1, y2 = np.broadcast_arrays(np.atleast_1d(x), y1, y2)

        polys = []
        for ind0, ind1 in cbook.contiguous_regions(where):
            xslice = x[ind0:ind1]
            y1slice = y1[ind0:ind1]
            y2slice = y2[ind0:ind1]
            if step is not None:
                step_func = STEP_LOOKUP_MAP["steps-" + step]
                xslice, y1slice, y2slice = step_func(xslice, y1slice, y2slice)

            if not len(xslice):
                continue

            N = len(xslice)
            X = np.zeros((2 * N + 2, 2), float)

            if interpolate:
                def get_interp_point(ind):
                    im1 = max(ind - 1, 0)
                    x_values = x[im1:ind + 1]
                    diff_values = y1[im1:ind + 1] - y2[im1:ind + 1]
                    y1_values = y1[im1:ind + 1]

                    if len(diff_values) == 2:
                        if np.ma.is_masked(diff_values[1]):
                            return x[im1], y1[im1]
                        elif np.ma.is_masked(diff_values[0]):
                            return x[ind], y1[ind]

                    diff_order = diff_values.argsort()
                    diff_root_x = np.interp(
                        0, diff_values[diff_order], x_values[diff_order])
                    x_order = x_values.argsort()
                    diff_root_y = np.interp(diff_root_x, x_values[x_order],
                                            y1_values[x_order])
                    return diff_root_x, diff_root_y

                start = get_interp_point(ind0)
                end = get_interp_point(ind1)
            else:
                # the purpose of the next two lines is for when y2 is a
                # scalar like 0 and we want the fill to go all the way
                # down to 0 even if none of the y1 sample points do
                start = xslice[0], y2slice[0]
                end = xslice[-1], y2slice[-1]

            X[0] = start
            X[N + 1] = end

            X[1:N + 1, 0] = xslice
            X[1:N + 1, 1] = y1slice
            X[N + 2:, 0] = xslice[::-1]
            X[N + 2:, 1] = y2slice[::-1]

            polys.append(X)

        collection = mcoll.PolyCollection(polys, **kwargs)

        # now update the datalim and autoscale
        XY1 = np.array([x[where], y1[where]]).T
        XY2 = np.array([x[where], y2[where]]).T
        self.dataLim.update_from_data_xy(XY1, self.ignore_existing_data_limits,
                                         updatex=True, updatey=True)
        self.ignore_existing_data_limits = False
        self.dataLim.update_from_data_xy(XY2, self.ignore_existing_data_limits,
                                         updatex=False, updatey=True)
        self.add_collection(collection, autolim=False)
        self.autoscale_view()
        return collection

    @_preprocess_data(replace_names=["y", "x1", "x2", "where"],
                      label_namer=None)
    @docstring.dedent_interpd
    def fill_betweenx(self, y, x1, x2=0, where=None,
                      step=None, interpolate=False, **kwargs):
        """
        Fill the area between two vertical curves.

        The curves are defined by the points (*x1*, *y*) and (*x2*, *y*). This
        creates one or multiple polygons describing the filled area.

        You may exclude some vertical sections from filling using *where*.

        By default, the edges connect the given points directly. Use *step* if
        the filling should be a step function, i.e. constant in between *y*.


        Parameters
        ----------
        y : array (length N)
            The y coordinates of the nodes defining the curves.

        x1 : array (length N) or scalar
            The x coordinates of the nodes defining the first curve.

        x2 : array (length N) or scalar, optional, default: 0
            The x coordinates of the nodes defining the second curve.

        where : array of bool (length N), optional, default: None
            Define *where* to exclude some vertical regions from being
            filled. The filled regions are defined by the coordinates
            ``y[where]``.  More precisely, fill between ``y[i]`` and ``y[i+1]``
            if ``where[i] and where[i+1]``.  Note that this definition implies
            that an isolated *True* value between two *False* values in
            *where* will not result in filling.  Both sides of the *True*
            position remain unfilled due to the adjacent *False* values.

        interpolate : bool, optional
            This option is only relvant if *where* is used and the two curves
            are crossing each other.

            Semantically, *where* is often used for *x1* > *x2* or similar.
            By default, the nodes of the polygon defining the filled region
            will only be placed at the positions in the *y* array.  Such a
            polygon cannot describe the above semantics close to the
            intersection.  The y-sections containing the intersecion are
            simply clipped.

            Setting *interpolate* to *True* will calculate the actual
            interscection point and extend the filled region up to this point.

        step : {'pre', 'post', 'mid'}, optional
            Define *step* if the filling should be a step function,
            i.e. constant in between *y*. The value determines where the
            step will occur:

            - 'pre': The y value is continued constantly to the left from
              every *x* position, i.e. the interval ``(x[i-1], x[i]]`` has the
              value ``y[i]``.
            - 'post': The y value is continued constantly to the right from
              every *x* position, i.e. the interval ``[x[i], x[i+1])`` has the
              value ``y[i]``.
            - 'mid': Steps occur half-way between the *x* positions.

        Other Parameters
        ----------------
        **kwargs
            All other keyword arguments are passed on to `.PolyCollection`.
            They control the `.Polygon` properties:

            %(PolyCollection)s

        Returns
        -------
        `.PolyCollection`
            A `.PolyCollection` containing the plotted polygons.

        See Also
        --------
        fill_between : Fill between two sets of y-values.

        Notes
        -----
        .. [notes section required to get data note injection right]

        """
        if not rcParams['_internal.classic_mode']:
            color_aliases = mcoll._color_aliases
            kwargs = cbook.normalize_kwargs(kwargs, color_aliases)

            if not any(c in kwargs for c in ('color', 'facecolors')):
                fc = self._get_patches_for_fill.get_next_color()
                kwargs['facecolors'] = fc
        # Handle united data, such as dates
        self._process_unit_info(ydata=y, xdata=x1, kwargs=kwargs)
        self._process_unit_info(xdata=x2)

        # Convert the arrays so we can work with them
        y = ma.masked_invalid(self.convert_yunits(y))
        x1 = ma.masked_invalid(self.convert_xunits(x1))
        x2 = ma.masked_invalid(self.convert_xunits(x2))

        for name, array in [('y', y), ('x1', x1), ('x2', x2)]:
            if array.ndim > 1:
                raise ValueError('Input passed into argument "%r"' % name +
                                 'is not 1-dimensional.')

        if where is None:
            where = True
        where = where & ~functools.reduce(np.logical_or,
                                          map(np.ma.getmask, [y, x1, x2]))

        y, x1, x2 = np.broadcast_arrays(np.atleast_1d(y), x1, x2)

        polys = []
        for ind0, ind1 in cbook.contiguous_regions(where):
            yslice = y[ind0:ind1]
            x1slice = x1[ind0:ind1]
            x2slice = x2[ind0:ind1]
            if step is not None:
                step_func = STEP_LOOKUP_MAP["steps-" + step]
                yslice, x1slice, x2slice = step_func(yslice, x1slice, x2slice)

            if not len(yslice):
                continue

            N = len(yslice)
            Y = np.zeros((2 * N + 2, 2), float)
            if interpolate:
                def get_interp_point(ind):
                    im1 = max(ind - 1, 0)
                    y_values = y[im1:ind + 1]
                    diff_values = x1[im1:ind + 1] - x2[im1:ind + 1]
                    x1_values = x1[im1:ind + 1]

                    if len(diff_values) == 2:
                        if np.ma.is_masked(diff_values[1]):
                            return x1[im1], y[im1]
                        elif np.ma.is_masked(diff_values[0]):
                            return x1[ind], y[ind]

                    diff_order = diff_values.argsort()
                    diff_root_y = np.interp(
                        0, diff_values[diff_order], y_values[diff_order])
                    y_order = y_values.argsort()
                    diff_root_x = np.interp(diff_root_y, y_values[y_order],
                                            x1_values[y_order])
                    return diff_root_x, diff_root_y

                start = get_interp_point(ind0)
                end = get_interp_point(ind1)
            else:
                # the purpose of the next two lines is for when x2 is a
                # scalar like 0 and we want the fill to go all the way
                # down to 0 even if none of the x1 sample points do
                start = x2slice[0], yslice[0]
                end = x2slice[-1], yslice[-1]

            Y[0] = start
            Y[N + 1] = end

            Y[1:N + 1, 0] = x1slice
            Y[1:N + 1, 1] = yslice
            Y[N + 2:, 0] = x2slice[::-1]
            Y[N + 2:, 1] = yslice[::-1]

            polys.append(Y)

        collection = mcoll.PolyCollection(polys, **kwargs)

        # now update the datalim and autoscale
        X1Y = np.array([x1[where], y[where]]).T
        X2Y = np.array([x2[where], y[where]]).T
        self.dataLim.update_from_data_xy(X1Y, self.ignore_existing_data_limits,
                                         updatex=True, updatey=True)
        self.ignore_existing_data_limits = False
        self.dataLim.update_from_data_xy(X2Y, self.ignore_existing_data_limits,
                                         updatex=True, updatey=False)
        self.add_collection(collection, autolim=False)
        self.autoscale_view()
        return collection

    #### plotting z(x,y): imshow, pcolor and relatives, contour
    @_preprocess_data(label_namer=None)
    def imshow(self, X, cmap=None, norm=None, aspect=None,
               interpolation=None, alpha=None, vmin=None, vmax=None,
               origin=None, extent=None, shape=None, filternorm=1,
               filterrad=4.0, imlim=None, resample=None, url=None, **kwargs):
        """
        Display an image on the axes.

        Parameters
        ----------
        X : array_like, shape (n, m) or (n, m, 3) or (n, m, 4)
            Display the image in `X` to current axes.  `X` may be an
            array or a PIL image. If `X` is an array, it
            can have the following shapes and types:

            - MxN -- values to be mapped (float or int)
            - MxNx3 -- RGB (float or uint8)
            - MxNx4 -- RGBA (float or uint8)

            MxN arrays are mapped to colors based on the `norm` (mapping
            scalar to scalar) and the `cmap` (mapping the normed scalar to
            a color).

            Elements of RGB and RGBA arrays represent pixels of an MxN image.
            All values should be in the range [0 .. 1] for floats or
            [0 .. 255] for integers.  Out-of-range values will be clipped to
            these bounds.

        cmap : `~matplotlib.colors.Colormap`, optional, default: None
            If None, default to rc `image.cmap` value. `cmap` is ignored
            if `X` is 3-D, directly specifying RGB(A) values.

        aspect : ['auto' | 'equal' | scalar], optional, default: None
            If 'auto', changes the image aspect ratio to match that of the
            axes.

            If 'equal', and `extent` is None, changes the axes aspect ratio to
            match that of the image. If `extent` is not `None`, the axes
            aspect ratio is changed to match that of the extent.

            If None, default to rc ``image.aspect`` value.

        interpolation : string, optional, default: None
            Acceptable values are 'none', 'nearest', 'bilinear', 'bicubic',
            'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser',
            'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc',
            'lanczos'

            If `interpolation` is None, default to rc `image.interpolation`.
            See also the `filternorm` and `filterrad` parameters.
            If `interpolation` is 'none', then no interpolation is performed
            on the Agg, ps and pdf backends. Other backends will fall back to
            'nearest'.

        norm : `~matplotlib.colors.Normalize`, optional, default: None
            A `~matplotlib.colors.Normalize` instance is used to scale
            a 2-D float `X` input to the (0, 1) range for input to the
            `cmap`. If `norm` is None, use the default func:`normalize`.
            If `norm` is an instance of `~matplotlib.colors.NoNorm`,
            `X` must be an array of integers that index directly into
            the lookup table of the `cmap`.

        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with norm to normalize
            luminance data.  Note if you pass a `norm` instance, your
            settings for `vmin` and `vmax` will be ignored.

        alpha : scalar, optional, default: None
            The alpha blending value, between 0 (transparent) and 1 (opaque).
            The ``alpha`` argument is ignored for RGBA input data.

        origin : ['upper' | 'lower'], optional, default: None
            Place the [0,0] index of the array in the upper left or lower left
            corner of the axes. If None, default to rc `image.origin`.

        extent : scalars (left, right, bottom, top), optional, default: None
            The location, in data-coordinates, of the lower-left and
            upper-right corners. If `None`, the image is positioned such that
            the pixel centers fall on zero-based (row, column) indices.

        shape : scalars (columns, rows), optional, default: None
            For raw buffer images

        filternorm : scalar, optional, default: 1
            A parameter for the antigrain image resize filter.  From the
            antigrain documentation, if `filternorm` = 1, the filter
            normalizes integer values and corrects the rounding errors. It
            doesn't do anything with the source floating point values, it
            corrects only integers according to the rule of 1.0 which means
            that any sum of pixel weights must be equal to 1.0.  So, the
            filter function must produce a graph of the proper shape.

        filterrad : scalar, optional, default: 4.0
            The filter radius for filters that have a radius parameter, i.e.
            when interpolation is one of: 'sinc', 'lanczos' or 'blackman'

        Returns
        -------
        image : `~matplotlib.image.AxesImage`

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.artist.Artist` properties.

        See also
        --------
        matshow : Plot a matrix or an array as an image.

        Notes
        -----
        Unless *extent* is used, pixel centers will be located at integer
        coordinates. In other words: the origin will coincide with the center
        of pixel (0, 0).

        Two typical representations are used for RGB images with an alpha
        channel:

        -   Straight (unassociated) alpha: R, G, and B channels represent the
            color of the pixel, disregarding its opacity.
        -   Premultiplied (associated) alpha: R, G, and B channels represent
            the color of the pixel, adjusted for its opacity by multiplication.

        `~matplotlib.pyplot.imshow` expects RGB images adopting the straight
        (unassociated) alpha representation.
        """

        if not self._hold:
            self.cla()

        if norm is not None and not isinstance(norm, mcolors.Normalize):
            raise ValueError(
                "'norm' must be an instance of 'mcolors.Normalize'")
        if aspect is None:
            aspect = rcParams['image.aspect']
        self.set_aspect(aspect)
        im = mimage.AxesImage(self, cmap, norm, interpolation, origin, extent,
                              filternorm=filternorm, filterrad=filterrad,
                              resample=resample, **kwargs)

        im.set_data(X)
        im.set_alpha(alpha)
        if im.get_clip_path() is None:
            # image does not already have clipping set, clip to axes patch
            im.set_clip_path(self.patch)
        #if norm is None and shape is None:
        #    im.set_clim(vmin, vmax)
        if vmin is not None or vmax is not None:
            im.set_clim(vmin, vmax)
        else:
            im.autoscale_None()
        im.set_url(url)

        # update ax.dataLim, and, if autoscaling, set viewLim
        # to tightly fit the image, regardless of dataLim.
        im.set_extent(im.get_extent())

        self.add_image(im)
        return im

    @staticmethod
    def _pcolorargs(funcname, *args, **kw):
        # This takes one kwarg, allmatch.
        # If allmatch is True, then the incoming X, Y, C must
        # have matching dimensions, taking into account that
        # X and Y can be 1-D rather than 2-D.  This perfect
        # match is required for Gouroud shading.  For flat
        # shading, X and Y specify boundaries, so we need
        # one more boundary than color in each direction.
        # For convenience, and consistent with Matlab, we
        # discard the last row and/or column of C if necessary
        # to meet this condition.  This is done if allmatch
        # is False.

        allmatch = kw.pop("allmatch", False)

        if len(args) == 1:
            C = np.asanyarray(args[0])
            numRows, numCols = C.shape
            if allmatch:
                X, Y = np.meshgrid(np.arange(numCols), np.arange(numRows))
            else:
                X, Y = np.meshgrid(np.arange(numCols + 1),
                                   np.arange(numRows + 1))
            C = cbook.safe_masked_invalid(C)
            return X, Y, C

        if len(args) == 3:
            # Check x and y for bad data...
            C = np.asanyarray(args[2])
            X, Y = [cbook.safe_masked_invalid(a) for a in args[:2]]
            if funcname == 'pcolormesh':
                if np.ma.is_masked(X) or np.ma.is_masked(Y):
                    raise ValueError(
                        'x and y arguments to pcolormesh cannot have '
                        'non-finite values or be of type '
                        'numpy.ma.core.MaskedArray with masked values')
                # safe_masked_invalid() returns an ndarray for dtypes other
                # than floating point.
                if isinstance(X, np.ma.core.MaskedArray):
                    X = X.data  # strip mask as downstream doesn't like it...
                if isinstance(Y, np.ma.core.MaskedArray):
                    Y = Y.data
            numRows, numCols = C.shape
        else:
            raise TypeError(
                'Illegal arguments to %s; see help(%s)' % (funcname, funcname))

        Nx = X.shape[-1]
        Ny = Y.shape[0]
        if X.ndim != 2 or X.shape[0] == 1:
            x = X.reshape(1, Nx)
            X = x.repeat(Ny, axis=0)
        if Y.ndim != 2 or Y.shape[1] == 1:
            y = Y.reshape(Ny, 1)
            Y = y.repeat(Nx, axis=1)
        if X.shape != Y.shape:
            raise TypeError(
                'Incompatible X, Y inputs to %s; see help(%s)' % (
                funcname, funcname))
        if allmatch:
            if not (Nx == numCols and Ny == numRows):
                raise TypeError('Dimensions of C %s are incompatible with'
                                ' X (%d) and/or Y (%d); see help(%s)' % (
                                    C.shape, Nx, Ny, funcname))
        else:
            if not (numCols in (Nx, Nx - 1) and numRows in (Ny, Ny - 1)):
                raise TypeError('Dimensions of C %s are incompatible with'
                                ' X (%d) and/or Y (%d); see help(%s)' % (
                                    C.shape, Nx, Ny, funcname))
            C = C[:Ny - 1, :Nx - 1]
        C = cbook.safe_masked_invalid(C)
        return X, Y, C

    @_preprocess_data(label_namer=None)
    @docstring.dedent_interpd
    def pcolor(self, *args, **kwargs):
        """
        Create a pseudocolor plot of a 2-D array.

        Call signatures::

            pcolor(C, **kwargs)
            pcolor(X, Y, C, **kwargs)

        pcolor can be very slow for large arrays; consider
        using the similar but much faster
        :func:`~matplotlib.pyplot.pcolormesh` instead.

        Parameters
        ----------
        C : array_like
            An array of color values.

        X, Y : array_like, optional
            If given, specify the (x, y) coordinates of the colored
            quadrilaterals; the quadrilateral for ``C[i,j]`` has corners at::

                (X[i,   j],   Y[i,   j]),
                (X[i,   j+1], Y[i,   j+1]),
                (X[i+1, j],   Y[i+1, j]),
                (X[i+1, j+1], Y[i+1, j+1])

            Ideally the dimensions of ``X`` and ``Y`` should be one greater
            than those of ``C``; if the dimensions are the same, then the last
            row and column of ``C`` will be ignored.

            Note that the column index corresponds to the
            x-coordinate, and the row index corresponds to y; for
            details, see the :ref:`Grid Orientation
            <axes-pcolor-grid-orientation>` section below.

            If either or both of ``X`` and ``Y`` are 1-D arrays or column
            vectors, they will be expanded as needed into the appropriate 2-D
            arrays, making a rectangular grid.

        cmap : `~matplotlib.colors.Colormap`, optional, default: None
            If `None`, default to rc settings.

        norm : `matplotlib.colors.Normalize`, optional, default: None
            An instance is used to scale luminance data to (0, 1).
            If `None`, defaults to :func:`normalize`.

        vmin, vmax : scalar, optional, default: None
            ``vmin`` and ``vmax`` are used in conjunction with ``norm`` to
            normalize luminance data.  If either is `None`, it is autoscaled to
            the respective min or max of the color array ``C``.  If not `None`,
            ``vmin`` or ``vmax`` passed in here override any pre-existing
            values supplied in the ``norm`` instance.

        edgecolors : {None, 'none', color, color sequence}
            If None, the rc setting is used by default.
            If 'none', edges will not be visible.
            An mpl color or sequence of colors will set the edge color.

        alpha : scalar, optional, default: None
            The alpha blending value, between 0 (transparent) and 1 (opaque).

        snap : bool, optional, default: False
            Whether to snap the mesh to pixel boundaries.

        Returns
        -------
        collection : `matplotlib.collections.Collection`

        Other Parameters
        ----------------
        antialiaseds : bool, optional, default: False
            The default ``antialiaseds`` is False if the default
            ``edgecolors="none"`` is used.  This eliminates artificial lines
            at patch boundaries, and works regardless of the value of alpha.
            If ``edgecolors`` is not "none", then the default ``antialiaseds``
            is taken from :rc:`patch.antialiased`, which defaults to True.
            Stroking the edges may be preferred if ``alpha`` is 1, but will
            cause artifacts otherwise.

        **kwargs :

            Any unused keyword arguments are passed along to the
            `~matplotlib.collections.PolyCollection` constructor:

        %(PolyCollection)s

        See Also
        --------
        pcolormesh : for an explanation of the differences between
            pcolor and pcolormesh.

        Notes
        -----
        .. _axes-pcolor-grid-orientation:

        ``X``, ``Y`` and ``C`` may be masked arrays. If either C[i, j], or one
        of the vertices surrounding C[i,j] (``X`` or ``Y`` at [i, j], [i+1, j],
        [i, j+1], [i+1, j+1]) is masked, nothing is plotted.

        The grid orientation follows the MATLAB convention: an array ``C`` with
        shape (nrows, ncolumns) is plotted with the column number as ``X`` and
        the row number as ``Y``, increasing up; hence it is plotted the way the
        array would be printed, except that the ``Y`` axis is reversed. That
        is, ``C`` is taken as ``C`` (y, x).

        Similarly for :func:`meshgrid`::

            x = np.arange(5)
            y = np.arange(3)
            X, Y = np.meshgrid(x, y)

        is equivalent to::

            X = array([[0, 1, 2, 3, 4],
                       [0, 1, 2, 3, 4],
                       [0, 1, 2, 3, 4]])

            Y = array([[0, 0, 0, 0, 0],
                       [1, 1, 1, 1, 1],
                       [2, 2, 2, 2, 2]])

        so if you have::

            C = rand(len(x), len(y))

        then you need to transpose C::

            pcolor(X, Y, C.T)

        or::

            pcolor(C.T)

        MATLAB :func:`pcolor` always discards the last row and column of ``C``,
        but Matplotlib displays the last row and column if ``X`` and ``Y`` are
        not specified, or if ``X`` and ``Y`` have one more row and column than
        ``C``.
        """

        if not self._hold:
            self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)

        X, Y, C = self._pcolorargs('pcolor', *args, allmatch=False)
        Ny, Nx = X.shape

        # unit conversion allows e.g. datetime objects as axis values
        self._process_unit_info(xdata=X, ydata=Y, kwargs=kwargs)
        X = self.convert_xunits(X)
        Y = self.convert_yunits(Y)

        # convert to MA, if necessary.
        C = ma.asarray(C)
        X = ma.asarray(X)
        Y = ma.asarray(Y)

        mask = ma.getmaskarray(X) + ma.getmaskarray(Y)
        xymask = (mask[0:-1, 0:-1] + mask[1:, 1:] +
                  mask[0:-1, 1:] + mask[1:, 0:-1])
        # don't plot if C or any of the surrounding vertices are masked.
        mask = ma.getmaskarray(C) + xymask

        newaxis = np.newaxis
        compress = np.compress

        ravelmask = (mask == 0).ravel()
        X1 = compress(ravelmask, ma.filled(X[0:-1, 0:-1]).ravel())
        Y1 = compress(ravelmask, ma.filled(Y[0:-1, 0:-1]).ravel())
        X2 = compress(ravelmask, ma.filled(X[1:, 0:-1]).ravel())
        Y2 = compress(ravelmask, ma.filled(Y[1:, 0:-1]).ravel())
        X3 = compress(ravelmask, ma.filled(X[1:, 1:]).ravel())
        Y3 = compress(ravelmask, ma.filled(Y[1:, 1:]).ravel())
        X4 = compress(ravelmask, ma.filled(X[0:-1, 1:]).ravel())
        Y4 = compress(ravelmask, ma.filled(Y[0:-1, 1:]).ravel())
        npoly = len(X1)

        xy = np.concatenate((X1[:, newaxis], Y1[:, newaxis],
                             X2[:, newaxis], Y2[:, newaxis],
                             X3[:, newaxis], Y3[:, newaxis],
                             X4[:, newaxis], Y4[:, newaxis],
                             X1[:, newaxis], Y1[:, newaxis]),
                            axis=1)
        verts = xy.reshape((npoly, 5, 2))

        C = compress(ravelmask, ma.filled(C[0:Ny - 1, 0:Nx - 1]).ravel())

        linewidths = (0.25,)
        if 'linewidth' in kwargs:
            kwargs['linewidths'] = kwargs.pop('linewidth')
        kwargs.setdefault('linewidths', linewidths)

        if 'edgecolor' in kwargs:
            kwargs['edgecolors'] = kwargs.pop('edgecolor')
        ec = kwargs.setdefault('edgecolors', 'none')

        # aa setting will default via collections to patch.antialiased
        # unless the boundary is not stroked, in which case the
        # default will be False; with unstroked boundaries, aa
        # makes artifacts that are often disturbing.
        if 'antialiased' in kwargs:
            kwargs['antialiaseds'] = kwargs.pop('antialiased')
        if 'antialiaseds' not in kwargs and (
                isinstance(ec, six.string_types) and ec.lower() == "none"):
            kwargs['antialiaseds'] = False

        kwargs.setdefault('snap', False)

        collection = mcoll.PolyCollection(verts, **kwargs)

        collection.set_alpha(alpha)
        collection.set_array(C)
        if norm is not None and not isinstance(norm, mcolors.Normalize):
            raise ValueError(
                "'norm' must be an instance of 'mcolors.Normalize'")
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_clim(vmin, vmax)
        collection.autoscale_None()
        self.grid(False)

        x = X.compressed()
        y = Y.compressed()

        # Transform from native to data coordinates?
        t = collection._transform
        if (not isinstance(t, mtransforms.Transform) and
            hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)

        if t and any(t.contains_branch_seperately(self.transData)):
            trans_to_data = t - self.transData
            pts = np.vstack([x, y]).T.astype(float)
            transformed_pts = trans_to_data.transform(pts)
            x = transformed_pts[..., 0]
            y = transformed_pts[..., 1]

        self.add_collection(collection, autolim=False)

        minx = np.min(x)
        maxx = np.max(x)
        miny = np.min(y)
        maxy = np.max(y)
        collection.sticky_edges.x[:] = [minx, maxx]
        collection.sticky_edges.y[:] = [miny, maxy]
        corners = (minx, miny), (maxx, maxy)
        self.update_datalim(corners)
        self.autoscale_view()
        return collection

    @_preprocess_data(label_namer=None)
    @docstring.dedent_interpd
    def pcolormesh(self, *args, **kwargs):
        """
        Plot a quadrilateral mesh.

        Call signatures::

          pcolormesh(C)
          pcolormesh(X, Y, C)
          pcolormesh(C, **kwargs)

        Create a pseudocolor plot of a 2-D array.

        pcolormesh is similar to :func:`~matplotlib.pyplot.pcolor`,
        but uses a different mechanism and returns a different
        object; pcolor returns a
        :class:`~matplotlib.collections.PolyCollection` but pcolormesh
        returns a
        :class:`~matplotlib.collections.QuadMesh`.  It is much faster,
        so it is almost always preferred for large arrays.

        *C* may be a masked array, but *X* and *Y* may not.  Masked
        array support is implemented via *cmap* and *norm*; in
        contrast, :func:`~matplotlib.pyplot.pcolor` simply does not
        draw quadrilaterals with masked colors or vertices.

        Other Parameters
        ----------------
        cmap : Colormap, optional
            A :class:`matplotlib.colors.Colormap` instance. If ``None``, use
            rc settings.

        norm : Normalize, optional
            A :class:`matplotlib.colors.Normalize` instance is used to
            scale luminance data to 0,1. If ``None``, defaults to
            :func:`normalize`.

        vmin, vmax : scalar, optional
            *vmin* and *vmax* are used in conjunction with *norm* to
            normalize luminance data. If either is ``None``, it is autoscaled
            to the respective min or max of the color array *C*.
            If not ``None``, *vmin* or *vmax* passed in here override any
            pre-existing values supplied in the *norm* instance.

        shading : [ 'flat' | 'gouraud' ], optional
            'flat' indicates a solid color for each quad. When
            'gouraud', each quad will be Gouraud shaded. When gouraud
            shading, *edgecolors* is ignored.

        edgecolors : string, color, color sequence, optional
            - If ``None``, the rc setting is used by default.
            - If ``'None'``, edges will not be visible.
            - If ``'face'``, edges will have the same color as the faces.

            An mpl color or sequence of colors will also set the edge color.

        alpha : scalar, optional
            Alpha blending value. Must be between 0 and 1.

        Returns
        -------
        matplotlib.collections.QuadMesh

        See Also
        --------
        matplotlib.pyplot.pcolor :
            For an explanation of the grid orientation
            (:ref:`Grid Orientation <axes-pcolor-grid-orientation>`)
            and the expansion of 1-D *X* and/or *Y* to 2-D arrays.

        Notes
        -----
        kwargs can be used to control the
        :class:`matplotlib.collections.QuadMesh` properties:

        %(QuadMesh)s
        """
        if not self._hold:
            self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        shading = kwargs.pop('shading', 'flat').lower()
        antialiased = kwargs.pop('antialiased', False)
        kwargs.setdefault('edgecolors', 'None')

        allmatch = (shading == 'gouraud')

        X, Y, C = self._pcolorargs('pcolormesh', *args, allmatch=allmatch)
        Ny, Nx = X.shape
        X = X.ravel()
        Y = Y.ravel()
        # unit conversion allows e.g. datetime objects as axis values
        self._process_unit_info(xdata=X, ydata=Y, kwargs=kwargs)
        X = self.convert_xunits(X)
        Y = self.convert_yunits(Y)

        # convert to one dimensional arrays
        C = C.ravel()
        coords = np.column_stack((X, Y)).astype(float, copy=False)
        collection = mcoll.QuadMesh(Nx - 1, Ny - 1, coords,
                                    antialiased=antialiased, shading=shading,
                                    **kwargs)
        collection.set_alpha(alpha)
        collection.set_array(C)
        if norm is not None and not isinstance(norm, mcolors.Normalize):
            raise ValueError(
                "'norm' must be an instance of 'mcolors.Normalize'")
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_clim(vmin, vmax)
        collection.autoscale_None()

        self.grid(False)

        # Transform from native to data coordinates?
        t = collection._transform
        if (not isinstance(t, mtransforms.Transform) and
            hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)

        if t and any(t.contains_branch_seperately(self.transData)):
            trans_to_data = t - self.transData
            coords = trans_to_data.transform(coords)

        self.add_collection(collection, autolim=False)

        minx, miny = np.min(coords, axis=0)
        maxx, maxy = np.max(coords, axis=0)
        collection.sticky_edges.x[:] = [minx, maxx]
        collection.sticky_edges.y[:] = [miny, maxy]
        corners = (minx, miny), (maxx, maxy)
        self.update_datalim(corners)
        self.autoscale_view()
        return collection

    @_preprocess_data(label_namer=None)
    @docstring.dedent_interpd
    def pcolorfast(self, *args, **kwargs):
        """
        pseudocolor plot of a 2-D array

        Experimental; this is a pcolor-type method that
        provides the fastest possible rendering with the Agg
        backend, and that can handle any quadrilateral grid.
        It supports only flat shading (no outlines), it lacks
        support for log scaling of the axes, and it does not
        have a pyplot wrapper.

        Call signatures::

          ax.pcolorfast(C, **kwargs)
          ax.pcolorfast(xr, yr, C, **kwargs)
          ax.pcolorfast(x, y, C, **kwargs)
          ax.pcolorfast(X, Y, C, **kwargs)

        C is the 2D array of color values corresponding to quadrilateral
        cells. Let (nr, nc) be its shape.  C may be a masked array.

        ``ax.pcolorfast(C, **kwargs)`` is equivalent to
        ``ax.pcolorfast([0,nc], [0,nr], C, **kwargs)``

        *xr*, *yr* specify the ranges of *x* and *y* corresponding to the
        rectangular region bounding *C*.  If::

            xr = [x0, x1]

        and::

            yr = [y0,y1]

        then *x* goes from *x0* to *x1* as the second index of *C* goes
        from 0 to *nc*, etc.  (*x0*, *y0*) is the outermost corner of
        cell (0,0), and (*x1*, *y1*) is the outermost corner of cell
        (*nr*-1, *nc*-1).  All cells are rectangles of the same size.
        This is the fastest version.

        *x*, *y* are monotonic 1D arrays of length *nc* +1 and *nr* +1,
        respectively, giving the x and y boundaries of the cells.  Hence
        the cells are rectangular but the grid may be nonuniform.  The
        speed is intermediate.  (The grid is checked, and if found to be
        uniform the fast version is used.)

        *X* and *Y* are 2D arrays with shape (*nr* +1, *nc* +1) that specify
        the (x,y) coordinates of the corners of the colored
        quadrilaterals; the quadrilateral for C[i,j] has corners at
        (X[i,j],Y[i,j]), (X[i,j+1],Y[i,j+1]), (X[i+1,j],Y[i+1,j]),
        (X[i+1,j+1],Y[i+1,j+1]).  The cells need not be rectangular.
        This is the most general, but the slowest to render.  It may
        produce faster and more compact output using ps, pdf, and
        svg backends, however.

        Note that the column index corresponds to the x-coordinate,
        and the row index corresponds to y; for details, see
        :ref:`Grid Orientation <axes-pcolor-grid-orientation>`.

        Optional keyword arguments:

          *cmap*: [ *None* | Colormap ]
            A :class:`matplotlib.colors.Colormap` instance from cm. If *None*,
            use rc settings.

          *norm*: [ *None* | Normalize ]
            A :class:`matplotlib.colors.Normalize` instance is used to scale
            luminance data to 0,1. If *None*, defaults to normalize()

          *vmin*/*vmax*: [ *None* | scalar ]
            *vmin* and *vmax* are used in conjunction with norm to normalize
            luminance data.  If either are *None*, the min and max
            of the color array *C* is used.  If you pass a norm instance,
            *vmin* and *vmax* will be *None*.

          *alpha*: ``0 <= scalar <= 1``  or *None*
            the alpha blending value

        Return value is an image if a regular or rectangular grid
        is specified, and a :class:`~matplotlib.collections.QuadMesh`
        collection in the general quadrilateral case.

        """

        if not self._hold:
            self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        if norm is not None and not isinstance(norm, mcolors.Normalize):
            raise ValueError(
                "'norm' must be an instance of 'mcolors.Normalize'")

        C = args[-1]
        nr, nc = C.shape
        if len(args) == 1:
            style = "image"
            x = [0, nc]
            y = [0, nr]
        elif len(args) == 3:
            x, y = args[:2]
            x = np.asarray(x)
            y = np.asarray(y)
            if x.ndim == 1 and y.ndim == 1:
                if x.size == 2 and y.size == 2:
                    style = "image"
                else:
                    dx = np.diff(x)
                    dy = np.diff(y)
                    if (np.ptp(dx) < 0.01 * np.abs(dx.mean()) and
                        np.ptp(dy) < 0.01 * np.abs(dy.mean())):
                        style = "image"
                    else:
                        style = "pcolorimage"
            elif x.ndim == 2 and y.ndim == 2:
                style = "quadmesh"
            else:
                raise TypeError("arguments do not match valid signatures")
        else:
            raise TypeError("need 1 argument or 3 arguments")

        if style == "quadmesh":

            # convert to one dimensional arrays
            # This should also be moved to the QuadMesh class

            # data point in each cell is value at lower left corner
            C = ma.ravel(C)
            X = x.ravel()
            Y = y.ravel()
            Nx = nc + 1
            Ny = nr + 1

            # The following needs to be cleaned up; the renderer
            # requires separate contiguous arrays for X and Y,
            # but the QuadMesh class requires the 2D array.
            coords = np.empty(((Nx * Ny), 2), np.float64)
            coords[:, 0] = X
            coords[:, 1] = Y

            # The QuadMesh class can also be changed to
            # handle relevant superclass kwargs; the initializer
            # should do much more than it does now.
            collection = mcoll.QuadMesh(nc, nr, coords, 0, edgecolors="None")
            collection.set_alpha(alpha)
            collection.set_array(C)
            collection.set_cmap(cmap)
            collection.set_norm(norm)
            self.add_collection(collection, autolim=False)
            xl, xr, yb, yt = X.min(), X.max(), Y.min(), Y.max()
            ret = collection

        else:  # It's one of the two image styles.
            xl, xr, yb, yt = x[0], x[-1], y[0], y[-1]

            if style == "image":
                im = mimage.AxesImage(self, cmap, norm,
                                      interpolation='nearest',
                                      origin='lower',
                                      extent=(xl, xr, yb, yt),
                                      **kwargs)
                im.set_data(C)
                im.set_alpha(alpha)
            elif style == "pcolorimage":
                im = mimage.PcolorImage(self, x, y, C,
                                        cmap=cmap,
                                        norm=norm,
                                        alpha=alpha,
                                        **kwargs)
                im.set_extent((xl, xr, yb, yt))
            self.add_image(im)
            ret = im

        if vmin is not None or vmax is not None:
            ret.set_clim(vmin, vmax)
        else:
            ret.autoscale_None()

        ret.sticky_edges.x[:] = [xl, xr]
        ret.sticky_edges.y[:] = [yb, yt]
        self.update_datalim(np.array([[xl, yb], [xr, yt]]))
        self.autoscale_view(tight=True)
        return ret

    @_preprocess_data()
    def contour(self, *args, **kwargs):
        if not self._hold:
            self.cla()
        kwargs['filled'] = False
        contours = mcontour.QuadContourSet(self, *args, **kwargs)
        self.autoscale_view()
        return contours
    contour.__doc__ = mcontour.QuadContourSet._contour_doc

    @_preprocess_data()
    def contourf(self, *args, **kwargs):
        if not self._hold:
            self.cla()
        kwargs['filled'] = True
        contours = mcontour.QuadContourSet(self, *args, **kwargs)
        self.autoscale_view()
        return contours
    contourf.__doc__ = mcontour.QuadContourSet._contour_doc

    def clabel(self, CS, *args, **kwargs):
        return CS.clabel(*args, **kwargs)
    clabel.__doc__ = mcontour.ContourSet.clabel.__doc__

    @docstring.dedent_interpd
    def table(self, **kwargs):
        """
        Add a table to the current axes.

        Call signature::

          table(cellText=None, cellColours=None,
                cellLoc='right', colWidths=None,
                rowLabels=None, rowColours=None, rowLoc='left',
                colLabels=None, colColours=None, colLoc='center',
                loc='bottom', bbox=None)

        Returns a :class:`matplotlib.table.Table` instance. Either `cellText`
        or `cellColours` must be provided. For finer grained control over
        tables, use the :class:`~matplotlib.table.Table` class and add it to
        the axes with :meth:`~matplotlib.axes.Axes.add_table`.

        Thanks to John Gill for providing the class and table.

        kwargs control the :class:`~matplotlib.table.Table`
        properties:

        %(Table)s
        """
        return mtable.table(self, **kwargs)

    #### Data analysis

    @_preprocess_data(replace_names=["x", 'weights'], label_namer="x")
    def hist(self, x, bins=None, range=None, density=None, weights=None,
             cumulative=False, bottom=None, histtype='bar', align='mid',
             orientation='vertical', rwidth=None, log=False,
             color=None, label=None, stacked=False, normed=None,
             **kwargs):
        """
        Plot a histogram.

        Compute and draw the histogram of *x*. The return value is a
        tuple (*n*, *bins*, *patches*) or ([*n0*, *n1*, ...], *bins*,
        [*patches0*, *patches1*,...]) if the input contains multiple
        data.

        Multiple data can be provided via *x* as a list of datasets
        of potentially different length ([*x0*, *x1*, ...]), or as
        a 2-D ndarray in which each column is a dataset.  Note that
        the ndarray form is transposed relative to the list form.

        Masked arrays are not supported at present.

        Parameters
        ----------
        x : (n,) array or sequence of (n,) arrays
            Input values, this takes either a single array or a sequence of
            arrays which are not required to be of the same length

        bins : integer or sequence or 'auto', optional
            If an integer is given, ``bins + 1`` bin edges are calculated and
            returned, consistent with :func:`numpy.histogram`.

            If `bins` is a sequence, gives bin edges, including left edge of
            first bin and right edge of last bin.  In this case, `bins` is
            returned unmodified.

            All but the last (righthand-most) bin is half-open.  In other
            words, if `bins` is::

                [1, 2, 3, 4]

            then the first bin is ``[1, 2)`` (including 1, but excluding 2) and
            the second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which
            *includes* 4.

            Unequally spaced bins are supported if *bins* is a sequence.

            If Numpy 1.11 is installed, may also be ``'auto'``.

            Default is taken from the rcParam ``hist.bins``.

        range : tuple or None, optional
            The lower and upper range of the bins. Lower and upper outliers
            are ignored. If not provided, *range* is ``(x.min(), x.max())``.
            Range has no effect if *bins* is a sequence.

            If *bins* is a sequence or *range* is specified, autoscaling
            is based on the specified bin range instead of the
            range of x.

            Default is ``None``

        density : boolean, optional
            If ``True``, the first element of the return tuple will
            be the counts normalized to form a probability density, i.e.,
            the area (or integral) under the histogram will sum to 1.
            This is achieved by dividing the count by the number of
            observations times the bin width and not dividing by the total
            number of observations. If *stacked* is also ``True``, the sum of
            the histograms is normalized to 1.

            Default is ``None`` for both *normed* and *density*. If either is
            set, then that value will be used. If neither are set, then the
            args will be treated as ``False``.

            If both *density* and *normed* are set an error is raised.

        weights : (n, ) array_like or None, optional
            An array of weights, of the same shape as *x*.  Each value in *x*
            only contributes its associated weight towards the bin count
            (instead of 1).  If *normed* or *density* is ``True``,
            the weights are normalized, so that the integral of the density
            over the range remains 1.

            Default is ``None``

        cumulative : boolean, optional
            If ``True``, then a histogram is computed where each bin gives the
            counts in that bin plus all bins for smaller values. The last bin
            gives the total number of datapoints. If *normed* or *density*
            is also ``True`` then the histogram is normalized such that the
            last bin equals 1. If *cumulative* evaluates to less than 0
            (e.g., -1), the direction of accumulation is reversed.
            In this case, if *normed* and/or *density* is also ``True``, then
            the histogram is normalized such that the first bin equals 1.

            Default is ``False``

        bottom : array_like, scalar, or None
            Location of the bottom baseline of each bin.  If a scalar,
            the base line for each bin is shifted by the same amount.
            If an array, each bin is shifted independently and the length
            of bottom must match the number of bins.  If None, defaults to 0.

            Default is ``None``

        histtype : {'bar', 'barstacked', 'step',  'stepfilled'}, optional
            The type of histogram to draw.

            - 'bar' is a traditional bar-type histogram.  If multiple data
              are given the bars are arranged side by side.

            - 'barstacked' is a bar-type histogram where multiple
              data are stacked on top of each other.

            - 'step' generates a lineplot that is by default
              unfilled.

            - 'stepfilled' generates a lineplot that is by default
              filled.

            Default is 'bar'

        align : {'left', 'mid', 'right'}, optional
            Controls how the histogram is plotted.

                - 'left': bars are centered on the left bin edges.

                - 'mid': bars are centered between the bin edges.

                - 'right': bars are centered on the right bin edges.

            Default is 'mid'

        orientation : {'horizontal', 'vertical'}, optional
            If 'horizontal', `~matplotlib.pyplot.barh` will be used for
            bar-type histograms and the *bottom* kwarg will be the left edges.

        rwidth : scalar or None, optional
            The relative width of the bars as a fraction of the bin width.  If
            ``None``, automatically compute the width.

            Ignored if *histtype* is 'step' or 'stepfilled'.

            Default is ``None``

        log : boolean, optional
            If ``True``, the histogram axis will be set to a log scale. If
            *log* is ``True`` and *x* is a 1D array, empty bins will be
            filtered out and only the non-empty ``(n, bins, patches)``
            will be returned.

            Default is ``False``

        color : color or array_like of colors or None, optional
            Color spec or sequence of color specs, one per dataset.  Default
            (``None``) uses the standard line color sequence.

            Default is ``None``

        label : string or None, optional
            String, or sequence of strings to match multiple datasets.  Bar
            charts yield multiple patches per dataset, but only the first gets
            the label, so that the legend command will work as expected.

            default is ``None``

        stacked : boolean, optional
            If ``True``, multiple data are stacked on top of each other If
            ``False`` multiple data are arranged side by side if histtype is
            'bar' or on top of each other if histtype is 'step'

            Default is ``False``

        normed : bool, optional
            Deprecated; use the density keyword argument instead.

        Returns
        -------
        n : array or list of arrays
            The values of the histogram bins. See *normed* or *density*
            and *weights* for a description of the possible semantics.
            If input *x* is an array, then this is an array of length
            *nbins*. If input is a sequence arrays
            ``[data1, data2,..]``, then this is a list of arrays with
            the values of the histograms for each of the arrays in the
            same order.

        bins : array
            The edges of the bins. Length nbins + 1 (nbins left edges and right
            edge of last bin).  Always a single array even when multiple data
            sets are passed in.

        patches : list or list of lists
            Silent list of individual patches used to create the histogram
            or list of such list if multiple input datasets.

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.patches.Patch` properties

        See also
        --------
        hist2d : 2D histograms

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        # Avoid shadowing the builtin.
        bin_range = range
        del range

        if not self._hold:
            self.cla()

        if np.isscalar(x):
            x = [x]

        if bins is None:
            bins = rcParams['hist.bins']

        # Validate string inputs here so we don't have to clutter
        # subsequent code.
        if histtype not in ['bar', 'barstacked', 'step', 'stepfilled']:
            raise ValueError("histtype %s is not recognized" % histtype)

        if align not in ['left', 'mid', 'right']:
            raise ValueError("align kwarg %s is not recognized" % align)

        if orientation not in ['horizontal', 'vertical']:
            raise ValueError(
                "orientation kwarg %s is not recognized" % orientation)

        if histtype == 'barstacked' and not stacked:
            stacked = True

        if density is not None and normed is not None:
            raise ValueError("kwargs 'density' and 'normed' cannot be used "
                             "simultaneously. "
                             "Please only use 'density', since 'normed'"
                             "is deprecated.")
        if normed is not None:
            warnings.warn("The 'normed' kwarg is deprecated, and has been "
                          "replaced by the 'density' kwarg.")

        # basic input validation
        input_empty = np.size(x) == 0
        # Massage 'x' for processing.
        if input_empty:
            x = [np.array([])]
        else:
            x = cbook._reshape_2D(x, 'x')
        nx = len(x)  # number of datasets

        # Process unit information
        # Unit conversion is done individually on each dataset
        self._process_unit_info(xdata=x[0], kwargs=kwargs)
        x = [self.convert_xunits(xi) for xi in x]

        if bin_range is not None:
            bin_range = self.convert_xunits(bin_range)

        # Check whether bins or range are given explicitly.
        binsgiven = (cbook.iterable(bins) or bin_range is not None)

        # We need to do to 'weights' what was done to 'x'
        if weights is not None:
            w = cbook._reshape_2D(weights, 'weights')
        else:
            w = [None] * nx

        if len(w) != nx:
            raise ValueError('weights should have the same shape as x')

        for xi, wi in zip(x, w):
            if wi is not None and len(wi) != len(xi):
                raise ValueError(
                    'weights should have the same shape as x')

        if color is None:
            color = [self._get_lines.get_next_color() for i in xrange(nx)]
        else:
            color = mcolors.to_rgba_array(color)
            if len(color) != nx:
                raise ValueError("color kwarg must have one color per dataset")

        # If bins are not specified either explicitly or via range,
        # we need to figure out the range required for all datasets,
        # and supply that to np.histogram.
        if not binsgiven and not input_empty:
            xmin = np.inf
            xmax = -np.inf
            for xi in x:
                if len(xi) > 0:
                    xmin = min(xmin, xi.min())
                    xmax = max(xmax, xi.max())
            bin_range = (xmin, xmax)
        density = bool(density) or bool(normed)
        if density and not stacked:
            hist_kwargs = dict(range=bin_range, density=density)
        else:
            hist_kwargs = dict(range=bin_range)

        # List to store all the top coordinates of the histograms
        tops = []
        mlast = None
        # Loop through datasets
        for i in xrange(nx):
            # this will automatically overwrite bins,
            # so that each histogram uses the same bins
            m, bins = np.histogram(x[i], bins, weights=w[i], **hist_kwargs)
            m = m.astype(float)  # causes problems later if it's an int
            if mlast is None:
                mlast = np.zeros(len(bins)-1, m.dtype)
            if stacked:
                m += mlast
                mlast[:] = m
            tops.append(m)

        # If a stacked density plot, normalize so the area of all the stacked
        # histograms together is 1
        if stacked and density:
            db = np.diff(bins)
            for m in tops:
                m[:] = (m / db) / tops[-1].sum()
        if cumulative:
            slc = slice(None)
            if cbook.is_numlike(cumulative) and cumulative < 0:
                slc = slice(None, None, -1)

            if density:
                tops = [(m * np.diff(bins))[slc].cumsum()[slc] for m in tops]
            else:
                tops = [m[slc].cumsum()[slc] for m in tops]

        patches = []

        # Save autoscale state for later restoration; turn autoscaling
        # off so we can do it all a single time at the end, instead
        # of having it done by bar or fill and then having to be redone.
        _saved_autoscalex = self.get_autoscalex_on()
        _saved_autoscaley = self.get_autoscaley_on()
        self.set_autoscalex_on(False)
        self.set_autoscaley_on(False)

        if histtype.startswith('bar'):

            totwidth = np.diff(bins)

            if rwidth is not None:
                dr = np.clip(rwidth, 0, 1)
            elif (len(tops) > 1 and
                  ((not stacked) or rcParams['_internal.classic_mode'])):
                dr = 0.8
            else:
                dr = 1.0

            if histtype == 'bar' and not stacked:
                width = dr * totwidth / nx
                dw = width
                boffset = -0.5 * dr * totwidth * (1 - 1 / nx)
            elif histtype == 'barstacked' or stacked:
                width = dr * totwidth
                boffset, dw = 0.0, 0.0

            if align == 'mid' or align == 'edge':
                boffset += 0.5 * totwidth
            elif align == 'right':
                boffset += totwidth

            if orientation == 'horizontal':
                _barfunc = self.barh
                bottom_kwarg = 'left'
            else:  # orientation == 'vertical'
                _barfunc = self.bar
                bottom_kwarg = 'bottom'

            for m, c in zip(tops, color):
                if bottom is None:
                    bottom = np.zeros(len(m))
                if stacked:
                    height = m - bottom
                else:
                    height = m
                patch = _barfunc(bins[:-1]+boffset, height, width,
                                 align='center', log=log,
                                 color=c, **{bottom_kwarg: bottom})
                patches.append(patch)
                if stacked:
                    bottom[:] = m
                boffset += dw

        elif histtype.startswith('step'):
            # these define the perimeter of the polygon
            x = np.zeros(4 * len(bins) - 3)
            y = np.zeros(4 * len(bins) - 3)

            x[0:2*len(bins)-1:2], x[1:2*len(bins)-1:2] = bins, bins[:-1]
            x[2*len(bins)-1:] = x[1:2*len(bins)-1][::-1]

            if bottom is None:
                bottom = np.zeros(len(bins) - 1)

            y[1:2*len(bins)-1:2], y[2:2*len(bins):2] = bottom, bottom
            y[2*len(bins)-1:] = y[1:2*len(bins)-1][::-1]

            if log:
                if orientation == 'horizontal':
                    self.set_xscale('log', nonposx='clip')
                    logbase = self.xaxis._scale.base
                else:  # orientation == 'vertical'
                    self.set_yscale('log', nonposy='clip')
                    logbase = self.yaxis._scale.base

                # Setting a minimum of 0 results in problems for log plots
                if np.min(bottom) > 0:
                    minimum = np.min(bottom)
                elif density or weights is not None:
                    # For data that is normed to form a probability density,
                    # set to minimum data value / logbase
                    # (gives 1 full tick-label unit for the lowest filled bin)
                    ndata = np.array(tops)
                    minimum = (np.min(ndata[ndata > 0])) / logbase
                else:
                    # For non-normed (density = False) data,
                    # set the min to 1 / log base,
                    # again so that there is 1 full tick-label unit
                    # for the lowest bin
                    minimum = 1.0 / logbase

                y[0], y[-1] = minimum, minimum
            else:
                minimum = 0

            if align == 'left' or align == 'center':
                x -= 0.5*(bins[1]-bins[0])
            elif align == 'right':
                x += 0.5*(bins[1]-bins[0])

            # If fill kwarg is set, it will be passed to the patch collection,
            # overriding this
            fill = (histtype == 'stepfilled')

            xvals, yvals = [], []
            for m in tops:
                if stacked:
                    # starting point for drawing polygon
                    y[0] = y[1]
                    # top of the previous polygon becomes the bottom
                    y[2*len(bins)-1:] = y[1:2*len(bins)-1][::-1]
                # set the top of this polygon
                y[1:2*len(bins)-1:2], y[2:2*len(bins):2] = (m + bottom,
                                                            m + bottom)
                if log:
                    y[y < minimum] = minimum
                if orientation == 'horizontal':
                    xvals.append(y.copy())
                    yvals.append(x.copy())
                else:
                    xvals.append(x.copy())
                    yvals.append(y.copy())

            # stepfill is closed, step is not
            split = -1 if fill else 2 * len(bins)
            # add patches in reverse order so that when stacking,
            # items lower in the stack are plotted on top of
            # items higher in the stack
            for x, y, c in reversed(list(zip(xvals, yvals, color))):
                patches.append(self.fill(
                    x[:split], y[:split],
                    closed=True if fill else None,
                    facecolor=c,
                    edgecolor=None if fill else c,
                    fill=fill if fill else None))
            for patch_list in patches:
                for patch in patch_list:
                    if orientation == 'vertical':
                        patch.sticky_edges.y.append(minimum)
                    elif orientation == 'horizontal':
                        patch.sticky_edges.x.append(minimum)

            # we return patches, so put it back in the expected order
            patches.reverse()

        self.set_autoscalex_on(_saved_autoscalex)
        self.set_autoscaley_on(_saved_autoscaley)
        self.autoscale_view()

        if label is None:
            labels = [None]
        elif isinstance(label, six.string_types):
            labels = [label]
        else:
            labels = [six.text_type(lab) for lab in label]

        for patch, lbl in zip_longest(patches, labels, fillvalue=None):
            if patch:
                p = patch[0]
                p.update(kwargs)
                if lbl is not None:
                    p.set_label(lbl)

                for p in patch[1:]:
                    p.update(kwargs)
                    p.set_label('_nolegend_')

        if nx == 1:
            return tops[0], bins, cbook.silent_list('Patch', patches[0])
        else:
            return tops, bins, cbook.silent_list('Lists of Patches', patches)

    @_preprocess_data(replace_names=["x", "y", "weights"], label_namer=None)
    def hist2d(self, x, y, bins=10, range=None, normed=False, weights=None,
               cmin=None, cmax=None, **kwargs):
        """
        Make a 2D histogram plot.

        Parameters
        ----------
        x, y: array_like, shape (n, )
            Input values

        bins: [None | int | [int, int] | array_like | [array, array]]

            The bin specification:

                - If int, the number of bins for the two dimensions
                  (nx=ny=bins).

                - If [int, int], the number of bins in each dimension
                  (nx, ny = bins).

                - If array_like, the bin edges for the two dimensions
                  (x_edges=y_edges=bins).

                - If [array, array], the bin edges in each dimension
                  (x_edges, y_edges = bins).

            The default value is 10.

        range : array_like shape(2, 2), optional, default: None
             The leftmost and rightmost edges of the bins along each dimension
             (if not specified explicitly in the bins parameters): [[xmin,
             xmax], [ymin, ymax]]. All values outside of this range will be
             considered outliers and not tallied in the histogram.

        normed : boolean, optional, default: False
             Normalize histogram.

        weights : array_like, shape (n, ), optional, default: None
            An array of values w_i weighing each sample (x_i, y_i).

        cmin : scalar, optional, default: None
             All bins that has count less than cmin will not be displayed and
             these count values in the return value count histogram will also
             be set to nan upon return

        cmax : scalar, optional, default: None
             All bins that has count more than cmax will not be displayed (set
             to none before passing to imshow) and these count values in the
             return value count histogram will also be set to nan upon return

        Returns
        -------
        h : 2D array
            The bi-dimensional histogram of samples x and y. Values in x are
            histogrammed along the first dimension and values in y are
            histogrammed along the second dimension.
        xedges : 1D array
            The bin edges along the x axis.
        yedges : 1D array
            The bin edges along the y axis.
        image : AxesImage

        Other Parameters
        ----------------
        cmap : {Colormap, string}, optional
            A :class:`matplotlib.colors.Colormap` instance.  If not set, use rc
            settings.

        norm : Normalize, optional
            A :class:`matplotlib.colors.Normalize` instance is used to
            scale luminance data to ``[0, 1]``. If not set, defaults to
            ``Normalize()``.

        vmin/vmax : {None, scalar}, optional
            Arguments passed to the `Normalize` instance.

        alpha : ``0 <= scalar <= 1`` or ``None``, optional
            The alpha blending value.

        See also
        --------
        hist : 1D histogram

        Notes
        -----
        Rendering the histogram with a logarithmic color scale is
        accomplished by passing a :class:`colors.LogNorm` instance to
        the *norm* keyword argument. Likewise, power-law normalization
        (similar in effect to gamma correction) can be accomplished with
        :class:`colors.PowerNorm`.
        """

        h, xedges, yedges = np.histogram2d(x, y, bins=bins, range=range,
                                           normed=normed, weights=weights)

        if cmin is not None:
            h[h < cmin] = None
        if cmax is not None:
            h[h > cmax] = None

        pc = self.pcolorfast(xedges, yedges, h.T, **kwargs)
        self.set_xlim(xedges[0], xedges[-1])
        self.set_ylim(yedges[0], yedges[-1])

        return h, xedges, yedges, pc

    @_preprocess_data(replace_names=["x"], label_namer=None)
    @docstring.dedent_interpd
    def psd(self, x, NFFT=None, Fs=None, Fc=None, detrend=None,
            window=None, noverlap=None, pad_to=None,
            sides=None, scale_by_freq=None, return_line=None, **kwargs):
        r"""
        Plot the power spectral density.

        Call signature::

          psd(x, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none,
              window=mlab.window_hanning, noverlap=0, pad_to=None,
              sides='default', scale_by_freq=None, return_line=None, **kwargs)

        The power spectral density :math:`P_{xx}` by Welch's average
        periodogram method.  The vector *x* is divided into *NFFT* length
        segments.  Each segment is detrended by function *detrend* and
        windowed by function *window*.  *noverlap* gives the length of
        the overlap between segments.  The :math:`|\mathrm{fft}(i)|^2`
        of each segment :math:`i` are averaged to compute :math:`P_{xx}`,
        with a scaling to correct for power loss due to windowing.

        If len(*x*) < *NFFT*, it will be zero padded to *NFFT*.

        Parameters
        ----------
        x : 1-D array or sequence
            Array or sequence containing the data

        %(Spectral)s

        %(PSD)s

        noverlap : integer
            The number of points of overlap between segments.
            The default value is 0 (no overlap).

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        return_line : bool
            Whether to include the line object plotted in the returned values.
            Default is False.

        Returns
        -------
        Pxx : 1-D array
            The values for the power spectrum `P_{xx}` before scaling
            (real valued)

        freqs : 1-D array
            The frequencies corresponding to the elements in *Pxx*

        line : a :class:`~matplotlib.lines.Line2D` instance
            The line created by this function.
            Only returned if *return_line* is True.

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        See Also
        --------
        :func:`specgram`
            :func:`specgram` differs in the default overlap; in not returning
            the mean of the segment periodograms; in returning the times of the
            segments; and in plotting a colormap instead of a line.

        :func:`magnitude_spectrum`
            :func:`magnitude_spectrum` plots the magnitude spectrum.

        :func:`csd`
            :func:`csd` plots the spectral density between two signals.

        Notes
        -----
        For plotting, the power is plotted as
        :math:`10\log_{10}(P_{xx})` for decibels, though *Pxx* itself
        is returned.

        References
        ----------
        Bendat & Piersol -- Random Data: Analysis and Measurement Procedures,
        John Wiley & Sons (1986)
        """
        if not self._hold:
            self.cla()

        if Fc is None:
            Fc = 0

        pxx, freqs = mlab.psd(x=x, NFFT=NFFT, Fs=Fs, detrend=detrend,
                              window=window, noverlap=noverlap, pad_to=pad_to,
                              sides=sides, scale_by_freq=scale_by_freq)
        freqs += Fc

        if scale_by_freq in (None, True):
            psd_units = 'dB/Hz'
        else:
            psd_units = 'dB'

        line = self.plot(freqs, 10 * np.log10(pxx), **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Power Spectral Density (%s)' % psd_units)
        self.grid(True)
        vmin, vmax = self.viewLim.intervaly
        intv = vmax - vmin
        logi = int(np.log10(intv))
        if logi == 0:
            logi = .1
        step = 10 * logi
        ticks = np.arange(math.floor(vmin), math.ceil(vmax) + 1, step)
        self.set_yticks(ticks)

        if return_line is None or not return_line:
            return pxx, freqs
        else:
            return pxx, freqs, line

    @_preprocess_data(replace_names=["x", "y"], label_namer="y")
    @docstring.dedent_interpd
    def csd(self, x, y, NFFT=None, Fs=None, Fc=None, detrend=None,
            window=None, noverlap=None, pad_to=None,
            sides=None, scale_by_freq=None, return_line=None, **kwargs):
        """
        Plot the cross-spectral density.

        Call signature::

          csd(x, y, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none,
              window=mlab.window_hanning, noverlap=0, pad_to=None,
              sides='default', scale_by_freq=None, return_line=None, **kwargs)

        The cross spectral density :math:`P_{xy}` by Welch's average
        periodogram method.  The vectors *x* and *y* are divided into
        *NFFT* length segments.  Each segment is detrended by function
        *detrend* and windowed by function *window*.  *noverlap* gives
        the length of the overlap between segments.  The product of
        the direct FFTs of *x* and *y* are averaged over each segment
        to compute :math:`P_{xy}`, with a scaling to correct for power
        loss due to windowing.

        If len(*x*) < *NFFT* or len(*y*) < *NFFT*, they will be zero
        padded to *NFFT*.

        Parameters
        ----------
        x, y : 1-D arrays or sequences
            Arrays or sequences containing the data

        %(Spectral)s

        %(PSD)s

        noverlap : integer
            The number of points of overlap between segments.
            The default value is 0 (no overlap).

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        return_line : bool
            Whether to include the line object plotted in the returned values.
            Default is False.

        Returns
        -------
        Pxy : 1-D array
            The values for the cross spectrum `P_{xy}` before scaling
            (complex valued)

        freqs : 1-D array
            The frequencies corresponding to the elements in *Pxy*

        line : a :class:`~matplotlib.lines.Line2D` instance
            The line created by this function.
            Only returned if *return_line* is True.

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        See Also
        --------
        :func:`psd`
            :func:`psd` is the equivalent to setting y=x.

        Notes
        -----
        For plotting, the power is plotted as
        :math:`10\\log_{10}(P_{xy})` for decibels, though `P_{xy}` itself
        is returned.

        References
        ----------
        Bendat & Piersol -- Random Data: Analysis and Measurement Procedures,
        John Wiley & Sons (1986)
        """
        if not self._hold:
            self.cla()

        if Fc is None:
            Fc = 0

        pxy, freqs = mlab.csd(x=x, y=y, NFFT=NFFT, Fs=Fs, detrend=detrend,
                              window=window, noverlap=noverlap, pad_to=pad_to,
                              sides=sides, scale_by_freq=scale_by_freq)
        # pxy is complex
        freqs += Fc

        line = self.plot(freqs, 10 * np.log10(np.abs(pxy)), **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Cross Spectrum Magnitude (dB)')
        self.grid(True)
        vmin, vmax = self.viewLim.intervaly

        intv = vmax - vmin
        step = 10 * int(np.log10(intv))

        ticks = np.arange(math.floor(vmin), math.ceil(vmax) + 1, step)
        self.set_yticks(ticks)

        if return_line is None or not return_line:
            return pxy, freqs
        else:
            return pxy, freqs, line

    @_preprocess_data(replace_names=["x"], label_namer=None)
    @docstring.dedent_interpd
    def magnitude_spectrum(self, x, Fs=None, Fc=None, window=None,
                           pad_to=None, sides=None, scale=None,
                           **kwargs):
        """
        Plot the magnitude spectrum.

        Call signature::

          magnitude_spectrum(x, Fs=2, Fc=0,  window=mlab.window_hanning,
                             pad_to=None, sides='default', **kwargs)

        Compute the magnitude spectrum of *x*.  Data is padded to a
        length of *pad_to* and the windowing function *window* is applied to
        the signal.

        Parameters
        ----------
        x : 1-D array or sequence
            Array or sequence containing the data

        %(Spectral)s

        %(Single_Spectrum)s

        scale : [ 'default' | 'linear' | 'dB' ]
            The scaling of the values in the *spec*.  'linear' is no scaling.
            'dB' returns the values in dB scale, i.e., the dB amplitude
            (20 * log10). 'default' is 'linear'.

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        Returns
        -------
        spectrum : 1-D array
            The values for the magnitude spectrum before scaling (real valued)

        freqs : 1-D array
            The frequencies corresponding to the elements in *spectrum*

        line : a :class:`~matplotlib.lines.Line2D` instance
            The line created by this function

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        See Also
        --------
        :func:`psd`
            :func:`psd` plots the power spectral density.`.

        :func:`angle_spectrum`
            :func:`angle_spectrum` plots the angles of the corresponding
            frequencies.

        :func:`phase_spectrum`
            :func:`phase_spectrum` plots the phase (unwrapped angle) of the
            corresponding frequencies.

        :func:`specgram`
            :func:`specgram` can plot the magnitude spectrum of segments within
            the signal in a colormap.

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        if not self._hold:
            self.cla()

        if Fc is None:
            Fc = 0

        if scale is None or scale == 'default':
            scale = 'linear'

        spec, freqs = mlab.magnitude_spectrum(x=x, Fs=Fs, window=window,
                                              pad_to=pad_to, sides=sides)
        freqs += Fc

        if scale == 'linear':
            Z = spec
            yunits = 'energy'
        elif scale == 'dB':
            Z = 20. * np.log10(spec)
            yunits = 'dB'
        else:
            raise ValueError('Unknown scale %s', scale)

        lines = self.plot(freqs, Z, **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Magnitude (%s)' % yunits)

        return spec, freqs, lines[0]

    @_preprocess_data(replace_names=["x"], label_namer=None)
    @docstring.dedent_interpd
    def angle_spectrum(self, x, Fs=None, Fc=None, window=None,
                       pad_to=None, sides=None, **kwargs):
        """
        Plot the angle spectrum.

        Call signature::

          angle_spectrum(x, Fs=2, Fc=0,  window=mlab.window_hanning,
                         pad_to=None, sides='default', **kwargs)

        Compute the angle spectrum (wrapped phase spectrum) of *x*.
        Data is padded to a length of *pad_to* and the windowing function
        *window* is applied to the signal.

        Parameters
        ----------
        x : 1-D array or sequence
            Array or sequence containing the data

        %(Spectral)s

        %(Single_Spectrum)s

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        Returns
        -------
        spectrum : 1-D array
            The values for the angle spectrum in radians (real valued)

        freqs : 1-D array
            The frequencies corresponding to the elements in *spectrum*

        line : a :class:`~matplotlib.lines.Line2D` instance
            The line created by this function

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        See Also
        --------
        :func:`magnitude_spectrum`
            :func:`angle_spectrum` plots the magnitudes of the corresponding
            frequencies.

        :func:`phase_spectrum`
            :func:`phase_spectrum` plots the unwrapped version of this
            function.

        :func:`specgram`
            :func:`specgram` can plot the angle spectrum of segments within the
            signal in a colormap.

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        if not self._hold:
            self.cla()

        if Fc is None:
            Fc = 0

        spec, freqs = mlab.angle_spectrum(x=x, Fs=Fs, window=window,
                                          pad_to=pad_to, sides=sides)
        freqs += Fc

        lines = self.plot(freqs, spec, **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Angle (radians)')

        return spec, freqs, lines[0]

    @_preprocess_data(replace_names=["x"], label_namer=None)
    @docstring.dedent_interpd
    def phase_spectrum(self, x, Fs=None, Fc=None, window=None,
                       pad_to=None, sides=None, **kwargs):
        """
        Plot the phase spectrum.

        Call signature::

          phase_spectrum(x, Fs=2, Fc=0,  window=mlab.window_hanning,
                         pad_to=None, sides='default', **kwargs)

        Compute the phase spectrum (unwrapped angle spectrum) of *x*.
        Data is padded to a length of *pad_to* and the windowing function
        *window* is applied to the signal.

        Parameters
        ----------
        x : 1-D array or sequence
            Array or sequence containing the data

        %(Spectral)s

        %(Single_Spectrum)s

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        Returns
        -------
        spectrum : 1-D array
            The values for the phase spectrum in radians (real valued)

        freqs : 1-D array
            The frequencies corresponding to the elements in *spectrum*

        line : a :class:`~matplotlib.lines.Line2D` instance
            The line created by this function

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        See Also
        --------
        :func:`magnitude_spectrum`
            :func:`magnitude_spectrum` plots the magnitudes of the
            corresponding frequencies.

        :func:`angle_spectrum`
            :func:`angle_spectrum` plots the wrapped version of this function.

        :func:`specgram`
            :func:`specgram` can plot the phase spectrum of segments within the
            signal in a colormap.

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """
        if not self._hold:
            self.cla()

        if Fc is None:
            Fc = 0

        spec, freqs = mlab.phase_spectrum(x=x, Fs=Fs, window=window,
                                          pad_to=pad_to, sides=sides)
        freqs += Fc

        lines = self.plot(freqs, spec, **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Phase (radians)')

        return spec, freqs, lines[0]

    @_preprocess_data(replace_names=["x", "y"], label_namer=None)
    @docstring.dedent_interpd
    def cohere(self, x, y, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none,
               window=mlab.window_hanning, noverlap=0, pad_to=None,
               sides='default', scale_by_freq=None, **kwargs):
        """
        Plot the coherence between *x* and *y*.

        Plot the coherence between *x* and *y*.  Coherence is the
        normalized cross spectral density:

        .. math::

          C_{xy} = \\frac{|P_{xy}|^2}{P_{xx}P_{yy}}

        Parameters
        ----------
        %(Spectral)s

        %(PSD)s

        noverlap : integer
            The number of points of overlap between blocks.  The
            default value is 0 (no overlap).

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.


        Returns
        -------
        The return value is a tuple (*Cxy*, *f*), where *f* are the
        frequencies of the coherence vector.

        kwargs are applied to the lines.

        Other Parameters
        ----------------
        **kwargs :
            Keyword arguments control the :class:`~matplotlib.lines.Line2D`
            properties:

            %(Line2D)s

        References
        ----------
        Bendat & Piersol -- Random Data: Analysis and Measurement Procedures,
        John Wiley & Sons (1986)
        """
        if not self._hold:
            self.cla()
        cxy, freqs = mlab.cohere(x=x, y=y, NFFT=NFFT, Fs=Fs, detrend=detrend,
                                 window=window, noverlap=noverlap,
                                 scale_by_freq=scale_by_freq)
        freqs += Fc

        self.plot(freqs, cxy, **kwargs)
        self.set_xlabel('Frequency')
        self.set_ylabel('Coherence')
        self.grid(True)

        return cxy, freqs

    @_preprocess_data(replace_names=["x"], label_namer=None)
    @docstring.dedent_interpd
    def specgram(self, x, NFFT=None, Fs=None, Fc=None, detrend=None,
                 window=None, noverlap=None,
                 cmap=None, xextent=None, pad_to=None, sides=None,
                 scale_by_freq=None, mode=None, scale=None,
                 vmin=None, vmax=None, **kwargs):
        """
        Plot a spectrogram.

        Call signature::

          specgram(x, NFFT=256, Fs=2, Fc=0, detrend=mlab.detrend_none,
                   window=mlab.window_hanning, noverlap=128,
                   cmap=None, xextent=None, pad_to=None, sides='default',
                   scale_by_freq=None, mode='default', scale='default',
                   **kwargs)

        Compute and plot a spectrogram of data in *x*.  Data are split into
        *NFFT* length segments and the spectrum of each section is
        computed.  The windowing function *window* is applied to each
        segment, and the amount of overlap of each segment is
        specified with *noverlap*. The spectrogram is plotted as a colormap
        (using imshow).

        Parameters
        ----------
        x : 1-D array or sequence
            Array or sequence containing the data.

        %(Spectral)s

        %(PSD)s

        mode : [ 'default' | 'psd' | 'magnitude' | 'angle' | 'phase' ]
            What sort of spectrum to use.  Default is 'psd', which takes
            the power spectral density.  'complex' returns the complex-valued
            frequency spectrum.  'magnitude' returns the magnitude spectrum.
            'angle' returns the phase spectrum without unwrapping.  'phase'
            returns the phase spectrum with unwrapping.

        noverlap : integer
            The number of points of overlap between blocks.  The
            default value is 128.

        scale : [ 'default' | 'linear' | 'dB' ]
            The scaling of the values in the *spec*.  'linear' is no scaling.
            'dB' returns the values in dB scale.  When *mode* is 'psd',
            this is dB power (10 * log10).  Otherwise this is dB amplitude
            (20 * log10). 'default' is 'dB' if *mode* is 'psd' or
            'magnitude' and 'linear' otherwise.  This must be 'linear'
            if *mode* is 'angle' or 'phase'.

        Fc : integer
            The center frequency of *x* (defaults to 0), which offsets
            the x extents of the plot to reflect the frequency range used
            when a signal is acquired and then filtered and downsampled to
            baseband.

        cmap :
            A :class:`matplotlib.colors.Colormap` instance; if *None*, use
            default determined by rc

        xextent : [None | (xmin, xmax)]
            The image extent along the x-axis. The default sets *xmin* to the
            left border of the first bin (*spectrum* column) and *xmax* to the
            right border of the last bin. Note that for *noverlap>0* the width
            of the bins is smaller than those of the segments.

        **kwargs :
            Additional kwargs are passed on to imshow which makes the
            specgram image

        Returns
        -------
        spectrum : 2-D array
            Columns are the periodograms of successive segments.

        freqs : 1-D array
            The frequencies corresponding to the rows in *spectrum*.

        t : 1-D array
            The times corresponding to midpoints of segments (i.e., the columns
            in *spectrum*).

        im : instance of class :class:`~matplotlib.image.AxesImage`
            The image created by imshow containing the spectrogram

        See Also
        --------
        :func:`psd`
            :func:`psd` differs in the default overlap; in returning the mean
            of the segment periodograms; in not returning times; and in
            generating a line plot instead of colormap.

        :func:`magnitude_spectrum`
            A single spectrum, similar to having a single segment when *mode*
            is 'magnitude'. Plots a line instead of a colormap.

        :func:`angle_spectrum`
            A single spectrum, similar to having a single segment when *mode*
            is 'angle'. Plots a line instead of a colormap.

        :func:`phase_spectrum`
            A single spectrum, similar to having a single segment when *mode*
            is 'phase'. Plots a line instead of a colormap.

        Notes
        -----
        The parameters *detrend* and *scale_by_freq* do only apply when *mode*
        is set to 'psd'.
        """
        if not self._hold:
            self.cla()

        if NFFT is None:
            NFFT = 256  # same default as in mlab.specgram()
        if Fc is None:
            Fc = 0  # same default as in mlab._spectral_helper()
        if noverlap is None:
            noverlap = 128  # same default as in mlab.specgram()

        if mode == 'complex':
            raise ValueError('Cannot plot a complex specgram')

        if scale is None or scale == 'default':
            if mode in ['angle', 'phase']:
                scale = 'linear'
            else:
                scale = 'dB'
        elif mode in ['angle', 'phase'] and scale == 'dB':
            raise ValueError('Cannot use dB scale with angle or phase mode')

        spec, freqs, t = mlab.specgram(x=x, NFFT=NFFT, Fs=Fs,
                                       detrend=detrend, window=window,
                                       noverlap=noverlap, pad_to=pad_to,
                                       sides=sides,
                                       scale_by_freq=scale_by_freq,
                                       mode=mode)

        if scale == 'linear':
            Z = spec
        elif scale == 'dB':
            if mode is None or mode == 'default' or mode == 'psd':
                Z = 10. * np.log10(spec)
            else:
                Z = 20. * np.log10(spec)
        else:
            raise ValueError('Unknown scale %s', scale)

        Z = np.flipud(Z)

        if xextent is None:
            # padding is needed for first and last segment:
            pad_xextent = (NFFT-noverlap) / Fs / 2
            xextent = np.min(t) - pad_xextent, np.max(t) + pad_xextent
        xmin, xmax = xextent
        freqs += Fc
        extent = xmin, xmax, freqs[0], freqs[-1]
        im = self.imshow(Z, cmap, extent=extent, vmin=vmin, vmax=vmax,
                         **kwargs)
        self.axis('auto')

        return spec, freqs, t, im

    def spy(self, Z, precision=0, marker=None, markersize=None,
            aspect='equal', origin="upper", **kwargs):
        """
        Plot the sparsity pattern on a 2-D array.

        ``spy(Z)`` plots the sparsity pattern of the 2-D array *Z*.

        Parameters
        ----------

        Z : sparse array (n, m)
            The array to be plotted.

        precision : float, optional, default: 0
            If *precision* is 0, any non-zero value will be plotted; else,
            values of :math:`|Z| > precision` will be plotted.

            For :class:`scipy.sparse.spmatrix` instances, there is a special
            case: if *precision* is 'present', any value present in the array
            will be plotted, even if it is identically zero.

        origin : ["upper", "lower"], optional, default: "upper"
            Place the [0,0] index of the array in the upper left or lower left
            corner of the axes.

        aspect : ['auto' | 'equal' | scalar], optional, default: "equal"

            If 'equal', and `extent` is None, changes the axes aspect ratio to
            match that of the image. If `extent` is not `None`, the axes
            aspect ratio is changed to match that of the extent.


            If 'auto', changes the image aspect ratio to match that of the
            axes.

            If None, default to rc ``image.aspect`` value.

        Two plotting styles are available: image or marker. Both
        are available for full arrays, but only the marker style
        works for :class:`scipy.sparse.spmatrix` instances.

        If *marker* and *markersize* are *None*, an image will be
        returned and any remaining kwargs are passed to
        :func:`~matplotlib.pyplot.imshow`; else, a
        :class:`~matplotlib.lines.Line2D` object will be returned with
        the value of marker determining the marker type, and any
        remaining kwargs passed to the
        :meth:`~matplotlib.axes.Axes.plot` method.

        If *marker* and *markersize* are *None*, useful kwargs include:

        * *cmap*
        * *alpha*

        See also
        --------
        imshow : for image options.
        plot : for plotting options
        """
        if marker is None and markersize is None and hasattr(Z, 'tocoo'):
            marker = 's'
        if marker is None and markersize is None:
            Z = np.asarray(Z)
            mask = np.abs(Z) > precision

            if 'cmap' not in kwargs:
                kwargs['cmap'] = mcolors.ListedColormap(['w', 'k'],
                                                        name='binary')
            nr, nc = Z.shape
            extent = [-0.5, nc - 0.5, nr - 0.5, -0.5]
            ret = self.imshow(mask, interpolation='nearest', aspect=aspect,
                                extent=extent, origin=origin, **kwargs)
        else:
            if hasattr(Z, 'tocoo'):
                c = Z.tocoo()
                if precision == 'present':
                    y = c.row
                    x = c.col
                else:
                    nonzero = np.abs(c.data) > precision
                    y = c.row[nonzero]
                    x = c.col[nonzero]
            else:
                Z = np.asarray(Z)
                nonzero = np.abs(Z) > precision
                y, x = np.nonzero(nonzero)
            if marker is None:
                marker = 's'
            if markersize is None:
                markersize = 10
            marks = mlines.Line2D(x, y, linestyle='None',
                         marker=marker, markersize=markersize, **kwargs)
            self.add_line(marks)
            nr, nc = Z.shape
            self.set_xlim(xmin=-0.5, xmax=nc - 0.5)
            self.set_ylim(ymin=nr - 0.5, ymax=-0.5)
            self.set_aspect(aspect)
            ret = marks
        self.title.set_y(1.05)
        self.xaxis.tick_top()
        self.xaxis.set_ticks_position('both')
        self.xaxis.set_major_locator(mticker.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
        self.yaxis.set_major_locator(mticker.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
        return ret

    def matshow(self, Z, **kwargs):
        """
        Plot a matrix or array as an image.

        The matrix will be shown the way it would be printed, with the first
        row at the top.  Row and column numbering is zero-based.

        Parameters
        ----------
        Z : array_like shape (n, m)
            The matrix to be displayed.

        Returns
        -------
        image : `~matplotlib.image.AxesImage`

        Other Parameters
        ----------------
        **kwargs : `~matplotlib.axes.Axes.imshow` arguments
            Sets `origin` to 'upper', 'interpolation' to 'nearest' and
            'aspect' to equal.

        See also
        --------
        imshow : plot an image

        """
        Z = np.asanyarray(Z)
        nr, nc = Z.shape
        kw = {'origin': 'upper',
              'interpolation': 'nearest',
              'aspect': 'equal'}          # (already the imshow default)
        kw.update(kwargs)
        im = self.imshow(Z, **kw)
        self.title.set_y(1.05)
        self.xaxis.tick_top()
        self.xaxis.set_ticks_position('both')
        self.xaxis.set_major_locator(mticker.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
        self.yaxis.set_major_locator(mticker.MaxNLocator(nbins=9,
                                                 steps=[1, 2, 5, 10],
                                                 integer=True))
        return im

    @_preprocess_data(replace_names=["dataset"], label_namer=None)
    def violinplot(self, dataset, positions=None, vert=True, widths=0.5,
                   showmeans=False, showextrema=True, showmedians=False,
                   points=100, bw_method=None):
        """
        Make a violin plot.

        Make a violin plot for each column of *dataset* or each vector in
        sequence *dataset*.  Each filled area extends to represent the
        entire data range, with optional lines at the mean, the median,
        the minimum, and the maximum.

        Parameters
        ----------
        dataset : Array or a sequence of vectors.
          The input data.

        positions : array-like, default = [1, 2, ..., n]
          Sets the positions of the violins. The ticks and limits are
          automatically set to match the positions.

        vert : bool, default = True.
          If true, creates a vertical violin plot.
          Otherwise, creates a horizontal violin plot.

        widths : array-like, default = 0.5
          Either a scalar or a vector that sets the maximal width of
          each violin. The default is 0.5, which uses about half of the
          available horizontal space.

        showmeans : bool, default = False
          If `True`, will toggle rendering of the means.

        showextrema : bool, default = True
          If `True`, will toggle rendering of the extrema.

        showmedians : bool, default = False
          If `True`, will toggle rendering of the medians.

        points : scalar, default = 100
          Defines the number of points to evaluate each of the
          gaussian kernel density estimations at.

        bw_method : str, scalar or callable, optional
          The method used to calculate the estimator bandwidth.  This can be
          'scott', 'silverman', a scalar constant or a callable.  If a
          scalar, this will be used directly as `kde.factor`.  If a
          callable, it should take a `GaussianKDE` instance as its only
          parameter and return a scalar. If None (default), 'scott' is used.

        Returns
        -------

        result : dict
          A dictionary mapping each component of the violinplot to a
          list of the corresponding collection instances created. The
          dictionary has the following keys:

            - ``bodies``: A list of the
              :class:`matplotlib.collections.PolyCollection` instances
              containing the filled area of each violin.

            - ``cmeans``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the mean values of each of the
              violin's distribution.

            - ``cmins``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the bottom of each violin's
              distribution.

            - ``cmaxes``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the top of each violin's
              distribution.

            - ``cbars``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the centers of each violin's
              distribution.

            - ``cmedians``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the median values of each of the
              violin's distribution.

        Notes
        -----
        .. [Notes section required for data comment. See #10189.]

        """

        def _kde_method(X, coords):
            # fallback gracefully if the vector contains only one value
            if np.all(X[0] == X):
                return (X[0] == coords).astype(float)
            kde = mlab.GaussianKDE(X, bw_method)
            return kde.evaluate(coords)

        vpstats = cbook.violin_stats(dataset, _kde_method, points=points)
        return self.violin(vpstats, positions=positions, vert=vert,
                           widths=widths, showmeans=showmeans,
                           showextrema=showextrema, showmedians=showmedians)

    def violin(self, vpstats, positions=None, vert=True, widths=0.5,
               showmeans=False, showextrema=True, showmedians=False):
        """Drawing function for violin plots.

        Draw a violin plot for each column of `vpstats`. Each filled area
        extends to represent the entire data range, with optional lines at the
        mean, the median, the minimum, and the maximum.

        Parameters
        ----------

        vpstats : list of dicts
          A list of dictionaries containing stats for each violin plot.
          Required keys are:

          - ``coords``: A list of scalars containing the coordinates that
            the violin's kernel density estimate were evaluated at.

          - ``vals``: A list of scalars containing the values of the
            kernel density estimate at each of the coordinates given
            in *coords*.

          - ``mean``: The mean value for this violin's dataset.

          - ``median``: The median value for this violin's dataset.

          - ``min``: The minimum value for this violin's dataset.

          - ``max``: The maximum value for this violin's dataset.

        positions : array-like, default = [1, 2, ..., n]
          Sets the positions of the violins. The ticks and limits are
          automatically set to match the positions.

        vert : bool, default = True.
          If true, plots the violins veritcally.
          Otherwise, plots the violins horizontally.

        widths : array-like, default = 0.5
          Either a scalar or a vector that sets the maximal width of
          each violin. The default is 0.5, which uses about half of the
          available horizontal space.

        showmeans : bool, default = False
          If true, will toggle rendering of the means.

        showextrema : bool, default = True
          If true, will toggle rendering of the extrema.

        showmedians : bool, default = False
          If true, will toggle rendering of the medians.

        Returns
        -------
        result : dict
          A dictionary mapping each component of the violinplot to a
          list of the corresponding collection instances created. The
          dictionary has the following keys:

            - ``bodies``: A list of the
              :class:`matplotlib.collections.PolyCollection` instances
              containing the filled area of each violin.

            - ``cmeans``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the mean values of each of the
              violin's distribution.

            - ``cmins``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the bottom of each violin's
              distribution.

            - ``cmaxes``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the top of each violin's
              distribution.

            - ``cbars``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the centers of each violin's
              distribution.

            - ``cmedians``: A
              :class:`matplotlib.collections.LineCollection` instance
              created to identify the median values of each of the
              violin's distribution.

        """

        # Statistical quantities to be plotted on the violins
        means = []
        mins = []
        maxes = []
        medians = []

        # Collections to be returned
        artists = {}

        N = len(vpstats)
        datashape_message = ("List of violinplot statistics and `{0}` "
                             "values must have the same length")

        # Validate positions
        if positions is None:
            positions = range(1, N + 1)
        elif len(positions) != N:
            raise ValueError(datashape_message.format("positions"))

        # Validate widths
        if np.isscalar(widths):
            widths = [widths] * N
        elif len(widths) != N:
            raise ValueError(datashape_message.format("widths"))

        # Calculate ranges for statistics lines
        pmins = -0.25 * np.array(widths) + positions
        pmaxes = 0.25 * np.array(widths) + positions

        # Check whether we are rendering vertically or horizontally
        if vert:
            fill = self.fill_betweenx
            perp_lines = self.hlines
            par_lines = self.vlines
        else:
            fill = self.fill_between
            perp_lines = self.vlines
            par_lines = self.hlines

        if rcParams['_internal.classic_mode']:
            fillcolor = 'y'
            edgecolor = 'r'
        else:
            fillcolor = edgecolor = self._get_lines.get_next_color()

        # Render violins
        bodies = []
        for stats, pos, width in zip(vpstats, positions, widths):
            # The 0.5 factor reflects the fact that we plot from v-p to
            # v+p
            vals = np.array(stats['vals'])
            vals = 0.5 * width * vals / vals.max()
            bodies += [fill(stats['coords'],
                            -vals + pos,
                            vals + pos,
                            facecolor=fillcolor,
                            alpha=0.3)]
            means.append(stats['mean'])
            mins.append(stats['min'])
            maxes.append(stats['max'])
            medians.append(stats['median'])
        artists['bodies'] = bodies

        # Render means
        if showmeans:
            artists['cmeans'] = perp_lines(means, pmins, pmaxes,
                                           colors=edgecolor)

        # Render extrema
        if showextrema:
            artists['cmaxes'] = perp_lines(maxes, pmins, pmaxes,
                                           colors=edgecolor)
            artists['cmins'] = perp_lines(mins, pmins, pmaxes,
                                          colors=edgecolor)
            artists['cbars'] = par_lines(positions, mins, maxes,
                                         colors=edgecolor)

        # Render medians
        if showmedians:
            artists['cmedians'] = perp_lines(medians,
                                             pmins,
                                             pmaxes,
                                             colors=edgecolor)

        return artists

    def tricontour(self, *args, **kwargs):
        return mtri.tricontour(self, *args, **kwargs)
    tricontour.__doc__ = mtri.tricontour.__doc__

    def tricontourf(self, *args, **kwargs):
        return mtri.tricontourf(self, *args, **kwargs)
    tricontourf.__doc__ = mtri.tricontour.__doc__

    def tripcolor(self, *args, **kwargs):
        return mtri.tripcolor(self, *args, **kwargs)
    tripcolor.__doc__ = mtri.tripcolor.__doc__

    def triplot(self, *args, **kwargs):
        return mtri.triplot(self, *args, **kwargs)
    triplot.__doc__ = mtri.triplot.__doc__
