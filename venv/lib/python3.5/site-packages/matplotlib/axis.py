"""
Classes for the ticks and x and y axis
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import logging

from matplotlib import rcParams
import matplotlib.artist as artist
from matplotlib.artist import allow_rasterization
import matplotlib.cbook as cbook
from matplotlib.cbook import _string_to_bool
import matplotlib.font_manager as font_manager
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.scale as mscale
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import matplotlib.units as munits
import numpy as np
import warnings

_log = logging.getLogger(__name__)

GRIDLINE_INTERPOLATION_STEPS = 180

# This list is being used for compatibility with Axes.grid, which
# allows all Line2D kwargs.
_line_AI = artist.ArtistInspector(mlines.Line2D)
_line_param_names = _line_AI.get_setters()
_line_param_aliases = [list(d.keys())[0] for d in _line_AI.aliasd.values()]
_gridline_param_names = ['grid_' + name
                         for name in _line_param_names + _line_param_aliases]


class Tick(artist.Artist):
    """
    Abstract base class for the axis ticks, grid lines and labels

    1 refers to the bottom of the plot for xticks and the left for yticks
    2 refers to the top of the plot for xticks and the right for yticks

    Attributes
    ----------
    tick1line : Line2D

    tick2line : Line2D

    gridline : Line2D

    label1 : Text

    label2 : Text

    gridOn : bool
        Determines whether to draw the tickline.

    tick1On : bool
        Determines whether to draw the first tickline.

    tick2On : bool
        Determines whether to draw the second tickline.

    label1On : bool
        Determines whether to draw the first tick label.

    label2On : bool
        Determines whether to draw the second tick label.
    """
    def __init__(self, axes, loc, label,
                 size=None,  # points
                 width=None,
                 color=None,
                 tickdir=None,
                 pad=None,
                 labelsize=None,
                 labelcolor=None,
                 zorder=None,
                 gridOn=None,  # defaults to axes.grid depending on
                               # axes.grid.which
                 tick1On=True,
                 tick2On=True,
                 label1On=True,
                 label2On=False,
                 major=True,
                 labelrotation=0,
                 grid_color=None,
                 grid_linestyle=None,
                 grid_linewidth=None,
                 grid_alpha=None,
                 **kw  # Other Line2D kwargs applied to gridlines.
                 ):
        """
        bbox is the Bound2D bounding box in display coords of the Axes
        loc is the tick location in data coords
        size is the tick size in points
        """
        artist.Artist.__init__(self)

        if gridOn is None:
            if major and (rcParams['axes.grid.which'] in ('both', 'major')):
                gridOn = rcParams['axes.grid']
            elif (not major) and (rcParams['axes.grid.which']
                                  in ('both', 'minor')):
                gridOn = rcParams['axes.grid']
            else:
                gridOn = False

        self.set_figure(axes.figure)
        self.axes = axes

        name = self.__name__.lower()
        self._name = name

        self._loc = loc

        if size is None:
            if major:
                size = rcParams['%s.major.size' % name]
            else:
                size = rcParams['%s.minor.size' % name]
        self._size = size

        if width is None:
            if major:
                width = rcParams['%s.major.width' % name]
            else:
                width = rcParams['%s.minor.width' % name]
        self._width = width

        if color is None:
            color = rcParams['%s.color' % name]
        self._color = color

        if pad is None:
            if major:
                pad = rcParams['%s.major.pad' % name]
            else:
                pad = rcParams['%s.minor.pad' % name]
        self._base_pad = pad

        if labelcolor is None:
            labelcolor = rcParams['%s.color' % name]
        self._labelcolor = labelcolor

        if labelsize is None:
            labelsize = rcParams['%s.labelsize' % name]
        self._labelsize = labelsize

        self._set_labelrotation(labelrotation)

        if zorder is None:
            if major:
                zorder = mlines.Line2D.zorder + 0.01
            else:
                zorder = mlines.Line2D.zorder
        self._zorder = zorder

        self._grid_color = (rcParams['grid.color']
                            if grid_color is None else grid_color)
        self._grid_linestyle = (rcParams['grid.linestyle']
                                if grid_linestyle is None else grid_linestyle)
        self._grid_linewidth = (rcParams['grid.linewidth']
                                if grid_linewidth is None else grid_linewidth)
        self._grid_alpha = (rcParams['grid.alpha']
                            if grid_alpha is None else grid_alpha)

        self._grid_kw = {k[5:]: v for k, v in kw.items()}

        self.apply_tickdir(tickdir)

        self.tick1line = self._get_tick1line()
        self.tick2line = self._get_tick2line()
        self.gridline = self._get_gridline()

        self.label1 = self._get_text1()
        self.label = self.label1  # legacy name
        self.label2 = self._get_text2()

        self.gridOn = gridOn
        self.tick1On = tick1On
        self.tick2On = tick2On
        self.label1On = label1On
        self.label2On = label2On

        self.update_position(loc)

    def _set_labelrotation(self, labelrotation):
        if isinstance(labelrotation, six.string_types):
            mode = labelrotation
            angle = 0
        elif isinstance(labelrotation, (tuple, list)):
            mode, angle = labelrotation
        else:
            mode = 'default'
            angle = labelrotation
        if mode not in ('auto', 'default'):
            raise ValueError("Label rotation mode must be 'default' or "
                             "'auto', not '{}'.".format(mode))
        self._labelrotation = (mode, angle)

    def apply_tickdir(self, tickdir):
        """
        Calculate self._pad and self._tickmarkers
        """
        pass

    def get_tickdir(self):
        return self._tickdir

    def get_tick_padding(self):
        """
        Get the length of the tick outside of the axes.
        """
        padding = {
            'in': 0.0,
            'inout': 0.5,
            'out': 1.0
        }
        return self._size * padding[self._tickdir]

    def get_children(self):
        children = [self.tick1line, self.tick2line,
                    self.gridline, self.label1, self.label2]
        return children

    def set_clip_path(self, clippath, transform=None):
        artist.Artist.set_clip_path(self, clippath, transform)
        self.gridline.set_clip_path(clippath, transform)
        self.stale = True

    set_clip_path.__doc__ = artist.Artist.set_clip_path.__doc__

    def get_pad_pixels(self):
        return self.figure.dpi * self._base_pad / 72.0

    def contains(self, mouseevent):
        """
        Test whether the mouse event occurred in the Tick marks.

        This function always returns false.  It is more useful to test if the
        axis as a whole contains the mouse rather than the set of tick marks.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)
        return False, {}

    def set_pad(self, val):
        """
        Set the tick label pad in points

        ACCEPTS: float
        """
        self._apply_params(pad=val)
        self.stale = True

    def get_pad(self):
        'Get the value of the tick label pad in points'
        return self._base_pad

    def _get_text1(self):
        'Get the default Text 1 instance'
        pass

    def _get_text2(self):
        'Get the default Text 2 instance'
        pass

    def _get_tick1line(self):
        'Get the default line2D instance for tick1'
        pass

    def _get_tick2line(self):
        'Get the default line2D instance for tick2'
        pass

    def _get_gridline(self):
        'Get the default grid Line2d instance for this tick'
        pass

    def get_loc(self):
        'Return the tick location (data coords) as a scalar'
        return self._loc

    @allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__)
        if self.gridOn:
            self.gridline.draw(renderer)
        if self.tick1On:
            self.tick1line.draw(renderer)
        if self.tick2On:
            self.tick2line.draw(renderer)

        if self.label1On:
            self.label1.draw(renderer)
        if self.label2On:
            self.label2.draw(renderer)
        renderer.close_group(self.__name__)

        self.stale = False

    def set_label1(self, s):
        """
        Set the text of ticklabel

        ACCEPTS: str
        """
        self.label1.set_text(s)
        self.stale = True

    set_label = set_label1

    def set_label2(self, s):
        """
        Set the text of ticklabel2

        ACCEPTS: str
        """
        self.label2.set_text(s)
        self.stale = True

    def _set_artist_props(self, a):
        a.set_figure(self.figure)

    def get_view_interval(self):
        'return the view Interval instance for the axis this tick is ticking'
        raise NotImplementedError('Derived must override')

    def _apply_params(self, **kw):
        switchkw = ['gridOn', 'tick1On', 'tick2On', 'label1On', 'label2On']
        switches = [k for k in kw if k in switchkw]
        for k in switches:
            setattr(self, k, kw.pop(k))
        newmarker = [k for k in kw if k in ['size', 'width', 'pad', 'tickdir']]
        if newmarker:
            self._size = kw.pop('size', self._size)
            # Width could be handled outside this block, but it is
            # convenient to leave it here.
            self._width = kw.pop('width', self._width)
            self._base_pad = kw.pop('pad', self._base_pad)
            # apply_tickdir uses _size and _base_pad to make _pad,
            # and also makes _tickmarkers.
            self.apply_tickdir(kw.pop('tickdir', self._tickdir))
            self.tick1line.set_marker(self._tickmarkers[0])
            self.tick2line.set_marker(self._tickmarkers[1])
            for line in (self.tick1line, self.tick2line):
                line.set_markersize(self._size)
                line.set_markeredgewidth(self._width)
            # _get_text1_transform uses _pad from apply_tickdir.
            trans = self._get_text1_transform()[0]
            self.label1.set_transform(trans)
            trans = self._get_text2_transform()[0]
            self.label2.set_transform(trans)
        tick_kw = {k: v for k, v in six.iteritems(kw)
                   if k in ['color', 'zorder']}
        if tick_kw:
            self.tick1line.set(**tick_kw)
            self.tick2line.set(**tick_kw)
            for k, v in six.iteritems(tick_kw):
                setattr(self, '_' + k, v)

        if 'labelrotation' in kw:
            self._set_labelrotation(kw.pop('labelrotation'))
            self.label1.set(rotation=self._labelrotation[1])
            self.label2.set(rotation=self._labelrotation[1])

        label_list = [k for k in six.iteritems(kw)
                      if k[0] in ['labelsize', 'labelcolor']]
        if label_list:
            label_kw = {k[5:]: v for k, v in label_list}
            self.label1.set(**label_kw)
            self.label2.set(**label_kw)
            for k, v in six.iteritems(label_kw):
                # for labelsize the text objects covert str ('small')
                # -> points. grab the integer from the `Text` object
                # instead of saving the string representation
                v = getattr(self.label1, 'get_' + k)()
                setattr(self, '_label' + k, v)

        grid_list = [k for k in six.iteritems(kw)
                     if k[0] in _gridline_param_names]
        if grid_list:
            grid_kw = {k[5:]: v for k, v in grid_list}
            self.gridline.set(**grid_kw)
            for k, v in six.iteritems(grid_kw):
                setattr(self, '_grid_' + k, v)

    def update_position(self, loc):
        'Set the location of tick in data coords with scalar *loc*'
        raise NotImplementedError('Derived must override')

    def _get_text1_transform(self):
        raise NotImplementedError('Derived must override')

    def _get_text2_transform(self):
        raise NotImplementedError('Derived must override')


class XTick(Tick):
    """
    Contains all the Artists needed to make an x tick - the tick line,
    the label text and the grid line
    """
    __name__ = 'xtick'

    def _get_text1_transform(self):
        return self.axes.get_xaxis_text1_transform(self._pad)

    def _get_text2_transform(self):
        return self.axes.get_xaxis_text2_transform(self._pad)

    def apply_tickdir(self, tickdir):
        if tickdir is None:
            tickdir = rcParams['%s.direction' % self._name]
        self._tickdir = tickdir

        if self._tickdir == 'in':
            self._tickmarkers = (mlines.TICKUP, mlines.TICKDOWN)
        elif self._tickdir == 'inout':
            self._tickmarkers = ('|', '|')
        else:
            self._tickmarkers = (mlines.TICKDOWN, mlines.TICKUP)
        self._pad = self._base_pad + self.get_tick_padding()
        self.stale = True

    def _get_text1(self):
        'Get the default Text instance'
        # the y loc is 3 points below the min of y axis
        # get the affine as an a,b,c,d,tx,ty list
        # x in data coords, y in axes coords
        trans, vert, horiz = self._get_text1_transform()
        t = mtext.Text(
            x=0, y=0,
            fontproperties=font_manager.FontProperties(size=self._labelsize),
            color=self._labelcolor,
            verticalalignment=vert,
            horizontalalignment=horiz,
            )
        t.set_transform(trans)
        self._set_artist_props(t)
        return t

    def _get_text2(self):

        'Get the default Text 2 instance'
        # x in data coords, y in axes coords
        trans, vert, horiz = self._get_text2_transform()
        t = mtext.Text(
            x=0, y=1,
            fontproperties=font_manager.FontProperties(size=self._labelsize),
            color=self._labelcolor,
            verticalalignment=vert,
            horizontalalignment=horiz,
            )
        t.set_transform(trans)
        self._set_artist_props(t)
        return t

    def _get_tick1line(self):
        'Get the default line2D instance'
        # x in data coords, y in axes coords
        l = mlines.Line2D(xdata=(0,), ydata=(0,), color=self._color,
                          linestyle='None', marker=self._tickmarkers[0],
                          markersize=self._size,
                          markeredgewidth=self._width, zorder=self._zorder)
        l.set_transform(self.axes.get_xaxis_transform(which='tick1'))
        self._set_artist_props(l)
        return l

    def _get_tick2line(self):
        'Get the default line2D instance'
        # x in data coords, y in axes coords
        l = mlines.Line2D(xdata=(0,), ydata=(1,),
                          color=self._color,
                          linestyle='None',
                          marker=self._tickmarkers[1],
                          markersize=self._size,
                          markeredgewidth=self._width,
                          zorder=self._zorder)

        l.set_transform(self.axes.get_xaxis_transform(which='tick2'))
        self._set_artist_props(l)
        return l

    def _get_gridline(self):
        'Get the default line2D instance'
        # x in data coords, y in axes coords
        l = mlines.Line2D(xdata=(0.0, 0.0), ydata=(0, 1.0),
                          color=self._grid_color,
                          linestyle=self._grid_linestyle,
                          linewidth=self._grid_linewidth,
                          alpha=self._grid_alpha,
                          markersize=0,
                          **self._grid_kw)
        l.set_transform(self.axes.get_xaxis_transform(which='grid'))
        l.get_path()._interpolation_steps = GRIDLINE_INTERPOLATION_STEPS
        self._set_artist_props(l)

        return l

    def update_position(self, loc):
        'Set the location of tick in data coords with scalar *loc*'
        if self.tick1On:
            self.tick1line.set_xdata((loc,))
        if self.tick2On:
            self.tick2line.set_xdata((loc,))
        if self.gridOn:
            self.gridline.set_xdata((loc,))
        if self.label1On:
            self.label1.set_x(loc)
        if self.label2On:
            self.label2.set_x(loc)

        self._loc = loc
        self.stale = True

    def get_view_interval(self):
        'return the Interval instance for this axis view limits'
        return self.axes.viewLim.intervalx


class YTick(Tick):
    """
    Contains all the Artists needed to make a Y tick - the tick line,
    the label text and the grid line
    """
    __name__ = 'ytick'

    def _get_text1_transform(self):
        return self.axes.get_yaxis_text1_transform(self._pad)

    def _get_text2_transform(self):
        return self.axes.get_yaxis_text2_transform(self._pad)

    def apply_tickdir(self, tickdir):
        if tickdir is None:
            tickdir = rcParams['%s.direction' % self._name]
        self._tickdir = tickdir

        if self._tickdir == 'in':
            self._tickmarkers = (mlines.TICKRIGHT, mlines.TICKLEFT)
        elif self._tickdir == 'inout':
            self._tickmarkers = ('_', '_')
        else:
            self._tickmarkers = (mlines.TICKLEFT, mlines.TICKRIGHT)
        self._pad = self._base_pad + self.get_tick_padding()
        self.stale = True

    # how far from the y axis line the right of the ticklabel are
    def _get_text1(self):
        'Get the default Text instance'
        # x in axes coords, y in data coords
        trans, vert, horiz = self._get_text1_transform()
        t = mtext.Text(
            x=0, y=0,
            fontproperties=font_manager.FontProperties(size=self._labelsize),
            color=self._labelcolor,
            verticalalignment=vert,
            horizontalalignment=horiz,
            )
        t.set_transform(trans)
        self._set_artist_props(t)
        return t

    def _get_text2(self):
        'Get the default Text instance'
        # x in axes coords, y in data coords
        trans, vert, horiz = self._get_text2_transform()
        t = mtext.Text(
            x=1, y=0,
            fontproperties=font_manager.FontProperties(size=self._labelsize),
            color=self._labelcolor,
            verticalalignment=vert,
            horizontalalignment=horiz,
            )
        t.set_transform(trans)
        self._set_artist_props(t)
        return t

    def _get_tick1line(self):
        'Get the default line2D instance'
        # x in axes coords, y in data coords

        l = mlines.Line2D((0,), (0,),
                          color=self._color,
                          marker=self._tickmarkers[0],
                          linestyle='None',
                          markersize=self._size,
                          markeredgewidth=self._width,
                          zorder=self._zorder)
        l.set_transform(self.axes.get_yaxis_transform(which='tick1'))
        self._set_artist_props(l)
        return l

    def _get_tick2line(self):
        'Get the default line2D instance'
        # x in axes coords, y in data coords
        l = mlines.Line2D((1,), (0,),
                          color=self._color,
                          marker=self._tickmarkers[1],
                          linestyle='None',
                          markersize=self._size,
                          markeredgewidth=self._width,
                          zorder=self._zorder)
        l.set_transform(self.axes.get_yaxis_transform(which='tick2'))
        self._set_artist_props(l)
        return l

    def _get_gridline(self):
        'Get the default line2D instance'
        # x in axes coords, y in data coords
        l = mlines.Line2D(xdata=(0, 1), ydata=(0, 0),
                          color=self._grid_color,
                          linestyle=self._grid_linestyle,
                          linewidth=self._grid_linewidth,
                          alpha=self._grid_alpha,
                          markersize=0,
                          **self._grid_kw)
        l.set_transform(self.axes.get_yaxis_transform(which='grid'))
        l.get_path()._interpolation_steps = GRIDLINE_INTERPOLATION_STEPS
        self._set_artist_props(l)
        return l

    def update_position(self, loc):
        'Set the location of tick in data coords with scalar *loc*'
        if self.tick1On:
            self.tick1line.set_ydata((loc,))
        if self.tick2On:
            self.tick2line.set_ydata((loc,))
        if self.gridOn:
            self.gridline.set_ydata((loc,))
        if self.label1On:
            self.label1.set_y(loc)
        if self.label2On:
            self.label2.set_y(loc)

        self._loc = loc
        self.stale = True

    def get_view_interval(self):
        'return the Interval instance for this axis view limits'
        return self.axes.viewLim.intervaly


class Ticker(object):
    locator = None
    formatter = None


class _LazyTickList(object):
    """
    A descriptor for lazy instantiation of tick lists.

    See comment above definition of the ``majorTicks`` and ``minorTicks``
    attributes.
    """

    def __init__(self, major):
        self._major = major

    def __get__(self, instance, cls):
        if instance is None:
            return self
        else:
            # instance._get_tick() can itself try to access the majorTicks
            # attribute (e.g. in certain projection classes which override
            # e.g. get_xaxis_text1_transform).  In order to avoid infinite
            # recursion, first set the majorTicks on the instance to an empty
            # list, then create the tick and append it.
            if self._major:
                instance.majorTicks = []
                tick = instance._get_tick(major=True)
                instance.majorTicks.append(tick)
                return instance.majorTicks
            else:
                instance.minorTicks = []
                tick = instance._get_tick(major=False)
                instance.minorTicks.append(tick)
                return instance.minorTicks


class Axis(artist.Artist):
    """
    Public attributes

    * :attr:`axes.transData` - transform data coords to display coords
    * :attr:`axes.transAxes` - transform axis coords to display coords
    * :attr:`labelpad` - number of points between the axis and its label
    """
    OFFSETTEXTPAD = 3

    def __str__(self):
        return self.__class__.__name__ \
            + "(%f,%f)" % tuple(self.axes.transAxes.transform_point((0, 0)))

    def __init__(self, axes, pickradius=15):
        """
        Init the axis with the parent Axes instance
        """
        artist.Artist.__init__(self)
        self.set_figure(axes.figure)

        self.isDefault_label = True

        self.axes = axes
        self.major = Ticker()
        self.minor = Ticker()
        self.callbacks = cbook.CallbackRegistry()

        self._autolabelpos = True
        self._smart_bounds = False

        self.label = self._get_label()
        self.labelpad = rcParams['axes.labelpad']
        self.offsetText = self._get_offset_text()

        self.pickradius = pickradius

        # Initialize here for testing; later add API
        self._major_tick_kw = dict()
        self._minor_tick_kw = dict()

        self.cla()
        self._set_scale('linear')

    # During initialization, Axis objects often create ticks that are later
    # unused; this turns out to be a very slow step.  Instead, use a custom
    # descriptor to make the tick lists lazy and instantiate them as needed.
    majorTicks = _LazyTickList(major=True)
    minorTicks = _LazyTickList(major=False)

    def set_label_coords(self, x, y, transform=None):
        """
        Set the coordinates of the label.  By default, the x
        coordinate of the y label is determined by the tick label
        bounding boxes, but this can lead to poor alignment of
        multiple ylabels if there are multiple axes.  Ditto for the y
        coordinate of the x label.

        You can also specify the coordinate system of the label with
        the transform.  If None, the default coordinate system will be
        the axes coordinate system (0,0) is (left,bottom), (0.5, 0.5)
        is middle, etc

        """

        self._autolabelpos = False
        if transform is None:
            transform = self.axes.transAxes

        self.label.set_transform(transform)
        self.label.set_position((x, y))
        self.stale = True

    def get_transform(self):
        return self._scale.get_transform()

    def get_scale(self):
        return self._scale.name

    def _set_scale(self, value, **kwargs):
        self._scale = mscale.scale_factory(value, self, **kwargs)
        self._scale.set_default_locators_and_formatters(self)

        self.isDefault_majloc = True
        self.isDefault_minloc = True
        self.isDefault_majfmt = True
        self.isDefault_minfmt = True

    def limit_range_for_scale(self, vmin, vmax):
        return self._scale.limit_range_for_scale(vmin, vmax, self.get_minpos())

    @property
    @cbook.deprecated("2.2.0")
    def unit_data(self):
        return self.units

    @unit_data.setter
    @cbook.deprecated("2.2.0")
    def unit_data(self, unit_data):
        self.set_units(unit_data)

    def get_children(self):
        children = [self.label, self.offsetText]
        majorticks = self.get_major_ticks()
        minorticks = self.get_minor_ticks()

        children.extend(majorticks)
        children.extend(minorticks)
        return children

    def cla(self):
        'clear the current axis'

        self.label.set_text('')  # self.set_label_text would change isDefault_

        self._set_scale('linear')

        # Clear the callback registry for this axis, or it may "leak"
        self.callbacks = cbook.CallbackRegistry()

        # whether the grids are on
        self._gridOnMajor = (rcParams['axes.grid'] and
                             rcParams['axes.grid.which'] in ('both', 'major'))
        self._gridOnMinor = (rcParams['axes.grid'] and
                             rcParams['axes.grid.which'] in ('both', 'minor'))

        self.reset_ticks()

        self.converter = None
        self.units = None
        self.set_units(None)
        self.stale = True

    def reset_ticks(self):
        """
        Re-initialize the major and minor Tick lists.

        Each list starts with a single fresh Tick.
        """
        # Restore the lazy tick lists.
        try:
            del self.majorTicks
        except AttributeError:
            pass
        try:
            del self.minorTicks
        except AttributeError:
            pass
        try:
            self.set_clip_path(self.axes.patch)
        except AttributeError:
            pass

    def set_tick_params(self, which='major', reset=False, **kw):
        """
        Set appearance parameters for ticks, ticklabels, and gridlines.

        For documentation of keyword arguments, see
        :meth:`matplotlib.axes.Axes.tick_params`.
        """
        dicts = []
        if which == 'major' or which == 'both':
            dicts.append(self._major_tick_kw)
        if which == 'minor' or which == 'both':
            dicts.append(self._minor_tick_kw)
        kwtrans = self._translate_tick_kw(kw, to_init_kw=True)
        for d in dicts:
            if reset:
                d.clear()
            d.update(kwtrans)

        if reset:
            self.reset_ticks()
        else:
            if which == 'major' or which == 'both':
                for tick in self.majorTicks:
                    tick._apply_params(**self._major_tick_kw)
            if which == 'minor' or which == 'both':
                for tick in self.minorTicks:
                    tick._apply_params(**self._minor_tick_kw)
            if 'labelcolor' in kwtrans:
                self.offsetText.set_color(kwtrans['labelcolor'])
        self.stale = True

    @staticmethod
    def _translate_tick_kw(kw, to_init_kw=True):
        # The following lists may be moved to a more
        # accessible location.
        kwkeys0 = ['size', 'width', 'color', 'tickdir', 'pad',
                   'labelsize', 'labelcolor', 'zorder', 'gridOn',
                   'tick1On', 'tick2On', 'label1On', 'label2On']
        kwkeys1 = ['length', 'direction', 'left', 'bottom', 'right', 'top',
                   'labelleft', 'labelbottom', 'labelright', 'labeltop',
                   'labelrotation']
        kwkeys2 = _gridline_param_names
        kwkeys = kwkeys0 + kwkeys1 + kwkeys2
        kwtrans = dict()
        if to_init_kw:
            if 'length' in kw:
                kwtrans['size'] = kw.pop('length')
            if 'direction' in kw:
                kwtrans['tickdir'] = kw.pop('direction')
            if 'rotation' in kw:
                kwtrans['labelrotation'] = kw.pop('rotation')
            if 'left' in kw:
                kwtrans['tick1On'] = _string_to_bool(kw.pop('left'))
            if 'bottom' in kw:
                kwtrans['tick1On'] = _string_to_bool(kw.pop('bottom'))
            if 'right' in kw:
                kwtrans['tick2On'] = _string_to_bool(kw.pop('right'))
            if 'top' in kw:
                kwtrans['tick2On'] = _string_to_bool(kw.pop('top'))

            if 'labelleft' in kw:
                kwtrans['label1On'] = _string_to_bool(kw.pop('labelleft'))
            if 'labelbottom' in kw:
                kwtrans['label1On'] = _string_to_bool(kw.pop('labelbottom'))
            if 'labelright' in kw:
                kwtrans['label2On'] = _string_to_bool(kw.pop('labelright'))
            if 'labeltop' in kw:
                kwtrans['label2On'] = _string_to_bool(kw.pop('labeltop'))
            if 'colors' in kw:
                c = kw.pop('colors')
                kwtrans['color'] = c
                kwtrans['labelcolor'] = c
            # Maybe move the checking up to the caller of this method.
            for key in kw:
                if key not in kwkeys:
                    raise ValueError(
                        "keyword %s is not recognized; valid keywords are %s"
                        % (key, kwkeys))
            kwtrans.update(kw)
        else:
            raise NotImplementedError("Inverse translation is deferred")
        return kwtrans

    def set_clip_path(self, clippath, transform=None):
        artist.Artist.set_clip_path(self, clippath, transform)
        for child in self.majorTicks + self.minorTicks:
            child.set_clip_path(clippath, transform)
        self.stale = True

    def get_view_interval(self):
        'return the Interval instance for this axis view limits'
        raise NotImplementedError('Derived must override')

    def set_view_interval(self, vmin, vmax, ignore=False):
        raise NotImplementedError('Derived must override')

    def get_data_interval(self):
        'return the Interval instance for this axis data limits'
        raise NotImplementedError('Derived must override')

    def set_data_interval(self):
        '''set the axis data limits'''
        raise NotImplementedError('Derived must override')

    def set_default_intervals(self):
        '''set the default limits for the axis data and view interval if they
        are not mutated'''

        # this is mainly in support of custom object plotting.  For
        # example, if someone passes in a datetime object, we do not
        # know automagically how to set the default min/max of the
        # data and view limits.  The unit conversion AxisInfo
        # interface provides a hook for custom types to register
        # default limits through the AxisInfo.default_limits
        # attribute, and the derived code below will check for that
        # and use it if is available (else just use 0..1)
        pass

    def _set_artist_props(self, a):
        if a is None:
            return
        a.set_figure(self.figure)

    def iter_ticks(self):
        """
        Iterate through all of the major and minor ticks.
        """
        majorLocs = self.major.locator()
        majorTicks = self.get_major_ticks(len(majorLocs))
        self.major.formatter.set_locs(majorLocs)
        majorLabels = [self.major.formatter(val, i)
                       for i, val in enumerate(majorLocs)]

        minorLocs = self.minor.locator()
        minorTicks = self.get_minor_ticks(len(minorLocs))
        self.minor.formatter.set_locs(minorLocs)
        minorLabels = [self.minor.formatter(val, i)
                       for i, val in enumerate(minorLocs)]

        major_minor = [
            (majorTicks, majorLocs, majorLabels),
            (minorTicks, minorLocs, minorLabels)]

        for group in major_minor:
            for tick in zip(*group):
                yield tick

    def get_ticklabel_extents(self, renderer):
        """
        Get the extents of the tick labels on either side
        of the axes.
        """

        ticks_to_draw = self._update_ticks(renderer)
        ticklabelBoxes, ticklabelBoxes2 = self._get_tick_bboxes(ticks_to_draw,
                                                                renderer)

        if len(ticklabelBoxes):
            bbox = mtransforms.Bbox.union(ticklabelBoxes)
        else:
            bbox = mtransforms.Bbox.from_extents(0, 0, 0, 0)
        if len(ticklabelBoxes2):
            bbox2 = mtransforms.Bbox.union(ticklabelBoxes2)
        else:
            bbox2 = mtransforms.Bbox.from_extents(0, 0, 0, 0)
        return bbox, bbox2

    def set_smart_bounds(self, value):
        """set the axis to have smart bounds"""
        self._smart_bounds = value
        self.stale = True

    def get_smart_bounds(self):
        """get whether the axis has smart bounds"""
        return self._smart_bounds

    def _update_ticks(self, renderer):
        """
        Update ticks (position and labels) using the current data
        interval of the axes. Returns a list of ticks that will be
        drawn.
        """

        interval = self.get_view_interval()
        tick_tups = list(self.iter_ticks())  # iter_ticks calls the locator
        if self._smart_bounds and tick_tups:
            # handle inverted limits
            view_low, view_high = sorted(interval)
            data_low, data_high = sorted(self.get_data_interval())
            locs = np.sort([ti[1] for ti in tick_tups])
            if data_low <= view_low:
                # data extends beyond view, take view as limit
                ilow = view_low
            else:
                # data stops within view, take best tick
                good_locs = locs[locs <= data_low]
                if len(good_locs):
                    # last tick prior or equal to first data point
                    ilow = good_locs[-1]
                else:
                    # No ticks (why not?), take first tick
                    ilow = locs[0]
            if data_high >= view_high:
                # data extends beyond view, take view as limit
                ihigh = view_high
            else:
                # data stops within view, take best tick
                good_locs = locs[locs >= data_high]
                if len(good_locs):
                    # first tick after or equal to last data point
                    ihigh = good_locs[0]
                else:
                    # No ticks (why not?), take last tick
                    ihigh = locs[-1]
            tick_tups = [ti for ti in tick_tups if ilow <= ti[1] <= ihigh]

        # so that we don't lose ticks on the end, expand out the interval ever
        # so slightly.  The "ever so slightly" is defined to be the width of a
        # half of a pixel.  We don't want to draw a tick that even one pixel
        # outside of the defined axis interval.
        if interval[0] <= interval[1]:
            interval_expanded = interval
        else:
            interval_expanded = interval[1], interval[0]

        if hasattr(self, '_get_pixel_distance_along_axis'):
            # normally, one does not want to catch all exceptions that
            # could possibly happen, but it is not clear exactly what
            # exceptions might arise from a user's projection (their
            # rendition of the Axis object).  So, we catch all, with
            # the idea that one would rather potentially lose a tick
            # from one side of the axis or another, rather than see a
            # stack trace.
            # We also catch users warnings here. These are the result of
            # invalid numpy calculations that may be the result of out of
            # bounds on axis with finite allowed intervals such as geo
            # projections i.e. Mollweide.
            with np.errstate(invalid='ignore'):
                try:
                    ds1 = self._get_pixel_distance_along_axis(
                        interval_expanded[0], -0.5)
                except:
                    warnings.warn("Unable to find pixel distance along axis "
                                  "for interval padding of ticks; assuming no "
                                  "interval padding needed.")
                    ds1 = 0.0
                if np.isnan(ds1):
                    ds1 = 0.0
                try:
                    ds2 = self._get_pixel_distance_along_axis(
                        interval_expanded[1], +0.5)
                except:
                    warnings.warn("Unable to find pixel distance along axis "
                                  "for interval padding of ticks; assuming no "
                                  "interval padding needed.")
                    ds2 = 0.0
                if np.isnan(ds2):
                    ds2 = 0.0
            interval_expanded = (interval_expanded[0] - ds1,
                                 interval_expanded[1] + ds2)

        ticks_to_draw = []
        for tick, loc, label in tick_tups:
            if tick is None:
                continue
            # NB: always update labels and position to avoid issues like #9397
            tick.update_position(loc)
            tick.set_label1(label)
            tick.set_label2(label)
            if not mtransforms.interval_contains(interval_expanded, loc):
                continue
            ticks_to_draw.append(tick)

        return ticks_to_draw

    def _get_tick_bboxes(self, ticks, renderer):
        """
        Given the list of ticks, return two lists of bboxes. One for
        tick lable1's and another for tick label2's.
        """

        ticklabelBoxes = []
        ticklabelBoxes2 = []

        for tick in ticks:
            if tick.label1On and tick.label1.get_visible():
                extent = tick.label1.get_window_extent(renderer)
                ticklabelBoxes.append(extent)
            if tick.label2On and tick.label2.get_visible():
                extent = tick.label2.get_window_extent(renderer)
                ticklabelBoxes2.append(extent)
        return ticklabelBoxes, ticklabelBoxes2

    def get_tightbbox(self, renderer):
        """
        Return a bounding box that encloses the axis. It only accounts
        tick labels, axis label, and offsetText.
        """
        if not self.get_visible():
            return

        ticks_to_draw = self._update_ticks(renderer)

        self._update_label_position(renderer)

        # go back to just this axis's tick labels
        ticklabelBoxes, ticklabelBoxes2 = self._get_tick_bboxes(
                    ticks_to_draw, renderer)

        self._update_offset_text_position(ticklabelBoxes, ticklabelBoxes2)
        self.offsetText.set_text(self.major.formatter.get_offset())

        bb = []

        for a in [self.label, self.offsetText]:
            if a.get_visible():
                bb.append(a.get_window_extent(renderer))

        bb.extend(ticklabelBoxes)
        bb.extend(ticklabelBoxes2)

        bb = [b for b in bb if b.width != 0 or b.height != 0]
        if bb:
            _bbox = mtransforms.Bbox.union(bb)
            return _bbox
        else:
            return None

    def get_tick_padding(self):
        values = []
        if len(self.majorTicks):
            values.append(self.majorTicks[0].get_tick_padding())
        if len(self.minorTicks):
            values.append(self.minorTicks[0].get_tick_padding())
        if len(values):
            return max(values)
        return 0.0

    @allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        'Draw the axis lines, grid lines, tick lines and labels'

        if not self.get_visible():
            return
        renderer.open_group(__name__)

        ticks_to_draw = self._update_ticks(renderer)
        ticklabelBoxes, ticklabelBoxes2 = self._get_tick_bboxes(ticks_to_draw,
                                                                renderer)

        for tick in ticks_to_draw:
            tick.draw(renderer)

        # scale up the axis label box to also find the neighbors, not
        # just the tick labels that actually overlap note we need a
        # *copy* of the axis label box because we don't wan't to scale
        # the actual bbox

        self._update_label_position(renderer)

        self.label.draw(renderer)

        self._update_offset_text_position(ticklabelBoxes, ticklabelBoxes2)
        self.offsetText.set_text(self.major.formatter.get_offset())
        self.offsetText.draw(renderer)

        if 0:  # draw the bounding boxes around the text for debug
            for tick in self.majorTicks:
                label = tick.label1
                mpatches.bbox_artist(label, renderer)
            mpatches.bbox_artist(self.label, renderer)

        renderer.close_group(__name__)
        self.stale = False

    def _get_label(self):
        raise NotImplementedError('Derived must override')

    def _get_offset_text(self):
        raise NotImplementedError('Derived must override')

    def get_gridlines(self):
        'Return the grid lines as a list of Line2D instance'
        ticks = self.get_major_ticks()
        return cbook.silent_list('Line2D gridline',
                                 [tick.gridline for tick in ticks])

    def get_label(self):
        'Return the axis label as a Text instance'
        return self.label

    def get_offset_text(self):
        'Return the axis offsetText as a Text instance'
        return self.offsetText

    def get_pickradius(self):
        'Return the depth of the axis used by the picker'
        return self.pickradius

    def get_majorticklabels(self):
        'Return a list of Text instances for the major ticklabels'
        ticks = self.get_major_ticks()
        labels1 = [tick.label1 for tick in ticks if tick.label1On]
        labels2 = [tick.label2 for tick in ticks if tick.label2On]
        return cbook.silent_list('Text major ticklabel', labels1 + labels2)

    def get_minorticklabels(self):
        'Return a list of Text instances for the minor ticklabels'
        ticks = self.get_minor_ticks()
        labels1 = [tick.label1 for tick in ticks if tick.label1On]
        labels2 = [tick.label2 for tick in ticks if tick.label2On]
        return cbook.silent_list('Text minor ticklabel', labels1 + labels2)

    def get_ticklabels(self, minor=False, which=None):
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

        if which is not None:
            if which == 'minor':
                return self.get_minorticklabels()
            elif which == 'major':
                return self.get_majorticklabels()
            elif which == 'both':
                return self.get_majorticklabels() + self.get_minorticklabels()
            else:
                raise ValueError("`which` must be one of ('minor', 'major', "
                                 "'both') not " + str(which))
        if minor:
            return self.get_minorticklabels()
        return self.get_majorticklabels()

    def get_majorticklines(self):
        'Return the major tick lines as a list of Line2D instances'
        lines = []
        ticks = self.get_major_ticks()
        for tick in ticks:
            lines.append(tick.tick1line)
            lines.append(tick.tick2line)
        return cbook.silent_list('Line2D ticklines', lines)

    def get_minorticklines(self):
        'Return the minor tick lines as a list of Line2D instances'
        lines = []
        ticks = self.get_minor_ticks()
        for tick in ticks:
            lines.append(tick.tick1line)
            lines.append(tick.tick2line)
        return cbook.silent_list('Line2D ticklines', lines)

    def get_ticklines(self, minor=False):
        'Return the tick lines as a list of Line2D instances'
        if minor:
            return self.get_minorticklines()
        return self.get_majorticklines()

    def get_majorticklocs(self):
        "Get the major tick locations in data coordinates as a numpy array"
        return self.major.locator()

    def get_minorticklocs(self):
        "Get the minor tick locations in data coordinates as a numpy array"
        return self.minor.locator()

    def get_ticklocs(self, minor=False):
        "Get the tick locations in data coordinates as a numpy array"
        if minor:
            return self.minor.locator()
        return self.major.locator()

    def get_ticks_direction(self, minor=False):
        """
        Get the tick directions as a numpy array

        Parameters
        ----------
        minor : boolean
            True to return the minor tick directions,
            False to return the major tick directions,
            Default is False

        Returns
        -------
        numpy array of tick directions
        """
        if minor:
            return np.array(
                [tick._tickdir for tick in self.get_minor_ticks()])
        else:
            return np.array(
                [tick._tickdir for tick in self.get_major_ticks()])

    def _get_tick(self, major):
        'return the default tick instance'
        raise NotImplementedError('derived must override')

    def _copy_tick_props(self, src, dest):
        'Copy the props from src tick to dest tick'
        if src is None or dest is None:
            return
        dest.label1.update_from(src.label1)
        dest.label2.update_from(src.label2)

        dest.tick1line.update_from(src.tick1line)
        dest.tick2line.update_from(src.tick2line)
        dest.gridline.update_from(src.gridline)

        dest.tick1On = src.tick1On
        dest.tick2On = src.tick2On
        dest.label1On = src.label1On
        dest.label2On = src.label2On

    def get_label_text(self):
        'Get the text of the label'
        return self.label.get_text()

    def get_major_locator(self):
        'Get the locator of the major ticker'
        return self.major.locator

    def get_minor_locator(self):
        'Get the locator of the minor ticker'
        return self.minor.locator

    def get_major_formatter(self):
        'Get the formatter of the major ticker'
        return self.major.formatter

    def get_minor_formatter(self):
        'Get the formatter of the minor ticker'
        return self.minor.formatter

    def get_major_ticks(self, numticks=None):
        'get the tick instances; grow as necessary'
        if numticks is None:
            numticks = len(self.get_major_locator()())

        while len(self.majorTicks) < numticks:
            # update the new tick label properties from the old
            tick = self._get_tick(major=True)
            self.majorTicks.append(tick)
            if self._gridOnMajor:
                tick.gridOn = True
            self._copy_tick_props(self.majorTicks[0], tick)

        return self.majorTicks[:numticks]

    def get_minor_ticks(self, numticks=None):
        'get the minor tick instances; grow as necessary'
        if numticks is None:
            numticks = len(self.get_minor_locator()())

        while len(self.minorTicks) < numticks:
            # update the new tick label properties from the old
            tick = self._get_tick(major=False)
            self.minorTicks.append(tick)
            if self._gridOnMinor:
                tick.gridOn = True
            self._copy_tick_props(self.minorTicks[0], tick)

        return self.minorTicks[:numticks]

    def grid(self, b=None, which='major', **kwargs):
        """
        Set the axis grid on or off; b is a boolean. Use *which* =
        'major' | 'minor' | 'both' to set the grid for major or minor ticks.

        If *b* is *None* and len(kwargs)==0, toggle the grid state.  If
        *kwargs* are supplied, it is assumed you want the grid on and *b*
        will be set to True.

        *kwargs* are used to set the line properties of the grids, e.g.,

          xax.grid(color='r', linestyle='-', linewidth=2)
        """
        if len(kwargs):
            b = True
        which = which.lower()
        gridkw = {'grid_' + item[0]: item[1] for item in kwargs.items()}
        if which in ['minor', 'both']:
            if b is None:
                self._gridOnMinor = not self._gridOnMinor
            else:
                self._gridOnMinor = b
            self.set_tick_params(which='minor', gridOn=self._gridOnMinor,
                                 **gridkw)
        if which in ['major', 'both']:
            if b is None:
                self._gridOnMajor = not self._gridOnMajor
            else:
                self._gridOnMajor = b
            self.set_tick_params(which='major', gridOn=self._gridOnMajor,
                                 **gridkw)
        self.stale = True

    def update_units(self, data):
        """
        introspect *data* for units converter and update the
        axis.converter instance if necessary. Return *True*
        if *data* is registered for unit conversion.
        """

        converter = munits.registry.get_converter(data)
        if converter is None:
            return False

        neednew = self.converter != converter
        self.converter = converter
        default = self.converter.default_units(data, self)
        if default is not None and self.units is None:
            self.set_units(default)

        if neednew:
            self._update_axisinfo()
        self.stale = True
        return True

    def _update_axisinfo(self):
        """
        check the axis converter for the stored units to see if the
        axis info needs to be updated
        """
        if self.converter is None:
            return

        info = self.converter.axisinfo(self.units, self)

        if info is None:
            return
        if info.majloc is not None and \
           self.major.locator != info.majloc and self.isDefault_majloc:
            self.set_major_locator(info.majloc)
            self.isDefault_majloc = True
        if info.minloc is not None and \
           self.minor.locator != info.minloc and self.isDefault_minloc:
            self.set_minor_locator(info.minloc)
            self.isDefault_minloc = True
        if info.majfmt is not None and \
           self.major.formatter != info.majfmt and self.isDefault_majfmt:
            self.set_major_formatter(info.majfmt)
            self.isDefault_majfmt = True
        if info.minfmt is not None and \
           self.minor.formatter != info.minfmt and self.isDefault_minfmt:
            self.set_minor_formatter(info.minfmt)
            self.isDefault_minfmt = True
        if info.label is not None and self.isDefault_label:
            self.set_label_text(info.label)
            self.isDefault_label = True

        self.set_default_intervals()

    def have_units(self):
        return self.converter is not None or self.units is not None

    def convert_units(self, x):
        # If x is already a number, doesn't need converting
        if munits.ConversionInterface.is_numlike(x):
            return x

        if self.converter is None:
            self.converter = munits.registry.get_converter(x)

        if self.converter is None:
            return x

        ret = self.converter.convert(x, self.units, self)
        return ret

    def set_units(self, u):
        """
        set the units for axis

        ACCEPTS: a units tag
        """
        pchanged = False
        if u is None:
            self.units = None
            pchanged = True
        else:
            if u != self.units:
                self.units = u
                pchanged = True
        if pchanged:
            self._update_axisinfo()
            self.callbacks.process('units')
            self.callbacks.process('units finalize')
        self.stale = True

    def get_units(self):
        'return the units for axis'
        return self.units

    def set_label_text(self, label, fontdict=None, **kwargs):
        """  Sets the text value of the axis label

        ACCEPTS: A string value for the label
        """
        self.isDefault_label = False
        self.label.set_text(label)
        if fontdict is not None:
            self.label.update(fontdict)
        self.label.update(kwargs)
        self.stale = True
        return self.label

    def set_major_formatter(self, formatter):
        """
        Set the formatter of the major ticker

        ACCEPTS: A :class:`~matplotlib.ticker.Formatter` instance
        """
        self.isDefault_majfmt = False
        self.major.formatter = formatter
        formatter.set_axis(self)
        self.stale = True

    def set_minor_formatter(self, formatter):
        """
        Set the formatter of the minor ticker

        ACCEPTS: A :class:`~matplotlib.ticker.Formatter` instance
        """
        self.isDefault_minfmt = False
        self.minor.formatter = formatter
        formatter.set_axis(self)
        self.stale = True

    def set_major_locator(self, locator):
        """
        Set the locator of the major ticker

        ACCEPTS: a :class:`~matplotlib.ticker.Locator` instance
        """
        self.isDefault_majloc = False
        self.major.locator = locator
        locator.set_axis(self)
        self.stale = True

    def set_minor_locator(self, locator):
        """
        Set the locator of the minor ticker

        ACCEPTS: a :class:`~matplotlib.ticker.Locator` instance
        """
        self.isDefault_minloc = False
        self.minor.locator = locator
        locator.set_axis(self)
        self.stale = True

    def set_pickradius(self, pickradius):
        """
        Set the depth of the axis used by the picker

        ACCEPTS: a distance in points
        """
        self.pickradius = pickradius

    def set_ticklabels(self, ticklabels, *args, **kwargs):
        """
        Set the text values of the tick labels. Return a list of Text
        instances.  Use *kwarg* *minor=True* to select minor ticks.
        All other kwargs are used to update the text object properties.
        As for get_ticklabels, label1 (left or bottom) is
        affected for a given tick only if its label1On attribute
        is True, and similarly for label2.  The list of returned
        label text objects consists of all such label1 objects followed
        by all such label2 objects.

        The input *ticklabels* is assumed to match the set of
        tick locations, regardless of the state of label1On and
        label2On.

        ACCEPTS: sequence of strings or Text objects
        """
        get_labels = []
        for t in ticklabels:
            # try calling get_text() to check whether it is Text object
            # if it is Text, get label content
            try:
                get_labels.append(t.get_text())
            # otherwise add the label to the list directly
            except AttributeError:
                get_labels.append(t)
        # replace the ticklabels list with the processed one
        ticklabels = get_labels

        minor = kwargs.pop('minor', False)
        if minor:
            self.set_minor_formatter(mticker.FixedFormatter(ticklabels))
            ticks = self.get_minor_ticks()
        else:
            self.set_major_formatter(mticker.FixedFormatter(ticklabels))
            ticks = self.get_major_ticks()
        ret = []
        for tick_label, tick in zip(ticklabels, ticks):
            # deal with label1
            tick.label1.set_text(tick_label)
            tick.label1.update(kwargs)
            # deal with label2
            tick.label2.set_text(tick_label)
            tick.label2.update(kwargs)
            # only return visible tick labels
            if tick.label1On:
                ret.append(tick.label1)
            if tick.label2On:
                ret.append(tick.label2)

        self.stale = True
        return ret

    def set_ticks(self, ticks, minor=False):
        """
        Set the locations of the tick marks from sequence ticks

        ACCEPTS: sequence of floats
        """
        # XXX if the user changes units, the information will be lost here
        ticks = self.convert_units(ticks)
        if len(ticks) > 1:
            xleft, xright = self.get_view_interval()
            if xright > xleft:
                self.set_view_interval(min(ticks), max(ticks))
            else:
                self.set_view_interval(max(ticks), min(ticks))
        if minor:
            self.set_minor_locator(mticker.FixedLocator(ticks))
            return self.get_minor_ticks(len(ticks))
        else:
            self.set_major_locator(mticker.FixedLocator(ticks))
            return self.get_major_ticks(len(ticks))

    def _get_tick_boxes_siblings(self, xdir, renderer):
        """
        Get the bounding boxes for this `.axis` and its siblings
        as set by `.Figure.align_xlabels` or  `.Figure.align_ylablels`.

        By default it just gets bboxes for self.
        """
        raise NotImplementedError('Derived must override')

    def _update_label_position(self, renderer):
        """
        Update the label position based on the bounding box enclosing
        all the ticklabels and axis spine
        """
        raise NotImplementedError('Derived must override')

    def _update_offset_text_position(self, bboxes, bboxes2):
        """
        Update the label position based on the sequence of bounding
        boxes of all the ticklabels
        """
        raise NotImplementedError('Derived must override')

    def pan(self, numsteps):
        'Pan *numsteps* (can be positive or negative)'
        self.major.locator.pan(numsteps)

    def zoom(self, direction):
        "Zoom in/out on axis; if *direction* is >0 zoom in, else zoom out"
        self.major.locator.zoom(direction)

    def axis_date(self, tz=None):
        """
        Sets up x-axis ticks and labels that treat the x data as dates.
        *tz* is a :class:`tzinfo` instance or a timezone string.
        This timezone is used to create date labels.
        """
        # By providing a sample datetime instance with the desired
        # timezone, the registered converter can be selected,
        # and the "units" attribute, which is the timezone, can
        # be set.
        import datetime
        if isinstance(tz, six.string_types):
            import pytz
            tz = pytz.timezone(tz)
        self.update_units(datetime.datetime(2009, 1, 1, 0, 0, 0, 0, tz))

    def get_tick_space(self):
        """
        Return the estimated number of ticks that can fit on the axis.
        """
        # Must be overridden in the subclass
        raise NotImplementedError()

    def get_label_position(self):
        """
        Return the label position (top or bottom)
        """
        return self.label_position

    def set_label_position(self, position):
        """
        Set the label position (top or bottom)

        ACCEPTS: [ 'top' | 'bottom' ]
        """
        raise NotImplementedError()

    def get_minpos(self):
        raise NotImplementedError()


class XAxis(Axis):
    __name__ = 'xaxis'
    axis_name = 'x'

    def contains(self, mouseevent):
        """Test whether the mouse event occurred in the x axis.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)

        x, y = mouseevent.x, mouseevent.y
        try:
            trans = self.axes.transAxes.inverted()
            xaxes, yaxes = trans.transform_point((x, y))
        except ValueError:
            return False, {}
        l, b = self.axes.transAxes.transform_point((0, 0))
        r, t = self.axes.transAxes.transform_point((1, 1))
        inaxis = xaxes >= 0 and xaxes <= 1 and (
            (y < b and y > b - self.pickradius) or
            (y > t and y < t + self.pickradius))
        return inaxis, {}

    def _get_tick(self, major):
        if major:
            tick_kw = self._major_tick_kw
        else:
            tick_kw = self._minor_tick_kw
        return XTick(self.axes, 0, '', major=major, **tick_kw)

    def _get_label(self):
        # x in axes coords, y in display coords (to be updated at draw
        # time by _update_label_positions)
        label = mtext.Text(x=0.5, y=0,
                           fontproperties=font_manager.FontProperties(
                               size=rcParams['axes.labelsize'],
                               weight=rcParams['axes.labelweight']),
                           color=rcParams['axes.labelcolor'],
                           verticalalignment='top',
                           horizontalalignment='center')

        label.set_transform(mtransforms.blended_transform_factory(
            self.axes.transAxes, mtransforms.IdentityTransform()))

        self._set_artist_props(label)
        self.label_position = 'bottom'
        return label

    def _get_offset_text(self):
        # x in axes coords, y in display coords (to be updated at draw time)
        offsetText = mtext.Text(x=1, y=0,
                                fontproperties=font_manager.FontProperties(
                                    size=rcParams['xtick.labelsize']),
                                color=rcParams['xtick.color'],
                                verticalalignment='top',
                                horizontalalignment='right')
        offsetText.set_transform(mtransforms.blended_transform_factory(
            self.axes.transAxes, mtransforms.IdentityTransform())
        )
        self._set_artist_props(offsetText)
        self.offset_text_position = 'bottom'
        return offsetText

    def _get_pixel_distance_along_axis(self, where, perturb):
        """
        Returns the amount, in data coordinates, that a single pixel
        corresponds to in the locality given by "where", which is also given
        in data coordinates, and is an x coordinate. "perturb" is the amount
        to perturb the pixel.  Usually +0.5 or -0.5.

        Implementing this routine for an axis is optional; if present, it will
        ensure that no ticks are lost due to round-off at the extreme ends of
        an axis.
        """

        # Note that this routine does not work for a polar axis, because of
        # the 1e-10 below.  To do things correctly, we need to use rmax
        # instead of 1e-10 for a polar axis.  But since we do not have that
        # kind of information at this point, we just don't try to pad anything
        # for the theta axis of a polar plot.
        if self.axes.name == 'polar':
            return 0.0

        #
        # first figure out the pixel location of the "where" point.  We use
        # 1e-10 for the y point, so that we remain compatible with log axes.

        # transformation from data coords to display coords
        trans = self.axes.transData
        # transformation from display coords to data coords
        transinv = trans.inverted()
        pix = trans.transform_point((where, 1e-10))
        # perturb the pixel
        ptp = transinv.transform_point((pix[0] + perturb, pix[1]))
        dx = abs(ptp[0] - where)

        return dx

    def set_label_position(self, position):
        """
        Set the label position (top or bottom)

        ACCEPTS: [ 'top' | 'bottom' ]
        """
        if position == 'top':
            self.label.set_verticalalignment('baseline')
        elif position == 'bottom':
            self.label.set_verticalalignment('top')
        else:
            raise ValueError("Position accepts only 'top' or 'bottom'")
        self.label_position = position
        self.stale = True

    def _get_tick_boxes_siblings(self, renderer):
        """
        Get the bounding boxes for this `.axis` and its siblings
        as set by `.Figure.align_xlabels` or  `.Figure.align_ylablels`.

        By default it just gets bboxes for self.
        """
        bboxes = []
        bboxes2 = []
        # get the Grouper that keeps track of x-label groups for this figure
        grp = self.figure._align_xlabel_grp
        # if we want to align labels from other axes:
        for nn, axx in enumerate(grp.get_siblings(self.axes)):
            ticks_to_draw = axx.xaxis._update_ticks(renderer)
            tlb, tlb2 = axx.xaxis._get_tick_bboxes(ticks_to_draw, renderer)
            bboxes.extend(tlb)
            bboxes2.extend(tlb2)
        return bboxes, bboxes2

    def _update_label_position(self, renderer):
        """
        Update the label position based on the bounding box enclosing
        all the ticklabels and axis spine
        """
        if not self._autolabelpos:
            return

        # get bounding boxes for this axis and any siblings
        # that have been set by `fig.align_xlabels()`
        bboxes, bboxes2 = self._get_tick_boxes_siblings(renderer=renderer)

        x, y = self.label.get_position()
        if self.label_position == 'bottom':
            try:
                spine = self.axes.spines['bottom']
                spinebbox = spine.get_transform().transform_path(
                    spine.get_path()).get_extents()
            except KeyError:
                # use axes if spine doesn't exist
                spinebbox = self.axes.bbox
            bbox = mtransforms.Bbox.union(bboxes + [spinebbox])
            bottom = bbox.y0

            self.label.set_position(
                (x, bottom - self.labelpad * self.figure.dpi / 72.0)
            )

        else:
            try:
                spine = self.axes.spines['top']
                spinebbox = spine.get_transform().transform_path(
                    spine.get_path()).get_extents()
            except KeyError:
                # use axes if spine doesn't exist
                spinebbox = self.axes.bbox
            bbox = mtransforms.Bbox.union(bboxes2 + [spinebbox])
            top = bbox.y1

            self.label.set_position(
                (x, top + self.labelpad * self.figure.dpi / 72.0)
            )

    def _update_offset_text_position(self, bboxes, bboxes2):
        """
        Update the offset_text position based on the sequence of bounding
        boxes of all the ticklabels
        """
        x, y = self.offsetText.get_position()
        if not len(bboxes):
            bottom = self.axes.bbox.ymin
        else:
            bbox = mtransforms.Bbox.union(bboxes)
            bottom = bbox.y0
        self.offsetText.set_position(
            (x, bottom - self.OFFSETTEXTPAD * self.figure.dpi / 72.0)
        )

    def get_text_heights(self, renderer):
        """
        Returns the amount of space one should reserve for text
        above and below the axes.  Returns a tuple (above, below)
        """
        bbox, bbox2 = self.get_ticklabel_extents(renderer)
        # MGDTODO: Need a better way to get the pad
        padPixels = self.majorTicks[0].get_pad_pixels()

        above = 0.0
        if bbox2.height:
            above += bbox2.height + padPixels
        below = 0.0
        if bbox.height:
            below += bbox.height + padPixels

        if self.get_label_position() == 'top':
            above += self.label.get_window_extent(renderer).height + padPixels
        else:
            below += self.label.get_window_extent(renderer).height + padPixels
        return above, below

    def set_ticks_position(self, position):
        """
        Set the ticks position (top, bottom, both, default or none)
        both sets the ticks to appear on both positions, but does not
        change the tick labels.  'default' resets the tick positions to
        the default: ticks on both positions, labels at bottom.  'none'
        can be used if you don't want any ticks. 'none' and 'both'
        affect only the ticks, not the labels.

        ACCEPTS: [ 'top' | 'bottom' | 'both' | 'default' | 'none' ]
        """
        if position == 'top':
            self.set_tick_params(which='both', top=True, labeltop=True,
                                 bottom=False, labelbottom=False)
        elif position == 'bottom':
            self.set_tick_params(which='both', top=False, labeltop=False,
                                 bottom=True, labelbottom=True)
        elif position == 'both':
            self.set_tick_params(which='both', top=True,
                                 bottom=True)
        elif position == 'none':
            self.set_tick_params(which='both', top=False,
                                 bottom=False)
        elif position == 'default':
            self.set_tick_params(which='both', top=True, labeltop=False,
                                 bottom=True, labelbottom=True)
        else:
            raise ValueError("invalid position: %s" % position)
        self.stale = True

    def tick_top(self):
        """
        Move ticks and ticklabels (if present) to the top of the axes.
        """
        label = True
        if 'label1On' in self._major_tick_kw:
            label = (self._major_tick_kw['label1On']
                     or self._major_tick_kw['label2On'])
        self.set_ticks_position('top')
        # if labels were turned off before this was called
        # leave them off
        self.set_tick_params(which='both', labeltop=label)

    def tick_bottom(self):
        """
        Move ticks and ticklabels (if present) to the bottom of the axes.
        """
        label = True
        if 'label1On' in self._major_tick_kw:
            label = (self._major_tick_kw['label1On']
                     or self._major_tick_kw['label2On'])
        self.set_ticks_position('bottom')
        # if labels were turned off before this was called
        # leave them off
        self.set_tick_params(which='both', labelbottom=label)

    def get_ticks_position(self):
        """
        Return the ticks position (top, bottom, default or unknown)
        """
        majt = self.majorTicks[0]
        mT = self.minorTicks[0]

        majorTop = ((not majt.tick1On) and majt.tick2On and
                    (not majt.label1On) and majt.label2On)
        minorTop = ((not mT.tick1On) and mT.tick2On and
                    (not mT.label1On) and mT.label2On)
        if majorTop and minorTop:
            return 'top'

        MajorBottom = (majt.tick1On and (not majt.tick2On) and
                       majt.label1On and (not majt.label2On))
        MinorBottom = (mT.tick1On and (not mT.tick2On) and
                       mT.label1On and (not mT.label2On))
        if MajorBottom and MinorBottom:
            return 'bottom'

        majorDefault = (majt.tick1On and majt.tick2On and
                        majt.label1On and (not majt.label2On))
        minorDefault = (mT.tick1On and mT.tick2On and
                        mT.label1On and (not mT.label2On))
        if majorDefault and minorDefault:
            return 'default'

        return 'unknown'

    def get_view_interval(self):
        'return the Interval instance for this axis view limits'
        return self.axes.viewLim.intervalx

    def set_view_interval(self, vmin, vmax, ignore=False):
        """
        If *ignore* is *False*, the order of vmin, vmax
        does not matter; the original axis orientation will
        be preserved. In addition, the view limits can be
        expanded, but will not be reduced.  This method is
        for mpl internal use; for normal use, see
        :meth:`~matplotlib.axes.Axes.set_xlim`.

        """
        if ignore:
            self.axes.viewLim.intervalx = vmin, vmax
        else:
            Vmin, Vmax = self.get_view_interval()
            if Vmin < Vmax:
                self.axes.viewLim.intervalx = (min(vmin, vmax, Vmin),
                                               max(vmin, vmax, Vmax))
            else:
                self.axes.viewLim.intervalx = (max(vmin, vmax, Vmin),
                                               min(vmin, vmax, Vmax))

    def get_minpos(self):
        return self.axes.dataLim.minposx

    def get_data_interval(self):
        'return the Interval instance for this axis data limits'
        return self.axes.dataLim.intervalx

    def set_data_interval(self, vmin, vmax, ignore=False):
        'set the axis data limits'
        if ignore:
            self.axes.dataLim.intervalx = vmin, vmax
        else:
            Vmin, Vmax = self.get_data_interval()
            self.axes.dataLim.intervalx = min(vmin, Vmin), max(vmax, Vmax)
        self.stale = True

    def set_default_intervals(self):
        'set the default limits for the axis interval if they are not mutated'
        xmin, xmax = 0., 1.
        dataMutated = self.axes.dataLim.mutatedx()
        viewMutated = self.axes.viewLim.mutatedx()
        if not dataMutated or not viewMutated:
            if self.converter is not None:
                info = self.converter.axisinfo(self.units, self)
                if info.default_limits is not None:
                    valmin, valmax = info.default_limits
                    xmin = self.converter.convert(valmin, self.units, self)
                    xmax = self.converter.convert(valmax, self.units, self)
            if not dataMutated:
                self.axes.dataLim.intervalx = xmin, xmax
            if not viewMutated:
                self.axes.viewLim.intervalx = xmin, xmax
        self.stale = True

    def get_tick_space(self):
        ends = self.axes.transAxes.transform([[0, 0], [1, 0]])
        length = ((ends[1][0] - ends[0][0]) / self.axes.figure.dpi) * 72.0
        tick = self._get_tick(True)
        # There is a heuristic here that the aspect ratio of tick text
        # is no more than 3:1
        size = tick.label1.get_size() * 3
        if size > 0:
            return int(np.floor(length / size))
        else:
            return 2**31 - 1


class YAxis(Axis):
    __name__ = 'yaxis'
    axis_name = 'y'

    def contains(self, mouseevent):
        """Test whether the mouse event occurred in the y axis.

        Returns *True* | *False*
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)

        x, y = mouseevent.x, mouseevent.y
        try:
            trans = self.axes.transAxes.inverted()
            xaxes, yaxes = trans.transform_point((x, y))
        except ValueError:
            return False, {}
        l, b = self.axes.transAxes.transform_point((0, 0))
        r, t = self.axes.transAxes.transform_point((1, 1))
        inaxis = yaxes >= 0 and yaxes <= 1 and (
            (x < l and x > l - self.pickradius) or
            (x > r and x < r + self.pickradius))
        return inaxis, {}

    def _get_tick(self, major):
        if major:
            tick_kw = self._major_tick_kw
        else:
            tick_kw = self._minor_tick_kw
        return YTick(self.axes, 0, '', major=major, **tick_kw)

    def _get_label(self):
        # x in display coords (updated by _update_label_position)
        # y in axes coords
        label = mtext.Text(x=0, y=0.5,
                           # todo: get the label position
                           fontproperties=font_manager.FontProperties(
                               size=rcParams['axes.labelsize'],
                               weight=rcParams['axes.labelweight']),
                           color=rcParams['axes.labelcolor'],
                           verticalalignment='bottom',
                           horizontalalignment='center',
                           rotation='vertical',
                           rotation_mode='anchor')
        label.set_transform(mtransforms.blended_transform_factory(
            mtransforms.IdentityTransform(), self.axes.transAxes))

        self._set_artist_props(label)
        self.label_position = 'left'
        return label

    def _get_offset_text(self):
        # x in display coords, y in axes coords (to be updated at draw time)
        offsetText = mtext.Text(x=0, y=0.5,
                                fontproperties=font_manager.FontProperties(
                                    size=rcParams['ytick.labelsize']
                                ),
                                color=rcParams['ytick.color'],
                                verticalalignment='baseline',
                                horizontalalignment='left')
        offsetText.set_transform(mtransforms.blended_transform_factory(
            self.axes.transAxes, mtransforms.IdentityTransform())
        )
        self._set_artist_props(offsetText)
        self.offset_text_position = 'left'
        return offsetText

    def _get_pixel_distance_along_axis(self, where, perturb):
        """
        Returns the amount, in data coordinates, that a single pixel
        corresponds to in the locality given by *where*, which is also given
        in data coordinates, and is a y coordinate.

        *perturb* is the amount to perturb the pixel.  Usually +0.5 or -0.5.

        Implementing this routine for an axis is optional; if present, it will
        ensure that no ticks are lost due to round-off at the extreme ends of
        an axis.
        """

        #
        # first figure out the pixel location of the "where" point.  We use
        # 1e-10 for the x point, so that we remain compatible with log axes.

        # transformation from data coords to display coords
        trans = self.axes.transData
        # transformation from display coords to data coords
        transinv = trans.inverted()
        pix = trans.transform_point((1e-10, where))
        # perturb the pixel
        ptp = transinv.transform_point((pix[0], pix[1] + perturb))
        dy = abs(ptp[1] - where)
        return dy

    def set_label_position(self, position):
        """
        Set the label position (left or right)

        ACCEPTS: [ 'left' | 'right' ]
        """
        self.label.set_rotation_mode('anchor')
        self.label.set_horizontalalignment('center')
        if position == 'left':
            self.label.set_verticalalignment('bottom')
        elif position == 'right':
            self.label.set_verticalalignment('top')
        else:
            raise ValueError("Position accepts only 'left' or 'right'")
        self.label_position = position
        self.stale = True

    def _get_tick_boxes_siblings(self, renderer):
        """
        Get the bounding boxes for this `.axis` and its siblings
        as set by `.Figure.align_xlabels` or  `.Figure.align_ylablels`.

        By default it just gets bboxes for self.
        """
        bboxes = []
        bboxes2 = []
        # get the Grouper that keeps track of y-label groups for this figure
        grp = self.figure._align_ylabel_grp
        # if we want to align labels from other axes:
        for axx in grp.get_siblings(self.axes):
            ticks_to_draw = axx.yaxis._update_ticks(renderer)
            tlb, tlb2 = axx.yaxis._get_tick_bboxes(ticks_to_draw, renderer)
            bboxes.extend(tlb)
            bboxes2.extend(tlb2)
        return bboxes, bboxes2

    def _update_label_position(self, renderer):
        """
        Update the label position based on the bounding box enclosing
        all the ticklabels and axis spine
        """
        if not self._autolabelpos:
            return

        # get bounding boxes for this axis and any siblings
        # that have been set by `fig.align_ylabels()`
        bboxes, bboxes2 = self._get_tick_boxes_siblings(renderer=renderer)

        x, y = self.label.get_position()
        if self.label_position == 'left':
            try:
                spine = self.axes.spines['left']
                spinebbox = spine.get_transform().transform_path(
                    spine.get_path()).get_extents()
            except KeyError:
                # use axes if spine doesn't exist
                spinebbox = self.axes.bbox
            bbox = mtransforms.Bbox.union(bboxes + [spinebbox])
            left = bbox.x0
            self.label.set_position(
                (left - self.labelpad * self.figure.dpi / 72.0, y)
            )

        else:
            try:
                spine = self.axes.spines['right']
                spinebbox = spine.get_transform().transform_path(
                    spine.get_path()).get_extents()
            except KeyError:
                # use axes if spine doesn't exist
                spinebbox = self.axes.bbox
            bbox = mtransforms.Bbox.union(bboxes2 + [spinebbox])
            right = bbox.x1

            self.label.set_position(
                (right + self.labelpad * self.figure.dpi / 72.0, y)
            )

    def _update_offset_text_position(self, bboxes, bboxes2):
        """
        Update the offset_text position based on the sequence of bounding
        boxes of all the ticklabels
        """
        x, y = self.offsetText.get_position()
        top = self.axes.bbox.ymax
        self.offsetText.set_position(
            (x, top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0)
        )

    def set_offset_position(self, position):
        """
        .. ACCEPTS: [ 'left' | 'right' ]
        """
        x, y = self.offsetText.get_position()
        if position == 'left':
            x = 0
        elif position == 'right':
            x = 1
        else:
            raise ValueError("Position accepts only [ 'left' | 'right' ]")

        self.offsetText.set_ha(position)
        self.offsetText.set_position((x, y))
        self.stale = True

    def get_text_widths(self, renderer):
        bbox, bbox2 = self.get_ticklabel_extents(renderer)
        # MGDTODO: Need a better way to get the pad
        padPixels = self.majorTicks[0].get_pad_pixels()

        left = 0.0
        if bbox.width:
            left += bbox.width + padPixels
        right = 0.0
        if bbox2.width:
            right += bbox2.width + padPixels

        if self.get_label_position() == 'left':
            left += self.label.get_window_extent(renderer).width + padPixels
        else:
            right += self.label.get_window_extent(renderer).width + padPixels
        return left, right

    def set_ticks_position(self, position):
        """
        Set the ticks position (left, right, both, default or none)
        'both' sets the ticks to appear on both positions, but does not
        change the tick labels.  'default' resets the tick positions to
        the default: ticks on both positions, labels at left.  'none'
        can be used if you don't want any ticks. 'none' and 'both'
        affect only the ticks, not the labels.

        ACCEPTS: [ 'left' | 'right' | 'both' | 'default' | 'none' ]
        """
        if position == 'right':
            self.set_tick_params(which='both', right=True, labelright=True,
                                 left=False, labelleft=False)
            self.set_offset_position(position)
        elif position == 'left':
            self.set_tick_params(which='both', right=False, labelright=False,
                                 left=True, labelleft=True)
            self.set_offset_position(position)
        elif position == 'both':
            self.set_tick_params(which='both', right=True,
                                 left=True)
        elif position == 'none':
            self.set_tick_params(which='both', right=False,
                                 left=False)
        elif position == 'default':
            self.set_tick_params(which='both', right=True, labelright=False,
                                 left=True, labelleft=True)
        else:
            raise ValueError("invalid position: %s" % position)
        self.stale = True

    def tick_right(self):
        """
        Move ticks and ticklabels (if present) to the right of the axes.
        """
        label = True
        if 'label1On' in self._major_tick_kw:
            label = (self._major_tick_kw['label1On']
                     or self._major_tick_kw['label2On'])
        self.set_ticks_position('right')
        # if labels were turned off before this was called
        # leave them off
        self.set_tick_params(which='both', labelright=label)

    def tick_left(self):
        """
        Move ticks and ticklabels (if present) to the left of the axes.
        """
        label = True
        if 'label1On' in self._major_tick_kw:
            label = (self._major_tick_kw['label1On']
                     or self._major_tick_kw['label2On'])
        self.set_ticks_position('left')
        # if labels were turned off before this was called
        # leave them off
        self.set_tick_params(which='both', labelleft=label)

    def get_ticks_position(self):
        """
        Return the ticks position (left, right, both or unknown)
        """
        majt = self.majorTicks[0]
        mT = self.minorTicks[0]

        majorRight = ((not majt.tick1On) and majt.tick2On and
                      (not majt.label1On) and majt.label2On)
        minorRight = ((not mT.tick1On) and mT.tick2On and
                      (not mT.label1On) and mT.label2On)
        if majorRight and minorRight:
            return 'right'

        majorLeft = (majt.tick1On and (not majt.tick2On) and
                     majt.label1On and (not majt.label2On))
        minorLeft = (mT.tick1On and (not mT.tick2On) and
                     mT.label1On and (not mT.label2On))
        if majorLeft and minorLeft:
            return 'left'

        majorDefault = (majt.tick1On and majt.tick2On and
                        majt.label1On and (not majt.label2On))
        minorDefault = (mT.tick1On and mT.tick2On and
                        mT.label1On and (not mT.label2On))
        if majorDefault and minorDefault:
            return 'default'

        return 'unknown'

    def get_view_interval(self):
        'return the Interval instance for this axis view limits'
        return self.axes.viewLim.intervaly

    def set_view_interval(self, vmin, vmax, ignore=False):
        """
        If *ignore* is *False*, the order of vmin, vmax
        does not matter; the original axis orientation will
        be preserved. In addition, the view limits can be
        expanded, but will not be reduced.  This method is
        for mpl internal use; for normal use, see
        :meth:`~matplotlib.axes.Axes.set_ylim`.

        """
        if ignore:
            self.axes.viewLim.intervaly = vmin, vmax
        else:
            Vmin, Vmax = self.get_view_interval()
            if Vmin < Vmax:
                self.axes.viewLim.intervaly = (min(vmin, vmax, Vmin),
                                               max(vmin, vmax, Vmax))
            else:
                self.axes.viewLim.intervaly = (max(vmin, vmax, Vmin),
                                               min(vmin, vmax, Vmax))
        self.stale = True

    def get_minpos(self):
        return self.axes.dataLim.minposy

    def get_data_interval(self):
        'return the Interval instance for this axis data limits'
        return self.axes.dataLim.intervaly

    def set_data_interval(self, vmin, vmax, ignore=False):
        'set the axis data limits'
        if ignore:
            self.axes.dataLim.intervaly = vmin, vmax
        else:
            Vmin, Vmax = self.get_data_interval()
            self.axes.dataLim.intervaly = min(vmin, Vmin), max(vmax, Vmax)
        self.stale = True

    def set_default_intervals(self):
        'set the default limits for the axis interval if they are not mutated'
        ymin, ymax = 0., 1.
        dataMutated = self.axes.dataLim.mutatedy()
        viewMutated = self.axes.viewLim.mutatedy()
        if not dataMutated or not viewMutated:
            if self.converter is not None:
                info = self.converter.axisinfo(self.units, self)
                if info.default_limits is not None:
                    valmin, valmax = info.default_limits
                    ymin = self.converter.convert(valmin, self.units, self)
                    ymax = self.converter.convert(valmax, self.units, self)
            if not dataMutated:
                self.axes.dataLim.intervaly = ymin, ymax
            if not viewMutated:
                self.axes.viewLim.intervaly = ymin, ymax
        self.stale = True

    def get_tick_space(self):
        ends = self.axes.transAxes.transform([[0, 0], [0, 1]])
        length = ((ends[1][1] - ends[0][1]) / self.axes.figure.dpi) * 72.0
        tick = self._get_tick(True)
        # Having a spacing of at least 2 just looks good.
        size = tick.label1.get_size() * 2.0
        if size > 0:
            return int(np.floor(length / size))
        else:
            return 2**31 - 1
