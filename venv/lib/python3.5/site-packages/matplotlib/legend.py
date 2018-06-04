"""
The legend module defines the Legend class, which is responsible for
drawing legends associated with axes and/or figures.

.. important::

    It is unlikely that you would ever create a Legend instance manually.
    Most users would normally create a legend via the
    :meth:`~matplotlib.axes.Axes.legend` function. For more details on legends
    there is also a :ref:`legend guide
    <sphx_glr_tutorials_intermediate_legend_guide.py>`.

The Legend class can be considered as a container of legend handles
and legend texts. Creation of corresponding legend handles from the
plot elements in the axes or figures (e.g., lines, patches, etc.) are
specified by the handler map, which defines the mapping between the
plot elements and the legend handlers to be used (the default legend
handlers are defined in the :mod:`~matplotlib.legend_handler` module).
Note that not all kinds of artist are supported by the legend yet by default
but it is possible to extend the legend handler's capabilities to support
arbitrary objects. See the :ref:`legend guide
<sphx_glr_tutorials_intermediate_legend_guide.py>` for more information.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import logging
import warnings

import numpy as np

from matplotlib import rcParams
from matplotlib import docstring
from matplotlib.artist import Artist, allow_rasterization
from matplotlib.cbook import silent_list, is_hashable
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Shadow, FancyBboxPatch
from matplotlib.collections import (LineCollection, RegularPolyCollection,
                                    CircleCollection, PathCollection,
                                    PolyCollection)
from matplotlib.transforms import Bbox, BboxBase, TransformedBbox
from matplotlib.transforms import BboxTransformTo, BboxTransformFrom

from matplotlib.offsetbox import HPacker, VPacker, TextArea, DrawingArea
from matplotlib.offsetbox import DraggableOffsetBox

from matplotlib.container import ErrorbarContainer, BarContainer, StemContainer
from . import legend_handler


class DraggableLegend(DraggableOffsetBox):
    def __init__(self, legend, use_blit=False, update="loc"):
        """
        Parameters
        ----------
        update : string
            If "loc", update *loc* parameter of legend upon finalizing.
            If "bbox", update *bbox_to_anchor* parameter.
        """
        self.legend = legend

        if update in ["loc", "bbox"]:
            self._update = update
        else:
            raise ValueError("update parameter '%s' is not supported." %
                             update)

        DraggableOffsetBox.__init__(self, legend, legend._legend_box,
                                    use_blit=use_blit)

    def artist_picker(self, legend, evt):
        return self.legend.contains(evt)

    def finalize_offset(self):
        loc_in_canvas = self.get_loc_in_canvas()

        if self._update == "loc":
            self._update_loc(loc_in_canvas)
        elif self._update == "bbox":
            self._update_bbox_to_anchor(loc_in_canvas)
        else:
            raise RuntimeError("update parameter '%s' is not supported." %
                               self.update)

    def _update_loc(self, loc_in_canvas):
        bbox = self.legend.get_bbox_to_anchor()

        # if bbox has zero width or height, the transformation is
        # ill-defined. Fall back to the defaul bbox_to_anchor.
        if bbox.width == 0 or bbox.height == 0:
            self.legend.set_bbox_to_anchor(None)
            bbox = self.legend.get_bbox_to_anchor()

        _bbox_transform = BboxTransformFrom(bbox)
        self.legend._loc = tuple(
            _bbox_transform.transform_point(loc_in_canvas)
        )

    def _update_bbox_to_anchor(self, loc_in_canvas):

        tr = self.legend.axes.transAxes
        loc_in_bbox = tr.transform_point(loc_in_canvas)

        self.legend.set_bbox_to_anchor(loc_in_bbox)


_legend_kw_doc = '''
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

'''
docstring.interpd.update(_legend_kw_doc=_legend_kw_doc)


class Legend(Artist):
    """
    Place a legend on the axes at location loc.

    """
    codes = {'best':         0,  # only implemented for axes legends
             'upper right':  1,
             'upper left':   2,
             'lower left':   3,
             'lower right':  4,
             'right':        5,
             'center left':  6,
             'center right': 7,
             'lower center': 8,
             'upper center': 9,
             'center':       10,
             }

    zorder = 5

    def __str__(self):
        return "Legend"

    @docstring.dedent_interpd
    def __init__(self, parent, handles, labels,
                 loc=None,
                 numpoints=None,    # the number of points in the legend line
                 markerscale=None,  # the relative size of legend markers
                                    # vs. original
                 markerfirst=True,  # controls ordering (left-to-right) of
                                    # legend marker and label
                 scatterpoints=None,    # number of scatter points
                 scatteryoffsets=None,
                 prop=None,          # properties for the legend texts
                 fontsize=None,        # keyword to set font size directly

                 # spacing & pad defined as a fraction of the font-size
                 borderpad=None,      # the whitespace inside the legend border
                 labelspacing=None,   # the vertical space between the legend
                                      # entries
                 handlelength=None,   # the length of the legend handles
                 handleheight=None,   # the height of the legend handles
                 handletextpad=None,  # the pad between the legend handle
                                      # and text
                 borderaxespad=None,  # the pad between the axes and legend
                                      # border
                 columnspacing=None,  # spacing between columns

                 ncol=1,     # number of columns
                 mode=None,  # mode for horizontal distribution of columns.
                             # None, "expand"

                 fancybox=None,  # True use a fancy box, false use a rounded
                                 # box, none use rc
                 shadow=None,
                 title=None,  # set a title for the legend

                 framealpha=None,  # set frame alpha
                 edgecolor=None,  # frame patch edgecolor
                 facecolor=None,  # frame patch facecolor

                 bbox_to_anchor=None,  # bbox that the legend will be anchored.
                 bbox_transform=None,  # transform for the bbox
                 frameon=None,  # draw frame
                 handler_map=None,
                 ):
        """

        Parameters
        ----------
        parent : `.Axes` or `.Figure`
            The artist that contains the legend.

        handles : sequence of `.Artist`
            A list of Artists (lines, patches) to be added to the legend.

        labels : sequence of strings
            A list of labels to show next to the artists. The length of handles
            and labels should be the same. If they are not, they are truncated
            to the smaller of both lengths.

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

        Notes
        -----

        Users can specify any arbitrary location for the legend using the
        *bbox_to_anchor* keyword argument. bbox_to_anchor can be an instance
        of BboxBase(or its derivatives) or a tuple of 2 or 4 floats.
        See :meth:`set_bbox_to_anchor` for more detail.

        The legend location can be specified by setting *loc* with a tuple of
        2 floats, which is interpreted as the lower-left corner of the legend
        in the normalized axes coordinate.
        """
        # local import only to avoid circularity
        from matplotlib.axes import Axes
        from matplotlib.figure import Figure

        Artist.__init__(self)

        if prop is None:
            if fontsize is not None:
                self.prop = FontProperties(size=fontsize)
            else:
                self.prop = FontProperties(size=rcParams["legend.fontsize"])
        elif isinstance(prop, dict):
            self.prop = FontProperties(**prop)
            if "size" not in prop:
                self.prop.set_size(rcParams["legend.fontsize"])
        else:
            self.prop = prop

        self._fontsize = self.prop.get_size_in_points()

        self.texts = []
        self.legendHandles = []
        self._legend_title_box = None

        #: A dictionary with the extra handler mappings for this Legend
        #: instance.
        self._custom_handler_map = handler_map

        locals_view = locals()
        for name in ["numpoints", "markerscale", "shadow", "columnspacing",
                     "scatterpoints", "handleheight", 'borderpad',
                     'labelspacing', 'handlelength', 'handletextpad',
                     'borderaxespad']:
            if locals_view[name] is None:
                value = rcParams["legend." + name]
            else:
                value = locals_view[name]
            setattr(self, name, value)
        del locals_view
        # trim handles and labels if illegal label...
        _lab, _hand = [], []
        for label, handle in zip(labels, handles):
            if (isinstance(label, six.string_types) and
                    label.startswith('_')):
                warnings.warn('The handle {!r} has a label of {!r} which '
                              'cannot be automatically added to the '
                              'legend.'.format(handle, label))
            else:
                _lab.append(label)
                _hand.append(handle)
        labels, handles = _lab, _hand

        handles = list(handles)
        if len(handles) < 2:
            ncol = 1
        self._ncol = ncol

        if self.numpoints <= 0:
            raise ValueError("numpoints must be > 0; it was %d" % numpoints)

        # introduce y-offset for handles of the scatter plot
        if scatteryoffsets is None:
            self._scatteryoffsets = np.array([3. / 8., 4. / 8., 2.5 / 8.])
        else:
            self._scatteryoffsets = np.asarray(scatteryoffsets)
        reps = self.scatterpoints // len(self._scatteryoffsets) + 1
        self._scatteryoffsets = np.tile(self._scatteryoffsets,
                                        reps)[:self.scatterpoints]

        # _legend_box is an OffsetBox instance that contains all
        # legend items and will be initialized from _init_legend_box()
        # method.
        self._legend_box = None

        if isinstance(parent, Axes):
            self.isaxes = True
            self.axes = parent
            self.set_figure(parent.figure)
        elif isinstance(parent, Figure):
            self.isaxes = False
            self.set_figure(parent)
        else:
            raise TypeError("Legend needs either Axes or Figure as parent")
        self.parent = parent

        if loc is None:
            loc = rcParams["legend.loc"]
            if not self.isaxes and loc in [0, 'best']:
                loc = 'upper right'
        if isinstance(loc, six.string_types):
            if loc not in self.codes:
                if self.isaxes:
                    warnings.warn('Unrecognized location "%s". Falling back '
                                  'on "best"; valid locations are\n\t%s\n'
                                  % (loc, '\n\t'.join(self.codes)))
                    loc = 0
                else:
                    warnings.warn('Unrecognized location "%s". Falling back '
                                  'on "upper right"; '
                                  'valid locations are\n\t%s\n'
                                  % (loc, '\n\t'.join(self.codes)))
                    loc = 1
            else:
                loc = self.codes[loc]
        if not self.isaxes and loc == 0:
            warnings.warn('Automatic legend placement (loc="best") not '
                          'implemented for figure legend. '
                          'Falling back on "upper right".')
            loc = 1

        self._mode = mode
        self.set_bbox_to_anchor(bbox_to_anchor, bbox_transform)

        # We use FancyBboxPatch to draw a legend frame. The location
        # and size of the box will be updated during the drawing time.

        if facecolor is None:
            facecolor = rcParams["legend.facecolor"]
        if facecolor == 'inherit':
            facecolor = rcParams["axes.facecolor"]

        if edgecolor is None:
            edgecolor = rcParams["legend.edgecolor"]
        if edgecolor == 'inherit':
            edgecolor = rcParams["axes.edgecolor"]

        self.legendPatch = FancyBboxPatch(
            xy=(0.0, 0.0), width=1., height=1.,
            facecolor=facecolor,
            edgecolor=edgecolor,
            mutation_scale=self._fontsize,
            snap=True
            )

        # The width and height of the legendPatch will be set (in the
        # draw()) to the length that includes the padding. Thus we set
        # pad=0 here.
        if fancybox is None:
            fancybox = rcParams["legend.fancybox"]

        if fancybox:
            self.legendPatch.set_boxstyle("round", pad=0,
                                          rounding_size=0.2)
        else:
            self.legendPatch.set_boxstyle("square", pad=0)

        self._set_artist_props(self.legendPatch)

        self._drawFrame = frameon
        if frameon is None:
            self._drawFrame = rcParams["legend.frameon"]

        # init with null renderer
        self._init_legend_box(handles, labels, markerfirst)

        # If shadow is activated use framealpha if not
        # explicitly passed. See Issue 8943
        if framealpha is None:
            if shadow:
                self.get_frame().set_alpha(1)
            else:
                self.get_frame().set_alpha(rcParams["legend.framealpha"])
        else:
            self.get_frame().set_alpha(framealpha)

        self._loc = loc
        self.set_title(title)
        self._last_fontsize_points = self._fontsize
        self._draggable = None

    def _set_artist_props(self, a):
        """
        Set the boilerplate props for artists added to axes.
        """
        a.set_figure(self.figure)
        if self.isaxes:
            # a.set_axes(self.axes)
            a.axes = self.axes

        a.set_transform(self.get_transform())

    def _set_loc(self, loc):
        # find_offset function will be provided to _legend_box and
        # _legend_box will draw itself at the location of the return
        # value of the find_offset.
        self._loc_real = loc
        self.stale = True
        self._legend_box.set_offset(self._findoffset)

    def _get_loc(self):
        return self._loc_real

    _loc = property(_get_loc, _set_loc)

    def _findoffset(self, width, height, xdescent, ydescent, renderer):
        "Helper function to locate the legend."

        if self._loc == 0:  # "best".
            x, y = self._find_best_position(width, height, renderer)
        elif self._loc in Legend.codes.values():  # Fixed location.
            bbox = Bbox.from_bounds(0, 0, width, height)
            x, y = self._get_anchored_bbox(self._loc, bbox,
                                           self.get_bbox_to_anchor(),
                                           renderer)
        else:  # Axes or figure coordinates.
            fx, fy = self._loc
            bbox = self.get_bbox_to_anchor()
            x, y = bbox.x0 + bbox.width * fx, bbox.y0 + bbox.height * fy

        return x + xdescent, y + ydescent

    @allow_rasterization
    def draw(self, renderer):
        "Draw everything that belongs to the legend."
        if not self.get_visible():
            return

        renderer.open_group('legend')

        fontsize = renderer.points_to_pixels(self._fontsize)

        # if mode == fill, set the width of the legend_box to the
        # width of the paret (minus pads)
        if self._mode in ["expand"]:
            pad = 2 * (self.borderaxespad + self.borderpad) * fontsize
            self._legend_box.set_width(self.get_bbox_to_anchor().width - pad)

        # update the location and size of the legend. This needs to
        # be done in any case to clip the figure right.
        bbox = self._legend_box.get_window_extent(renderer)
        self.legendPatch.set_bounds(bbox.x0, bbox.y0,
                                    bbox.width, bbox.height)
        self.legendPatch.set_mutation_scale(fontsize)

        if self._drawFrame:
            if self.shadow:
                shadow = Shadow(self.legendPatch, 2, -2)
                shadow.draw(renderer)

            self.legendPatch.draw(renderer)

        self._legend_box.draw(renderer)

        renderer.close_group('legend')
        self.stale = False

    def _approx_text_height(self, renderer=None):
        """
        Return the approximate height of the text. This is used to place
        the legend handle.
        """
        if renderer is None:
            return self._fontsize
        else:
            return renderer.points_to_pixels(self._fontsize)

    # _default_handler_map defines the default mapping between plot
    # elements and the legend handlers.

    _default_handler_map = {
        StemContainer: legend_handler.HandlerStem(),
        ErrorbarContainer: legend_handler.HandlerErrorbar(),
        Line2D: legend_handler.HandlerLine2D(),
        Patch: legend_handler.HandlerPatch(),
        LineCollection: legend_handler.HandlerLineCollection(),
        RegularPolyCollection: legend_handler.HandlerRegularPolyCollection(),
        CircleCollection: legend_handler.HandlerCircleCollection(),
        BarContainer: legend_handler.HandlerPatch(
            update_func=legend_handler.update_from_first_child),
        tuple: legend_handler.HandlerTuple(),
        PathCollection: legend_handler.HandlerPathCollection(),
        PolyCollection: legend_handler.HandlerPolyCollection()
        }

    # (get|set|update)_default_handler_maps are public interfaces to
    # modify the default handler map.

    @classmethod
    def get_default_handler_map(cls):
        """
        A class method that returns the default handler map.
        """
        return cls._default_handler_map

    @classmethod
    def set_default_handler_map(cls, handler_map):
        """
        A class method to set the default handler map.
        """
        cls._default_handler_map = handler_map

    @classmethod
    def update_default_handler_map(cls, handler_map):
        """
        A class method to update the default handler map.
        """
        cls._default_handler_map.update(handler_map)

    def get_legend_handler_map(self):
        """
        Return the handler map.
        """

        default_handler_map = self.get_default_handler_map()

        if self._custom_handler_map:
            hm = default_handler_map.copy()
            hm.update(self._custom_handler_map)
            return hm
        else:
            return default_handler_map

    @staticmethod
    def get_legend_handler(legend_handler_map, orig_handle):
        """
        Return a legend handler from *legend_handler_map* that
        corresponds to *orig_handler*.

        *legend_handler_map* should be a dictionary object (that is
        returned by the get_legend_handler_map method).

        It first checks if the *orig_handle* itself is a key in the
        *legend_hanler_map* and return the associated value.
        Otherwise, it checks for each of the classes in its
        method-resolution-order. If no matching key is found, it
        returns ``None``.
        """
        if is_hashable(orig_handle):
            try:
                return legend_handler_map[orig_handle]
            except KeyError:
                pass

        for handle_type in type(orig_handle).mro():
            try:
                return legend_handler_map[handle_type]
            except KeyError:
                pass

        return None

    def _init_legend_box(self, handles, labels, markerfirst=True):
        """
        Initialize the legend_box. The legend_box is an instance of
        the OffsetBox, which is packed with legend handles and
        texts. Once packed, their location is calculated during the
        drawing time.
        """

        fontsize = self._fontsize

        # legend_box is a HPacker, horizontally packed with
        # columns. Each column is a VPacker, vertically packed with
        # legend items. Each legend item is HPacker packed with
        # legend handleBox and labelBox. handleBox is an instance of
        # offsetbox.DrawingArea which contains legend handle. labelBox
        # is an instance of offsetbox.TextArea which contains legend
        # text.

        text_list = []  # the list of text instances
        handle_list = []  # the list of text instances
        handles_and_labels = []

        label_prop = dict(verticalalignment='baseline',
                          horizontalalignment='left',
                          fontproperties=self.prop,
                          )

        # The approximate height and descent of text. These values are
        # only used for plotting the legend handle.
        descent = 0.35 * self._approx_text_height() * (self.handleheight - 0.7)
        # 0.35 and 0.7 are just heuristic numbers and may need to be improved.
        height = self._approx_text_height() * self.handleheight - descent
        # each handle needs to be drawn inside a box of (x, y, w, h) =
        # (0, -descent, width, height).  And their coordinates should
        # be given in the display coordinates.

        # The transformation of each handle will be automatically set
        # to self.get_trasnform(). If the artist does not use its
        # default transform (e.g., Collections), you need to
        # manually set their transform to the self.get_transform().
        legend_handler_map = self.get_legend_handler_map()

        for orig_handle, lab in zip(handles, labels):
            handler = self.get_legend_handler(legend_handler_map, orig_handle)
            if handler is None:
                warnings.warn(
                    "Legend does not support {!r} instances.\nA proxy artist "
                    "may be used instead.\nSee: "
                    "http://matplotlib.org/users/legend_guide.html"
                    "#creating-artists-specifically-for-adding-to-the-legend-"
                    "aka-proxy-artists".format(orig_handle)
                )
                # We don't have a handle for this artist, so we just defer
                # to None.
                handle_list.append(None)
            else:
                textbox = TextArea(lab, textprops=label_prop,
                                   multilinebaseline=True,
                                   minimumdescent=True)
                handlebox = DrawingArea(width=self.handlelength * fontsize,
                                        height=height,
                                        xdescent=0., ydescent=descent)

                text_list.append(textbox._text)
                # Create the artist for the legend which represents the
                # original artist/handle.
                handle_list.append(handler.legend_artist(self, orig_handle,
                                                         fontsize, handlebox))
                handles_and_labels.append((handlebox, textbox))

        if handles_and_labels:
            # We calculate number of rows in each column. The first
            # (num_largecol) columns will have (nrows+1) rows, and remaining
            # (num_smallcol) columns will have (nrows) rows.
            ncol = min(self._ncol, len(handles_and_labels))
            nrows, num_largecol = divmod(len(handles_and_labels), ncol)
            num_smallcol = ncol - num_largecol
            # starting index of each column and number of rows in it.
            rows_per_col = [nrows + 1] * num_largecol + [nrows] * num_smallcol
            start_idxs = np.concatenate([[0], np.cumsum(rows_per_col)[:-1]])
            cols = zip(start_idxs, rows_per_col)
        else:
            cols = []

        columnbox = []
        for i0, di in cols:
            # pack handleBox and labelBox into itemBox
            itemBoxes = [HPacker(pad=0,
                                 sep=self.handletextpad * fontsize,
                                 children=[h, t] if markerfirst else [t, h],
                                 align="baseline")
                         for h, t in handles_and_labels[i0:i0 + di]]
            # minimumdescent=False for the text of the last row of the column
            if markerfirst:
                itemBoxes[-1].get_children()[1].set_minimumdescent(False)
            else:
                itemBoxes[-1].get_children()[0].set_minimumdescent(False)

            # pack columnBox
            alignment = "baseline" if markerfirst else "right"
            columnbox.append(VPacker(pad=0,
                                     sep=self.labelspacing * fontsize,
                                     align=alignment,
                                     children=itemBoxes))

        mode = "expand" if self._mode == "expand" else "fixed"
        sep = self.columnspacing * fontsize
        self._legend_handle_box = HPacker(pad=0,
                                          sep=sep, align="baseline",
                                          mode=mode,
                                          children=columnbox)
        self._legend_title_box = TextArea("")
        self._legend_box = VPacker(pad=self.borderpad * fontsize,
                                   sep=self.labelspacing * fontsize,
                                   align="center",
                                   children=[self._legend_title_box,
                                             self._legend_handle_box])
        self._legend_box.set_figure(self.figure)
        self.texts = text_list
        self.legendHandles = handle_list

    def _auto_legend_data(self):
        """
        Returns list of vertices and extents covered by the plot.

        Returns a two long list.

        First element is a list of (x, y) vertices (in
        display-coordinates) covered by all the lines and line
        collections, in the legend's handles.

        Second element is a list of bounding boxes for all the patches in
        the legend's handles.
        """
        # should always hold because function is only called internally
        assert self.isaxes

        ax = self.parent
        bboxes = []
        lines = []
        offsets = []

        for handle in ax.lines:
            assert isinstance(handle, Line2D)
            path = handle.get_path()
            trans = handle.get_transform()
            tpath = trans.transform_path(path)
            lines.append(tpath)

        for handle in ax.patches:
            assert isinstance(handle, Patch)

            if isinstance(handle, Rectangle):
                transform = handle.get_data_transform()
                bboxes.append(handle.get_bbox().transformed(transform))
            else:
                transform = handle.get_transform()
                bboxes.append(handle.get_path().get_extents(transform))

        for handle in ax.collections:
            transform, transOffset, hoffsets, paths = handle._prepare_points()

            if len(hoffsets):
                for offset in transOffset.transform(hoffsets):
                    offsets.append(offset)

        try:
            vertices = np.concatenate([l.vertices for l in lines])
        except ValueError:
            vertices = np.array([])

        return [vertices, bboxes, lines, offsets]

    def draw_frame(self, b):
        '''
        Set draw frame to b.

        Parameters
        ----------
        b : bool
        '''
        self.set_frame_on(b)

    def get_children(self):
        'Return a list of child artists.'
        children = []
        if self._legend_box:
            children.append(self._legend_box)
        children.append(self.get_frame())

        return children

    def get_frame(self):
        '''
        Return the `~.patches.Rectangle` instances used to frame the legend.
        '''
        return self.legendPatch

    def get_lines(self):
        'Return a list of `~.lines.Line2D` instances in the legend.'
        return [h for h in self.legendHandles if isinstance(h, Line2D)]

    def get_patches(self):
        'Return a list of `~.patches.Patch` instances in the legend.'
        return silent_list('Patch',
                           [h for h in self.legendHandles
                            if isinstance(h, Patch)])

    def get_texts(self):
        'Return a list of `~.text.Text` instances in the legend.'
        return silent_list('Text', self.texts)

    def set_title(self, title, prop=None):
        """
        Set the legend title. Fontproperties can be optionally set
        with *prop* parameter.
        """
        self._legend_title_box._text.set_text(title)

        if prop is not None:
            if isinstance(prop, dict):
                prop = FontProperties(**prop)
            self._legend_title_box._text.set_fontproperties(prop)

        if title:
            self._legend_title_box.set_visible(True)
        else:
            self._legend_title_box.set_visible(False)
        self.stale = True

    def get_title(self):
        'Return the `.Text` instance for the legend title.'
        return self._legend_title_box._text

    def get_window_extent(self, *args, **kwargs):
        'Return extent of the legend.'
        return self.legendPatch.get_window_extent(*args, **kwargs)

    def get_frame_on(self):
        """Get whether the legend box patch is drawn."""
        return self._drawFrame

    def set_frame_on(self, b):
        """
        Set whether the legend box patch is drawn.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._drawFrame = b
        self.stale = True

    def get_bbox_to_anchor(self):
        """Return the bbox that the legend will be anchored to."""
        if self._bbox_to_anchor is None:
            return self.parent.bbox
        else:
            return self._bbox_to_anchor

    def set_bbox_to_anchor(self, bbox, transform=None):
        """
        Set the bbox that the legend will be anchored to.

        *bbox* can be

        - A `.BboxBase` instance
        - A tuple of ``(left, bottom, width, height)`` in the given transform
          (normalized axes coordinate if None)
        - A tuple of ``(left, bottom)`` where the width and height will be
          assumed to be zero.
        """
        if bbox is None:
            self._bbox_to_anchor = None
            return
        elif isinstance(bbox, BboxBase):
            self._bbox_to_anchor = bbox
        else:
            try:
                l = len(bbox)
            except TypeError:
                raise ValueError("Invalid argument for bbox : %s" % str(bbox))

            if l == 2:
                bbox = [bbox[0], bbox[1], 0, 0]

            self._bbox_to_anchor = Bbox.from_bounds(*bbox)

        if transform is None:
            transform = BboxTransformTo(self.parent.bbox)

        self._bbox_to_anchor = TransformedBbox(self._bbox_to_anchor,
                                               transform)
        self.stale = True

    def _get_anchored_bbox(self, loc, bbox, parentbbox, renderer):
        """
        Place the *bbox* inside the *parentbbox* according to a given
        location code. Return the (x,y) coordinate of the bbox.

        - loc: a location code in range(1, 11).
          This corresponds to the possible values for self._loc, excluding
          "best".

        - bbox: bbox to be placed, display coordinate units.
        - parentbbox: a parent box which will contain the bbox. In
            display coordinates.
        """
        assert loc in range(1, 11)  # called only internally

        BEST, UR, UL, LL, LR, R, CL, CR, LC, UC, C = range(11)

        anchor_coefs = {UR: "NE",
                        UL: "NW",
                        LL: "SW",
                        LR: "SE",
                        R: "E",
                        CL: "W",
                        CR: "E",
                        LC: "S",
                        UC: "N",
                        C: "C"}

        c = anchor_coefs[loc]

        fontsize = renderer.points_to_pixels(self._fontsize)
        container = parentbbox.padded(-(self.borderaxespad) * fontsize)
        anchored_box = bbox.anchored(c, container=container)
        return anchored_box.x0, anchored_box.y0

    def _find_best_position(self, width, height, renderer, consider=None):
        """
        Determine the best location to place the legend.

        *consider* is a list of ``(x, y)`` pairs to consider as a potential
        lower-left corner of the legend. All are display coords.
        """
        # should always hold because function is only called internally
        assert self.isaxes

        verts, bboxes, lines, offsets = self._auto_legend_data()

        bbox = Bbox.from_bounds(0, 0, width, height)
        if consider is None:
            consider = [self._get_anchored_bbox(x, bbox,
                                                self.get_bbox_to_anchor(),
                                                renderer)
                        for x in range(1, len(self.codes))]

        candidates = []
        for idx, (l, b) in enumerate(consider):
            legendBox = Bbox.from_bounds(l, b, width, height)
            badness = 0
            # XXX TODO: If markers are present, it would be good to
            # take them into account when checking vertex overlaps in
            # the next line.
            badness = (legendBox.count_contains(verts)
                       + legendBox.count_contains(offsets)
                       + legendBox.count_overlaps(bboxes)
                       + sum(line.intersects_bbox(legendBox, filled=False)
                             for line in lines))
            if badness == 0:
                return l, b
            # Include the index to favor lower codes in case of a tie.
            candidates.append((badness, idx, (l, b)))

        _, _, (l, b) = min(candidates)
        return l, b

    def contains(self, event):
        return self.legendPatch.contains(event)

    def draggable(self, state=None, use_blit=False, update="loc"):
        """
        Set the draggable state -- if state is

          * None : toggle the current state

          * True : turn draggable on

          * False : turn draggable off

        If draggable is on, you can drag the legend on the canvas with
        the mouse. The `.DraggableLegend` helper instance is returned if
        draggable is on.

        The update parameter control which parameter of the legend changes
        when dragged. If update is "loc", the *loc* parameter of the legend
        is changed. If "bbox", the *bbox_to_anchor* parameter is changed.
        """
        is_draggable = self._draggable is not None

        # if state is None we'll toggle
        if state is None:
            state = not is_draggable

        if state:
            if self._draggable is None:
                self._draggable = DraggableLegend(self,
                                                  use_blit,
                                                  update=update)
        else:
            if self._draggable is not None:
                self._draggable.disconnect()
            self._draggable = None

        return self._draggable


# Helper functions to parse legend arguments for both `figure.legend` and
# `axes.legend`:
def _get_legend_handles(axs, legend_handler_map=None):
    """
    Return a generator of artists that can be used as handles in
    a legend.

    """
    handles_original = []
    for ax in axs:
        handles_original += (ax.lines + ax.patches +
                             ax.collections + ax.containers)
        # support parasite axes:
        if hasattr(ax, 'parasites'):
            for axx in ax.parasites:
                handles_original += (axx.lines + axx.patches +
                                     axx.collections + axx.containers)

    handler_map = Legend.get_default_handler_map()

    if legend_handler_map is not None:
        handler_map = handler_map.copy()
        handler_map.update(legend_handler_map)

    has_handler = Legend.get_legend_handler

    for handle in handles_original:
        label = handle.get_label()
        if label != '_nolegend_' and has_handler(handler_map, handle):
            yield handle


def _get_legend_handles_labels(axs, legend_handler_map=None):
    """
    Return handles and labels for legend, internal method.

    """
    handles = []
    labels = []

    for handle in _get_legend_handles(axs, legend_handler_map):
        label = handle.get_label()
        if (label and not label.startswith('_')):
            handles.append(handle)
            labels.append(label)
    return handles, labels


def _parse_legend_args(axs, *args, **kwargs):
    """
    Get the handles and labels from the calls to either ``figure.legend``
    or ``axes.legend``.

    ``axs`` is a list of axes (to get legend artists from)
    """
    log = logging.getLogger(__name__)

    handlers = kwargs.get('handler_map', {}) or {}

    # Support handles and labels being passed as keywords.
    handles = kwargs.pop('handles', None)
    labels = kwargs.pop('labels', None)

    extra_args = ()

    if (handles is not None or labels is not None) and len(args):
        warnings.warn("You have mixed positional and keyword "
                      "arguments, some input may be "
                      "discarded.")

    # if got both handles and labels as kwargs, make same length
    if handles and labels:
        handles, labels = zip(*zip(handles, labels))

    elif handles is not None and labels is None:
        labels = [handle.get_label() for handle in handles]

    elif labels is not None and handles is None:
        # Get as many handles as there are labels.
        handles = [handle for handle, label
                   in zip(_get_legend_handles(axs, handlers), labels)]

    # No arguments - automatically detect labels and handles.
    elif len(args) == 0:
        handles, labels = _get_legend_handles_labels(axs, handlers)
        if not handles:
            log.warning('No handles with labels found to put in legend.')

    # One argument. User defined labels - automatic handle detection.
    elif len(args) == 1:
        labels, = args
        # Get as many handles as there are labels.
        handles = [handle for handle, label
                   in zip(_get_legend_handles(axs, handlers), labels)]

    # Two arguments:
    #   * user defined handles and labels
    elif len(args) >= 2:
        handles, labels = args[:2]
        extra_args = args[2:]

    else:
        raise TypeError('Invalid arguments to legend.')

    return handles, labels, extra_args, kwargs
