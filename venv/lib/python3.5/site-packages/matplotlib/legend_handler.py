"""
This module defines default legend handlers.

It is strongly encouraged to have read the :ref:`legend guide
<sphx_glr_tutorials_intermediate_legend_guide.py>` before this documentation.

Legend handlers are expected to be a callable object with a following
signature. ::

    legend_handler(legend, orig_handle, fontsize, handlebox)

Where *legend* is the legend itself, *orig_handle* is the original
plot, *fontsize* is the fontsize in pixels, and *handlebox* is a
OffsetBox instance. Within the call, you should create relevant
artists (using relevant properties from the *legend* and/or
*orig_handle*) and add them into the handlebox. The artists needs to
be scaled according to the fontsize (note that the size is in pixel,
i.e., this is dpi-scaled value).

This module includes definition of several legend handler classes
derived from the base class (HandlerBase) with the following method::

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip
from itertools import cycle

import numpy as np

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors


def update_from_first_child(tgt, src):
    tgt.update_from(src.get_children()[0])


class HandlerBase(object):
    """
    A Base class for default legend handlers.

    The derived classes are meant to override *create_artists* method, which
    has a following signature.::

      def create_artists(self, legend, orig_handle,
                         xdescent, ydescent, width, height, fontsize,
                         trans):

    The overridden method needs to create artists of the given
    transform that fits in the given dimension (xdescent, ydescent,
    width, height) that are scaled by fontsize if necessary.

    """
    def __init__(self, xpad=0., ypad=0., update_func=None):
        self._xpad, self._ypad = xpad, ypad
        self._update_prop_func = update_func

    def _update_prop(self, legend_handle, orig_handle):
        if self._update_prop_func is None:
            self._default_update_prop(legend_handle, orig_handle)
        else:
            self._update_prop_func(legend_handle, orig_handle)

    def _default_update_prop(self, legend_handle, orig_handle):
        legend_handle.update_from(orig_handle)

    def update_prop(self, legend_handle, orig_handle, legend):

        self._update_prop(legend_handle, orig_handle)

        legend._set_artist_props(legend_handle)
        legend_handle.set_clip_box(None)
        legend_handle.set_clip_path(None)

    def adjust_drawing_area(self, legend, orig_handle,
                            xdescent, ydescent, width, height, fontsize,
                            ):
        xdescent = xdescent - self._xpad * fontsize
        ydescent = ydescent - self._ypad * fontsize
        width = width - self._xpad * fontsize
        height = height - self._ypad * fontsize
        return xdescent, ydescent, width, height

    def legend_artist(self, legend, orig_handle,
                      fontsize, handlebox):
        """
        Return the artist that this HandlerBase generates for the given
        original artist/handle.

        Parameters
        ----------
        legend : :class:`matplotlib.legend.Legend` instance
            The legend for which these legend artists are being created.
        orig_handle : :class:`matplotlib.artist.Artist` or similar
            The object for which these legend artists are being created.
        fontsize : float or int
            The fontsize in pixels. The artists being created should
            be scaled according to the given fontsize.
        handlebox : :class:`matplotlib.offsetbox.OffsetBox` instance
            The box which has been created to hold this legend entry's
            artists. Artists created in the `legend_artist` method must
            be added to this handlebox inside this method.

        """
        xdescent, ydescent, width, height = self.adjust_drawing_area(
                 legend, orig_handle,
                 handlebox.xdescent, handlebox.ydescent,
                 handlebox.width, handlebox.height,
                 fontsize)
        artists = self.create_artists(legend, orig_handle,
                                      xdescent, ydescent, width, height,
                                      fontsize, handlebox.get_transform())

        # create_artists will return a list of artists.
        for a in artists:
            handlebox.add_artist(a)

        # we only return the first artist
        return artists[0]

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        raise NotImplementedError('Derived must override')


class HandlerNpoints(HandlerBase):
    """
    A legend handler that shows *numpoints* points in the legend entry.
    """
    def __init__(self, marker_pad=0.3, numpoints=None, **kw):
        """
        Parameters
        ----------
        marker_pad : float
            Padding between points in legend entry.

        numpoints : int
            Number of points to show in legend entry.

        Notes
        -----
        Any other keyword arguments are given to `HandlerBase`.
        """
        HandlerBase.__init__(self, **kw)

        self._numpoints = numpoints
        self._marker_pad = marker_pad

    def get_numpoints(self, legend):
        if self._numpoints is None:
            return legend.numpoints
        else:
            return self._numpoints

    def get_xdata(self, legend, xdescent, ydescent, width, height, fontsize):
        numpoints = self.get_numpoints(legend)
        if numpoints > 1:
            # we put some pad here to compensate the size of the marker
            pad = self._marker_pad * fontsize
            xdata = np.linspace(-xdescent + pad,
                                -xdescent + width - pad,
                                numpoints)
            xdata_marker = xdata
        else:
            xdata = np.linspace(-xdescent, -xdescent + width, 2)
            xdata_marker = [-xdescent + 0.5 * width]
        return xdata, xdata_marker


class HandlerNpointsYoffsets(HandlerNpoints):
    """
    A legend handler that shows *numpoints* in the legend, and allows them to
    be individually offest in the y-direction.
    """
    def __init__(self, numpoints=None, yoffsets=None, **kw):
        """
        Parameters
        ----------
        numpoints : int
            Number of points to show in legend entry.

        yoffsets : array of floats
            Length *numpoints* list of y offsets for each point in
            legend entry.

        Notes
        -----
        Any other keyword arguments are given to `HandlerNpoints`.
        """
        HandlerNpoints.__init__(self, numpoints=numpoints, **kw)
        self._yoffsets = yoffsets

    def get_ydata(self, legend, xdescent, ydescent, width, height, fontsize):
        if self._yoffsets is None:
            ydata = height * legend._scatteryoffsets
        else:
            ydata = height * np.asarray(self._yoffsets)

        return ydata


class HandlerLine2D(HandlerNpoints):
    """
    Handler for `.Line2D` instances.
    """
    def __init__(self, marker_pad=0.3, numpoints=None, **kw):
        """
        Parameters
        ----------
        marker_pad : float
            Padding between points in legend entry.

        numpoints : int
            Number of points to show in legend entry.

        Notes
        -----
        Any other keyword arguments are given to `HandlerNpoints`.
        """
        HandlerNpoints.__init__(self, marker_pad=marker_pad,
                                numpoints=numpoints, **kw)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)

        ydata = ((height - ydescent) / 2.) * np.ones(xdata.shape, float)
        legline = Line2D(xdata, ydata)

        self.update_prop(legline, orig_handle, legend)
        legline.set_drawstyle('default')
        legline.set_marker("")

        legline_marker = Line2D(xdata_marker, ydata[:len(xdata_marker)])
        self.update_prop(legline_marker, orig_handle, legend)
        legline_marker.set_linestyle('None')
        if legend.markerscale != 1:
            newsz = legline_marker.get_markersize() * legend.markerscale
            legline_marker.set_markersize(newsz)
        # we don't want to add this to the return list because
        # the texts and handles are assumed to be in one-to-one
        # correspondence.
        legline._legmarker = legline_marker

        legline.set_transform(trans)
        legline_marker.set_transform(trans)

        return [legline, legline_marker]


class HandlerPatch(HandlerBase):
    """
    Handler for `.Patch` instances.
    """
    def __init__(self, patch_func=None, **kw):
        """
        Parameters
        ----------
        patch_func : callable, optional
            The function that creates the legend key artist.
            *patch_func* should have the signature::

                def patch_func(legend=legend, orig_handle=orig_handle,
                               xdescent=xdescent, ydescent=ydescent,
                               width=width, height=height, fontsize=fontsize)

            Subsequently the created artist will have its ``update_prop`` method
            called and the appropriate transform will be applied.

        Notes
        -----
        Any other keyword arguments are given to `HandlerBase`.
        """
        HandlerBase.__init__(self, **kw)
        self._patch_func = patch_func

    def _create_patch(self, legend, orig_handle,
                      xdescent, ydescent, width, height, fontsize):
        if self._patch_func is None:
            p = Rectangle(xy=(-xdescent, -ydescent),
                          width=width, height=height)
        else:
            p = self._patch_func(legend=legend, orig_handle=orig_handle,
                                 xdescent=xdescent, ydescent=ydescent,
                                 width=width, height=height, fontsize=fontsize)
        return p

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        p = self._create_patch(legend, orig_handle,
                               xdescent, ydescent, width, height, fontsize)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


class HandlerLineCollection(HandlerLine2D):
    """
    Handler for `.LineCollection` instances.
    """
    def get_numpoints(self, legend):
        if self._numpoints is None:
            return legend.scatterpoints
        else:
            return self._numpoints

    def _default_update_prop(self, legend_handle, orig_handle):
        lw = orig_handle.get_linewidths()[0]
        dashes = orig_handle._us_linestyles[0]
        color = orig_handle.get_colors()[0]
        legend_handle.set_color(color)
        legend_handle.set_linestyle(dashes)
        legend_handle.set_linewidth(lw)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):

        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        ydata = ((height - ydescent) / 2.) * np.ones(xdata.shape, float)
        legline = Line2D(xdata, ydata)

        self.update_prop(legline, orig_handle, legend)
        legline.set_transform(trans)

        return [legline]


class HandlerRegularPolyCollection(HandlerNpointsYoffsets):
    """
    Handler for `.RegularPolyCollections`.
    """
    def __init__(self, yoffsets=None, sizes=None, **kw):
        HandlerNpointsYoffsets.__init__(self, yoffsets=yoffsets, **kw)

        self._sizes = sizes

    def get_numpoints(self, legend):
        if self._numpoints is None:
            return legend.scatterpoints
        else:
            return self._numpoints

    def get_sizes(self, legend, orig_handle,
                  xdescent, ydescent, width, height, fontsize):
        if self._sizes is None:
            handle_sizes = orig_handle.get_sizes()
            if not len(handle_sizes):
                handle_sizes = [1]
            size_max = max(handle_sizes) * legend.markerscale ** 2
            size_min = min(handle_sizes) * legend.markerscale ** 2

            numpoints = self.get_numpoints(legend)
            if numpoints < 4:
                sizes = [.5 * (size_max + size_min), size_max,
                         size_min][:numpoints]
            else:
                rng = (size_max - size_min)
                sizes = rng * np.linspace(0, 1, numpoints) + size_min
        else:
            sizes = self._sizes

        return sizes

    def update_prop(self, legend_handle, orig_handle, legend):

        self._update_prop(legend_handle, orig_handle)

        legend_handle.set_figure(legend.figure)
        #legend._set_artist_props(legend_handle)
        legend_handle.set_clip_box(None)
        legend_handle.set_clip_path(None)

    def create_collection(self, orig_handle, sizes, offsets, transOffset):
        p = type(orig_handle)(orig_handle.get_numsides(),
                              rotation=orig_handle.get_rotation(),
                              sizes=sizes,
                              offsets=offsets,
                              transOffset=transOffset,
                              )
        return p

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)

        ydata = self.get_ydata(legend, xdescent, ydescent,
                               width, height, fontsize)

        sizes = self.get_sizes(legend, orig_handle, xdescent, ydescent,
                               width, height, fontsize)

        p = self.create_collection(orig_handle, sizes,
                                   offsets=list(zip(xdata_marker, ydata)),
                                   transOffset=trans)

        self.update_prop(p, orig_handle, legend)
        p._transOffset = trans
        return [p]


class HandlerPathCollection(HandlerRegularPolyCollection):
    """
    Handler for `.PathCollections`, which are used by `~.Axes.scatter`.
    """
    def create_collection(self, orig_handle, sizes, offsets, transOffset):
        p = type(orig_handle)([orig_handle.get_paths()[0]],
                              sizes=sizes,
                              offsets=offsets,
                              transOffset=transOffset,
                              )
        return p


class HandlerCircleCollection(HandlerRegularPolyCollection):
    """
    Handler for `.CircleCollections`.
    """
    def create_collection(self, orig_handle, sizes, offsets, transOffset):
        p = type(orig_handle)(sizes,
                              offsets=offsets,
                              transOffset=transOffset,
                              )
        return p


class HandlerErrorbar(HandlerLine2D):
    """
    Handler for Errorbars.
    """
    def __init__(self, xerr_size=0.5, yerr_size=None,
                 marker_pad=0.3, numpoints=None, **kw):

        self._xerr_size = xerr_size
        self._yerr_size = yerr_size

        HandlerLine2D.__init__(self, marker_pad=marker_pad,
                               numpoints=numpoints, **kw)

    def get_err_size(self, legend, xdescent, ydescent,
                     width, height, fontsize):
        xerr_size = self._xerr_size * fontsize

        if self._yerr_size is None:
            yerr_size = xerr_size
        else:
            yerr_size = self._yerr_size * fontsize

        return xerr_size, yerr_size

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        plotlines, caplines, barlinecols = orig_handle

        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)

        ydata = ((height - ydescent) / 2.) * np.ones(xdata.shape, float)
        legline = Line2D(xdata, ydata)

        xdata_marker = np.asarray(xdata_marker)
        ydata_marker = np.asarray(ydata[:len(xdata_marker)])

        xerr_size, yerr_size = self.get_err_size(legend, xdescent, ydescent,
                                                 width, height, fontsize)

        legline_marker = Line2D(xdata_marker, ydata_marker)

        # when plotlines are None (only errorbars are drawn), we just
        # make legline invisible.
        if plotlines is None:
            legline.set_visible(False)
            legline_marker.set_visible(False)
        else:
            self.update_prop(legline, plotlines, legend)

            legline.set_drawstyle('default')
            legline.set_marker('None')

            self.update_prop(legline_marker, plotlines, legend)
            legline_marker.set_linestyle('None')

            if legend.markerscale != 1:
                newsz = legline_marker.get_markersize() * legend.markerscale
                legline_marker.set_markersize(newsz)

        handle_barlinecols = []
        handle_caplines = []

        if orig_handle.has_xerr:
            verts = [ ((x - xerr_size, y), (x + xerr_size, y))
                      for x, y in zip(xdata_marker, ydata_marker)]
            coll = mcoll.LineCollection(verts)
            self.update_prop(coll, barlinecols[0], legend)
            handle_barlinecols.append(coll)

            if caplines:
                capline_left = Line2D(xdata_marker - xerr_size, ydata_marker)
                capline_right = Line2D(xdata_marker + xerr_size, ydata_marker)
                self.update_prop(capline_left, caplines[0], legend)
                self.update_prop(capline_right, caplines[0], legend)
                capline_left.set_marker("|")
                capline_right.set_marker("|")

                handle_caplines.append(capline_left)
                handle_caplines.append(capline_right)

        if orig_handle.has_yerr:
            verts = [ ((x, y - yerr_size), (x, y + yerr_size))
                      for x, y in zip(xdata_marker, ydata_marker)]
            coll = mcoll.LineCollection(verts)
            self.update_prop(coll, barlinecols[0], legend)
            handle_barlinecols.append(coll)

            if caplines:
                capline_left = Line2D(xdata_marker, ydata_marker - yerr_size)
                capline_right = Line2D(xdata_marker, ydata_marker + yerr_size)
                self.update_prop(capline_left, caplines[0], legend)
                self.update_prop(capline_right, caplines[0], legend)
                capline_left.set_marker("_")
                capline_right.set_marker("_")

                handle_caplines.append(capline_left)
                handle_caplines.append(capline_right)

        artists = []
        artists.extend(handle_barlinecols)
        artists.extend(handle_caplines)
        artists.append(legline)
        artists.append(legline_marker)

        for artist in artists:
            artist.set_transform(trans)

        return artists


class HandlerStem(HandlerNpointsYoffsets):
    """
    Handler for plots produced by `~.Axes.stem`.
    """
    def __init__(self, marker_pad=0.3, numpoints=None,
                 bottom=None, yoffsets=None, **kw):
        """
        Parameters
        ----------
        marker_pad : float
            Padding between points in legend entry. Default is 0.3.

        numpoints : int, optional
            Number of points to show in legend entry.

        bottom : float, optional

        yoffsets : array of floats, optional
            Length *numpoints* list of y offsets for each point in
            legend entry.

        Notes
        -----
        Any other keyword arguments are given to `HandlerNpointsYoffsets`.
        """

        HandlerNpointsYoffsets.__init__(self, marker_pad=marker_pad,
                                        numpoints=numpoints,
                                        yoffsets=yoffsets,
                                        **kw)
        self._bottom = bottom

    def get_ydata(self, legend, xdescent, ydescent, width, height, fontsize):
        if self._yoffsets is None:
            ydata = height * (0.5 * legend._scatteryoffsets + 0.5)
        else:
            ydata = height * np.asarray(self._yoffsets)

        return ydata

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        markerline, stemlines, baseline = orig_handle

        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)

        ydata = self.get_ydata(legend, xdescent, ydescent,
                               width, height, fontsize)

        if self._bottom is None:
            bottom = 0.
        else:
            bottom = self._bottom

        leg_markerline = Line2D(xdata_marker, ydata[:len(xdata_marker)])
        self.update_prop(leg_markerline, markerline, legend)

        leg_stemlines = []
        for thisx, thisy in zip(xdata_marker, ydata):
            l = Line2D([thisx, thisx], [bottom, thisy])
            leg_stemlines.append(l)

        for lm, m in zip(leg_stemlines, stemlines):
            self.update_prop(lm, m, legend)

        leg_baseline = Line2D([np.min(xdata), np.max(xdata)],
                              [bottom, bottom])

        self.update_prop(leg_baseline, baseline, legend)

        artists = [leg_markerline]
        artists.extend(leg_stemlines)
        artists.append(leg_baseline)

        for artist in artists:
            artist.set_transform(trans)

        return artists


class HandlerTuple(HandlerBase):
    """
    Handler for Tuple.

    Additional kwargs are passed through to `HandlerBase`.

    Parameters
    ----------
    ndivide : int, optional
        The number of sections to divide the legend area into. If None,
        use the length of the input tuple. Default is 1.


    pad : float, optional
        If None, fall back to ``legend.borderpad`` as the default.
        In units of fraction of font size. Default is None.
    """
    def __init__(self, ndivide=1, pad=None, **kwargs):

        self._ndivide = ndivide
        self._pad = pad
        HandlerBase.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        handler_map = legend.get_legend_handler_map()

        if self._ndivide is None:
            ndivide = len(orig_handle)
        else:
            ndivide = self._ndivide

        if self._pad is None:
            pad = legend.borderpad * fontsize
        else:
            pad = self._pad * fontsize

        if ndivide > 1:
            width = (width - pad * (ndivide - 1)) / ndivide

        xds_cycle = cycle(xdescent - (width + pad) * np.arange(ndivide))

        a_list = []
        for handle1 in orig_handle:
            handler = legend.get_legend_handler(handler_map, handle1)
            _a_list = handler.create_artists(
                legend, handle1,
                next(xds_cycle), ydescent, width, height, fontsize, trans)
            a_list.extend(_a_list)

        return a_list


class HandlerPolyCollection(HandlerBase):
    """
    Handler for `.PolyCollection` used in `~.Axes.fill_between` and `~.Axes.stackplot`.
    """
    def _update_prop(self, legend_handle, orig_handle):
        def first_color(colors):
            if colors is None:
                return None
            colors = mcolors.to_rgba_array(colors)
            if len(colors):
                return colors[0]
            else:
                return "none"

        def get_first(prop_array):
            if len(prop_array):
                return prop_array[0]
            else:
                return None
        edgecolor = getattr(orig_handle, '_original_edgecolor',
                            orig_handle.get_edgecolor())
        legend_handle.set_edgecolor(first_color(edgecolor))
        facecolor = getattr(orig_handle, '_original_facecolor',
                            orig_handle.get_facecolor())
        legend_handle.set_facecolor(first_color(facecolor))
        legend_handle.set_fill(orig_handle.get_fill())
        legend_handle.set_hatch(orig_handle.get_hatch())
        legend_handle.set_linewidth(get_first(orig_handle.get_linewidths()))
        legend_handle.set_linestyle(get_first(orig_handle.get_linestyles()))
        legend_handle.set_transform(get_first(orig_handle.get_transforms()))
        legend_handle.set_figure(orig_handle.get_figure())
        legend_handle.set_alpha(orig_handle.get_alpha())

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        p = Rectangle(xy=(-xdescent, -ydescent),
                      width=width, height=height)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]
