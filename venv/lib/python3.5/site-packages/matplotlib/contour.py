"""
These are classes to support contour plotting and labelling for the Axes class.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

import warnings
import matplotlib as mpl
import numpy as np
from numpy import ma
import matplotlib._contour as _contour
import matplotlib.path as mpath
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.collections as mcoll
import matplotlib.font_manager as font_manager
import matplotlib.text as text
import matplotlib.cbook as cbook
import matplotlib.mathtext as mathtext
import matplotlib.patches as mpatches
import matplotlib.texmanager as texmanager
import matplotlib.transforms as mtransforms

# Import needed for adding manual selection capability to clabel
from matplotlib.blocking_input import BlockingContourLabeler

# We can't use a single line collection for contour because a line
# collection can have only a single line style, and we want to be able to have
# dashed negative contours, for example, and solid positive contours.
# We could use a single polygon collection for filled contours, but it
# seems better to keep line and filled contours similar, with one collection
# per level.


class ClabelText(text.Text):
    """
    Unlike the ordinary text, the get_rotation returns an updated
    angle in the pixel coordinate assuming that the input rotation is
    an angle in data coordinate (or whatever transform set).
    """
    def get_rotation(self):
        angle = text.Text.get_rotation(self)
        trans = self.get_transform()
        x, y = self.get_position()
        new_angles = trans.transform_angles(np.array([angle]),
                                            np.array([[x, y]]))
        return new_angles[0]


class ContourLabeler(object):
    """Mixin to provide labelling capability to ContourSet"""

    def clabel(self, *args, **kwargs):
        """
        Label a contour plot.

        Call signature::

          clabel(cs, **kwargs)

        Adds labels to line contours in *cs*, where *cs* is a
        :class:`~matplotlib.contour.ContourSet` object returned by
        contour.

        ::

          clabel(cs, v, **kwargs)

        only labels contours listed in *v*.

        Parameters
        ----------
        fontsize : string or float, optional
            Size in points or relative size e.g., 'smaller', 'x-large'.
            See `Text.set_size` for accepted string values.

        colors :
            Color of each label

            - if *None*, the color of each label matches the color of
              the corresponding contour

            - if one string color, e.g., *colors* = 'r' or *colors* =
              'red', all labels will be plotted in this color

            - if a tuple of matplotlib color args (string, float, rgb, etc),
              different labels will be plotted in different colors in the order
              specified

        inline : bool, optional
            If ``True`` the underlying contour is removed where the label is
            placed. Default is ``True``.

        inline_spacing : float, optional
            Space in pixels to leave on each side of label when
            placing inline. Defaults to 5.

            This spacing will be exact for labels at locations where the
            contour is straight, less so for labels on curved contours.

        fmt : string or dict, optional
            A format string for the label. Default is '%1.3f'

            Alternatively, this can be a dictionary matching contour
            levels with arbitrary strings to use for each contour level
            (i.e., fmt[level]=string), or it can be any callable, such
            as a :class:`~matplotlib.ticker.Formatter` instance, that
            returns a string when called with a numeric contour level.

        manual : bool or iterable, optional
            If ``True``, contour labels will be placed manually using
            mouse clicks. Click the first button near a contour to
            add a label, click the second button (or potentially both
            mouse buttons at once) to finish adding labels. The third
            button can be used to remove the last label added, but
            only if labels are not inline. Alternatively, the keyboard
            can be used to select label locations (enter to end label
            placement, delete or backspace act like the third mouse button,
            and any other key will select a label location).

            *manual* can also be an iterable object of x,y tuples.
            Contour labels will be created as if mouse is clicked at each
            x,y positions.

        rightside_up : bool, optional
            If ``True``, label rotations will always be plus
            or minus 90 degrees from level. Default is ``True``.

        use_clabeltext : bool, optional
            If ``True``, `ClabelText` class (instead of `Text`) is used to
            create labels. `ClabelText` recalculates rotation angles
            of texts during the drawing time, therefore this can be used if
            aspect of the axes changes. Default is ``False``.
        """

        """
        NOTES on how this all works:

        clabel basically takes the input arguments and uses them to
        add a list of "label specific" attributes to the ContourSet
        object.  These attributes are all of the form label* and names
        should be fairly self explanatory.

        Once these attributes are set, clabel passes control to the
        labels method (case of automatic label placement) or
        `BlockingContourLabeler` (case of manual label placement).
        """

        fontsize = kwargs.get('fontsize', None)
        inline = kwargs.get('inline', 1)
        inline_spacing = kwargs.get('inline_spacing', 5)
        self.labelFmt = kwargs.get('fmt', '%1.3f')
        _colors = kwargs.get('colors', None)

        self._use_clabeltext = kwargs.get('use_clabeltext', False)

        # Detect if manual selection is desired and remove from argument list
        self.labelManual = kwargs.get('manual', False)

        self.rightside_up = kwargs.get('rightside_up', True)
        if len(args) == 0:
            levels = self.levels
            indices = list(xrange(len(self.cvalues)))
        elif len(args) == 1:
            levlabs = list(args[0])
            indices, levels = [], []
            for i, lev in enumerate(self.levels):
                if lev in levlabs:
                    indices.append(i)
                    levels.append(lev)
            if len(levels) < len(levlabs):
                raise ValueError("Specified levels {} don't match available "
                                 "levels {}".format(levlabs, self.levels))
        else:
            raise TypeError("Illegal arguments to clabel, see help(clabel)")
        self.labelLevelList = levels
        self.labelIndiceList = indices

        self.labelFontProps = font_manager.FontProperties()
        self.labelFontProps.set_size(fontsize)
        font_size_pts = self.labelFontProps.get_size_in_points()
        self.labelFontSizeList = [font_size_pts] * len(levels)

        if _colors is None:
            self.labelMappable = self
            self.labelCValueList = np.take(self.cvalues, self.labelIndiceList)
        else:
            cmap = colors.ListedColormap(_colors, N=len(self.labelLevelList))
            self.labelCValueList = list(xrange(len(self.labelLevelList)))
            self.labelMappable = cm.ScalarMappable(cmap=cmap,
                                                   norm=colors.NoNorm())

        self.labelXYs = []

        if cbook.iterable(self.labelManual):
            for x, y in self.labelManual:
                self.add_label_near(x, y, inline,
                                    inline_spacing)

        elif self.labelManual:
            print('Select label locations manually using first mouse button.')
            print('End manual selection with second mouse button.')
            if not inline:
                print('Remove last label by clicking third mouse button.')

            blocking_contour_labeler = BlockingContourLabeler(self)
            blocking_contour_labeler(inline, inline_spacing)
        else:
            self.labels(inline, inline_spacing)

        # Hold on to some old attribute names.  These are deprecated and will
        # be removed in the near future (sometime after 2008-08-01), but
        # keeping for now for backwards compatibility
        self.cl = self.labelTexts
        self.cl_xy = self.labelXYs
        self.cl_cvalues = self.labelCValues

        self.labelTextsList = cbook.silent_list('text.Text', self.labelTexts)
        return self.labelTextsList

    def print_label(self, linecontour, labelwidth):
        "Return *False* if contours are too short for a label."
        return (len(linecontour) > 10 * labelwidth
                or (np.ptp(linecontour, axis=0) > 1.2 * labelwidth).any())

    def too_close(self, x, y, lw):
        "Return *True* if a label is already near this location."
        for loc in self.labelXYs:
            d = np.sqrt((x - loc[0]) ** 2 + (y - loc[1]) ** 2)
            if d < 1.2 * lw:
                return True
        return False

    def get_label_coords(self, distances, XX, YY, ysize, lw):
        """
        Return x, y, and the index of a label location.

        Labels are plotted at a location with the smallest
        deviation of the contour from a straight line
        unless there is another label nearby, in which case
        the next best place on the contour is picked up.
        If all such candidates are rejected, the beginning
        of the contour is chosen.
        """
        hysize = int(ysize / 2)
        adist = np.argsort(distances)

        for ind in adist:
            x, y = XX[ind][hysize], YY[ind][hysize]
            if self.too_close(x, y, lw):
                continue
            return x, y, ind

        ind = adist[0]
        x, y = XX[ind][hysize], YY[ind][hysize]
        return x, y, ind

    def get_label_width(self, lev, fmt, fsize):
        """
        Return the width of the label in points.
        """
        if not isinstance(lev, six.string_types):
            lev = self.get_text(lev, fmt)

        lev, ismath = text.Text.is_math_text(lev)
        if ismath == 'TeX':
            if not hasattr(self, '_TeX_manager'):
                self._TeX_manager = texmanager.TexManager()
            lw, _, _ = self._TeX_manager.get_text_width_height_descent(lev,
                                                                       fsize)
        elif ismath:
            if not hasattr(self, '_mathtext_parser'):
                self._mathtext_parser = mathtext.MathTextParser('bitmap')
            img, _ = self._mathtext_parser.parse(lev, dpi=72,
                                                 prop=self.labelFontProps)
            lw = img.get_width()  # at dpi=72, the units are PostScript points
        else:
            # width is much less than "font size"
            lw = (len(lev)) * fsize * 0.6

        return lw

    @cbook.deprecated("2.2")
    def get_real_label_width(self, lev, fmt, fsize):
        """
        This computes actual onscreen label width.
        This uses some black magic to determine onscreen extent of non-drawn
        label.  This magic may not be very robust.

        This method is not being used, and may be modified or removed.
        """
        # Find middle of axes
        xx = np.mean(np.asarray(self.ax.axis()).reshape(2, 2), axis=1)

        # Temporarily create text object
        t = text.Text(xx[0], xx[1])
        self.set_label_props(t, self.get_text(lev, fmt), 'k')

        # Some black magic to get onscreen extent
        # NOTE: This will only work for already drawn figures, as the canvas
        # does not have a renderer otherwise.  This is the reason this function
        # can't be integrated into the rest of the code.
        bbox = t.get_window_extent(renderer=self.ax.figure.canvas.renderer)

        # difference in pixel extent of image
        lw = np.diff(bbox.corners()[0::2, 0])[0]

        return lw

    def set_label_props(self, label, text, color):
        "set the label properties - color, fontsize, text"
        label.set_text(text)
        label.set_color(color)
        label.set_fontproperties(self.labelFontProps)
        label.set_clip_box(self.ax.bbox)

    def get_text(self, lev, fmt):
        "get the text of the label"
        if isinstance(lev, six.string_types):
            return lev
        else:
            if isinstance(fmt, dict):
                return fmt.get(lev, '%1.3f')
            elif callable(fmt):
                return fmt(lev)
            else:
                return fmt % lev

    def locate_label(self, linecontour, labelwidth):
        """
        Find good place to draw a label (relatively flat part of the contour).
        """

        # Number of contour points
        nsize = len(linecontour)
        if labelwidth > 1:
            xsize = int(np.ceil(nsize / labelwidth))
        else:
            xsize = 1
        if xsize == 1:
            ysize = nsize
        else:
            ysize = int(labelwidth)

        XX = np.resize(linecontour[:, 0], (xsize, ysize))
        YY = np.resize(linecontour[:, 1], (xsize, ysize))
        # I might have fouled up the following:
        yfirst = YY[:, :1]
        ylast = YY[:, -1:]
        xfirst = XX[:, :1]
        xlast = XX[:, -1:]
        s = (yfirst - YY) * (xlast - xfirst) - (xfirst - XX) * (ylast - yfirst)
        L = np.hypot(xlast - xfirst, ylast - yfirst)
        # Ignore warning that divide by zero throws, as this is a valid option
        with np.errstate(divide='ignore', invalid='ignore'):
            dist = np.sum(np.abs(s) / L, axis=-1)
        x, y, ind = self.get_label_coords(dist, XX, YY, ysize, labelwidth)

        # There must be a more efficient way...
        lc = [tuple(l) for l in linecontour]
        dind = lc.index((x, y))

        return x, y, dind

    def calc_label_rot_and_inline(self, slc, ind, lw, lc=None, spacing=5):
        """
        This function calculates the appropriate label rotation given
        the linecontour coordinates in screen units, the index of the
        label location and the label width.

        It will also break contour and calculate inlining if *lc* is
        not empty (lc defaults to the empty list if None).  *spacing*
        is the space around the label in pixels to leave empty.

        Do both of these tasks at once to avoid calculating path lengths
        multiple times, which is relatively costly.

        The method used here involves calculating the path length
        along the contour in pixel coordinates and then looking
        approximately label width / 2 away from central point to
        determine rotation and then to break contour if desired.
        """

        if lc is None:
            lc = []
        # Half the label width
        hlw = lw / 2.0

        # Check if closed and, if so, rotate contour so label is at edge
        closed = _is_closed_polygon(slc)
        if closed:
            slc = np.r_[slc[ind:-1], slc[:ind + 1]]

            if len(lc):  # Rotate lc also if not empty
                lc = np.r_[lc[ind:-1], lc[:ind + 1]]

            ind = 0

        # Calculate path lengths
        pl = np.zeros(slc.shape[0], dtype=float)
        dx = np.diff(slc, axis=0)
        pl[1:] = np.cumsum(np.hypot(dx[:, 0], dx[:, 1]))
        pl = pl - pl[ind]

        # Use linear interpolation to get points around label
        xi = np.array([-hlw, hlw])
        if closed:  # Look at end also for closed contours
            dp = np.array([pl[-1], 0])
        else:
            dp = np.zeros_like(xi)

        # Get angle of vector between the two ends of the label - must be
        # calculated in pixel space for text rotation to work correctly.
        (dx,), (dy,) = (np.diff(np.interp(dp + xi, pl, slc_col))
                        for slc_col in slc.T)
        rotation = np.rad2deg(np.arctan2(dy, dx))

        if self.rightside_up:
            # Fix angle so text is never upside-down
            rotation = (rotation + 90) % 180 - 90

        # Break contour if desired
        nlc = []
        if len(lc):
            # Expand range by spacing
            xi = dp + xi + np.array([-spacing, spacing])

            # Get (integer) indices near points of interest; use -1 as marker
            # for out of bounds.
            I = np.interp(xi, pl, np.arange(len(pl)), left=-1, right=-1)
            I = [np.floor(I[0]).astype(int), np.ceil(I[1]).astype(int)]
            if I[0] != -1:
                xy1 = [np.interp(xi[0], pl, lc_col) for lc_col in lc.T]
            if I[1] != -1:
                xy2 = [np.interp(xi[1], pl, lc_col) for lc_col in lc.T]

            # Actually break contours
            if closed:
                # This will remove contour if shorter than label
                if all(i != -1 for i in I):
                    nlc.append(np.row_stack([xy2, lc[I[1]:I[0]+1], xy1]))
            else:
                # These will remove pieces of contour if they have length zero
                if I[0] != -1:
                    nlc.append(np.row_stack([lc[:I[0]+1], xy1]))
                if I[1] != -1:
                    nlc.append(np.row_stack([xy2, lc[I[1]:]]))

            # The current implementation removes contours completely
            # covered by labels.  Uncomment line below to keep
            # original contour if this is the preferred behavior.
            # if not len(nlc): nlc = [ lc ]

        return rotation, nlc

    def _get_label_text(self, x, y, rotation):
        dx, dy = self.ax.transData.inverted().transform_point((x, y))
        t = text.Text(dx, dy, rotation=rotation,
                      horizontalalignment='center',
                      verticalalignment='center')
        return t

    def _get_label_clabeltext(self, x, y, rotation):
        # x, y, rotation is given in pixel coordinate. Convert them to
        # the data coordinate and create a label using ClabelText
        # class. This way, the roation of the clabel is along the
        # contour line always.
        transDataInv = self.ax.transData.inverted()
        dx, dy = transDataInv.transform_point((x, y))
        drotation = transDataInv.transform_angles(np.array([rotation]),
                                                  np.array([[x, y]]))
        t = ClabelText(dx, dy, rotation=drotation[0],
                       horizontalalignment='center',
                       verticalalignment='center')

        return t

    def _add_label(self, t, x, y, lev, cvalue):
        color = self.labelMappable.to_rgba(cvalue, alpha=self.alpha)

        _text = self.get_text(lev, self.labelFmt)
        self.set_label_props(t, _text, color)
        self.labelTexts.append(t)
        self.labelCValues.append(cvalue)
        self.labelXYs.append((x, y))

        # Add label to plot here - useful for manual mode label selection
        self.ax.add_artist(t)

    def add_label(self, x, y, rotation, lev, cvalue):
        """
        Add contour label using :class:`~matplotlib.text.Text` class.
        """

        t = self._get_label_text(x, y, rotation)
        self._add_label(t, x, y, lev, cvalue)

    def add_label_clabeltext(self, x, y, rotation, lev, cvalue):
        """
        Add contour label using :class:`ClabelText` class.
        """
        # x, y, rotation is given in pixel coordinate. Convert them to
        # the data coordinate and create a label using ClabelText
        # class. This way, the roation of the clabel is along the
        # contour line always.

        t = self._get_label_clabeltext(x, y, rotation)
        self._add_label(t, x, y, lev, cvalue)

    def add_label_near(self, x, y, inline=True, inline_spacing=5,
                       transform=None):
        """
        Add a label near the point (x, y). If transform is None
        (default), (x, y) is in data coordinates; if transform is
        False, (x, y) is in display coordinates; otherwise, the
        specified transform will be used to translate (x, y) into
        display coordinates.

        *inline*:
          controls whether the underlying contour is removed or
          not. Default is *True*.

        *inline_spacing*:
          space in pixels to leave on each side of label when
          placing inline.  Defaults to 5.  This spacing will be
          exact for labels at locations where the contour is
          straight, less so for labels on curved contours.
        """

        if transform is None:
            transform = self.ax.transData

        if transform:
            x, y = transform.transform_point((x, y))

        # find the nearest contour _in screen units_
        conmin, segmin, imin, xmin, ymin = self.find_nearest_contour(
            x, y, self.labelIndiceList)[:5]

        # The calc_label_rot_and_inline routine requires that (xmin,ymin)
        # be a vertex in the path. So, if it isn't, add a vertex here

        # grab the paths from the collections
        paths = self.collections[conmin].get_paths()
        # grab the correct segment
        active_path = paths[segmin]
        # grab its vertices
        lc = active_path.vertices
        # sort out where the new vertex should be added data-units
        xcmin = self.ax.transData.inverted().transform_point([xmin, ymin])
        # if there isn't a vertex close enough
        if not np.allclose(xcmin, lc[imin]):
            # insert new data into the vertex list
            lc = np.r_[lc[:imin], np.array(xcmin)[None, :], lc[imin:]]
            # replace the path with the new one
            paths[segmin] = mpath.Path(lc)

        # Get index of nearest level in subset of levels used for labeling
        lmin = self.labelIndiceList.index(conmin)

        # Coordinates of contour
        paths = self.collections[conmin].get_paths()
        lc = paths[segmin].vertices

        # In pixel/screen space
        slc = self.ax.transData.transform(lc)

        # Get label width for rotating labels and breaking contours
        lw = self.get_label_width(self.labelLevelList[lmin],
                                  self.labelFmt, self.labelFontSizeList[lmin])
        # lw is in points.
        lw *= self.ax.figure.dpi / 72.0  # scale to screen coordinates
        # now lw in pixels

        # Figure out label rotation.
        if inline:
            lcarg = lc
        else:
            lcarg = None
        rotation, nlc = self.calc_label_rot_and_inline(
            slc, imin, lw, lcarg,
            inline_spacing)

        self.add_label(xmin, ymin, rotation, self.labelLevelList[lmin],
                       self.labelCValueList[lmin])

        if inline:
            # Remove old, not looping over paths so we can do this up front
            paths.pop(segmin)

            # Add paths if not empty or single point
            for n in nlc:
                if len(n) > 1:
                    paths.append(mpath.Path(n))

    def pop_label(self, index=-1):
        """Defaults to removing last label, but any index can be supplied"""
        self.labelCValues.pop(index)
        t = self.labelTexts.pop(index)
        t.remove()

    def labels(self, inline, inline_spacing):

        if self._use_clabeltext:
            add_label = self.add_label_clabeltext
        else:
            add_label = self.add_label

        for icon, lev, fsize, cvalue in zip(
                self.labelIndiceList, self.labelLevelList,
                self.labelFontSizeList, self.labelCValueList):

            con = self.collections[icon]
            trans = con.get_transform()
            lw = self.get_label_width(lev, self.labelFmt, fsize)
            lw *= self.ax.figure.dpi / 72.0  # scale to screen coordinates
            additions = []
            paths = con.get_paths()
            for segNum, linepath in enumerate(paths):
                lc = linepath.vertices  # Line contour
                slc0 = trans.transform(lc)  # Line contour in screen coords

                # For closed polygons, add extra point to avoid division by
                # zero in print_label and locate_label.  Other than these
                # functions, this is not necessary and should probably be
                # eventually removed.
                if _is_closed_polygon(lc):
                    slc = np.r_[slc0, slc0[1:2, :]]
                else:
                    slc = slc0

                # Check if long enough for a label
                if self.print_label(slc, lw):
                    x, y, ind = self.locate_label(slc, lw)

                    if inline:
                        lcarg = lc
                    else:
                        lcarg = None
                    rotation, new = self.calc_label_rot_and_inline(
                        slc0, ind, lw, lcarg,
                        inline_spacing)

                    # Actually add the label
                    add_label(x, y, rotation, lev, cvalue)

                    # If inline, add new contours
                    if inline:
                        for n in new:
                            # Add path if not empty or single point
                            if len(n) > 1:
                                additions.append(mpath.Path(n))
                else:  # If not adding label, keep old path
                    additions.append(linepath)

            # After looping over all segments on a contour, remove old
            # paths and add new ones if inlining
            if inline:
                del paths[:]
                paths.extend(additions)


def _find_closest_point_on_leg(p1, p2, p0):
    """find closest point to p0 on line segment connecting p1 and p2"""

    # handle degenerate case
    if np.all(p2 == p1):
        d = np.sum((p0 - p1)**2)
        return d, p1

    d21 = p2 - p1
    d01 = p0 - p1

    # project on to line segment to find closest point
    proj = np.dot(d01, d21) / np.dot(d21, d21)
    if proj < 0:
        proj = 0
    if proj > 1:
        proj = 1
    pc = p1 + proj * d21

    # find squared distance
    d = np.sum((pc-p0)**2)

    return d, pc


def _is_closed_polygon(X):
    """
    Tests whether first and last object in a sequence are the same.  These are
    presumably coordinates on a polygonal curve, in which case this function
    tests if that curve is closed.
    """
    return np.all(X[0] == X[-1])


def _find_closest_point_on_path(lc, point):
    """
    lc: coordinates of vertices
    point: coordinates of test point
    """

    # find index of closest vertex for this segment
    ds = np.sum((lc - point[None, :])**2, 1)
    imin = np.argmin(ds)

    dmin = np.inf
    xcmin = None
    legmin = (None, None)

    closed = _is_closed_polygon(lc)

    # build list of legs before and after this vertex
    legs = []
    if imin > 0 or closed:
        legs.append(((imin-1) % len(lc), imin))
    if imin < len(lc) - 1 or closed:
        legs.append((imin, (imin+1) % len(lc)))

    for leg in legs:
        d, xc = _find_closest_point_on_leg(lc[leg[0]], lc[leg[1]], point)
        if d < dmin:
            dmin = d
            xcmin = xc
            legmin = leg

    return (dmin, xcmin, legmin)


class ContourSet(cm.ScalarMappable, ContourLabeler):
    """
    Store a set of contour lines or filled regions.

    User-callable method: clabel

    Attributes
    ----------
    ax:
        The axes object in which the contours are drawn.

    collections:
        A silent_list of LineCollections or PolyCollections.

    levels:
        Contour levels.

    layers:
        Same as levels for line contours; half-way between
        levels for filled contours.  See :meth:`_process_colors`.
    """

    def __init__(self, ax, *args, **kwargs):
        """
        Draw contour lines or filled regions, depending on
        whether keyword arg *filled* is ``False`` (default) or ``True``.

        The first three arguments must be:

          *ax*: axes object.

          *levels*: [level0, level1, ..., leveln]
            A list of floating point numbers indicating the contour
            levels.

          *allsegs*: [level0segs, level1segs, ...]
            List of all the polygon segments for all the *levels*.
            For contour lines ``len(allsegs) == len(levels)``, and for
            filled contour regions ``len(allsegs) = len(levels)-1``. The lists
            should look like::

                level0segs = [polygon0, polygon1, ...]
                polygon0 = array_like [[x0,y0], [x1,y1], ...]

          *allkinds*: *None* or [level0kinds, level1kinds, ...]
            Optional list of all the polygon vertex kinds (code types), as
            described and used in Path. This is used to allow multiply-
            connected paths such as holes within filled polygons.
            If not ``None``, ``len(allkinds) == len(allsegs)``. The lists
            should look like::

                level0kinds = [polygon0kinds, ...]
                polygon0kinds = [vertexcode0, vertexcode1, ...]

            If *allkinds* is not ``None``, usually all polygons for a
            particular contour level are grouped together so that
            ``level0segs = [polygon0]`` and ``level0kinds = [polygon0kinds]``.

        Keyword arguments are as described in the docstring of
        `~.Axes.contour`.
        """
        self.ax = ax
        self.levels = kwargs.pop('levels', None)
        self.filled = kwargs.pop('filled', False)
        self.linewidths = kwargs.pop('linewidths', None)
        self.linestyles = kwargs.pop('linestyles', None)

        self.hatches = kwargs.pop('hatches', [None])

        self.alpha = kwargs.pop('alpha', None)
        self.origin = kwargs.pop('origin', None)
        self.extent = kwargs.pop('extent', None)
        cmap = kwargs.pop('cmap', None)
        self.colors = kwargs.pop('colors', None)
        norm = kwargs.pop('norm', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        self.extend = kwargs.pop('extend', 'neither')
        self.antialiased = kwargs.pop('antialiased', None)
        if self.antialiased is None and self.filled:
            self.antialiased = False  # eliminate artifacts; we are not
                                      # stroking the boundaries.
            # The default for line contours will be taken from
            # the LineCollection default, which uses the
            # rcParams['lines.antialiased']

        self.nchunk = kwargs.pop('nchunk', 0)
        self.locator = kwargs.pop('locator', None)
        if (isinstance(norm, colors.LogNorm)
                or isinstance(self.locator, ticker.LogLocator)):
            self.logscale = True
            if norm is None:
                norm = colors.LogNorm()
            if self.extend is not 'neither':
                raise ValueError('extend kwarg does not work yet with log '
                                 ' scale')
        else:
            self.logscale = False

        if self.origin not in [None, 'lower', 'upper', 'image']:
            raise ValueError("If given, *origin* must be one of [ 'lower' |"
                             " 'upper' | 'image']")
        if self.extent is not None and len(self.extent) != 4:
            raise ValueError("If given, *extent* must be '[ *None* |"
                             " (x0,x1,y0,y1) ]'")
        if self.colors is not None and cmap is not None:
            raise ValueError('Either colors or cmap must be None')
        if self.origin == 'image':
            self.origin = mpl.rcParams['image.origin']

        self._transform = kwargs.pop('transform', None)

        kwargs = self._process_args(*args, **kwargs)
        self._process_levels()

        if self.colors is not None:
            ncolors = len(self.levels)
            if self.filled:
                ncolors -= 1
            i0 = 0

            # Handle the case where colors are given for the extended
            # parts of the contour.
            extend_min = self.extend in ['min', 'both']
            extend_max = self.extend in ['max', 'both']
            use_set_under_over = False
            # if we are extending the lower end, and we've been given enough
            # colors then skip the first color in the resulting cmap. For the
            # extend_max case we don't need to worry about passing more colors
            # than ncolors as ListedColormap will clip.
            total_levels = ncolors + int(extend_min) + int(extend_max)
            if (len(self.colors) == total_levels and
                    any([extend_min, extend_max])):
                use_set_under_over = True
                if extend_min:
                    i0 = 1

            cmap = colors.ListedColormap(self.colors[i0:None], N=ncolors)

            if use_set_under_over:
                if extend_min:
                    cmap.set_under(self.colors[0])
                if extend_max:
                    cmap.set_over(self.colors[-1])

        if self.filled:
            self.collections = cbook.silent_list('mcoll.PathCollection')
        else:
            self.collections = cbook.silent_list('mcoll.LineCollection')
        # label lists must be initialized here
        self.labelTexts = []
        self.labelCValues = []

        kw = {'cmap': cmap}
        if norm is not None:
            kw['norm'] = norm
        # sets self.cmap, norm if needed;
        cm.ScalarMappable.__init__(self, **kw)
        if vmin is not None:
            self.norm.vmin = vmin
        if vmax is not None:
            self.norm.vmax = vmax
        self._process_colors()

        self.allsegs, self.allkinds = self._get_allsegs_and_allkinds()

        if self.filled:
            if self.linewidths is not None:
                warnings.warn('linewidths is ignored by contourf')

            # Lower and upper contour levels.
            lowers, uppers = self._get_lowers_and_uppers()

            # Ensure allkinds can be zipped below.
            if self.allkinds is None:
                self.allkinds = [None] * len(self.allsegs)

            # Default zorder taken from Collection
            zorder = kwargs.pop('zorder', 1)
            for level, level_upper, segs, kinds in \
                    zip(lowers, uppers, self.allsegs, self.allkinds):
                paths = self._make_paths(segs, kinds)

                col = mcoll.PathCollection(
                    paths,
                    antialiaseds=(self.antialiased,),
                    edgecolors='none',
                    alpha=self.alpha,
                    transform=self.get_transform(),
                    zorder=zorder)
                self.ax.add_collection(col, autolim=False)
                self.collections.append(col)
        else:
            tlinewidths = self._process_linewidths()
            self.tlinewidths = tlinewidths
            tlinestyles = self._process_linestyles()
            aa = self.antialiased
            if aa is not None:
                aa = (self.antialiased,)
            # Default zorder taken from LineCollection
            zorder = kwargs.pop('zorder', 2)
            for level, width, lstyle, segs in \
                    zip(self.levels, tlinewidths, tlinestyles, self.allsegs):
                col = mcoll.LineCollection(
                    segs,
                    antialiaseds=aa,
                    linewidths=width,
                    linestyles=[lstyle],
                    alpha=self.alpha,
                    transform=self.get_transform(),
                    zorder=zorder)
                col.set_label('_nolegend_')
                self.ax.add_collection(col, autolim=False)
                self.collections.append(col)

        for col in self.collections:
            col.sticky_edges.x[:] = [self._mins[0], self._maxs[0]]
            col.sticky_edges.y[:] = [self._mins[1], self._maxs[1]]
        self.ax.update_datalim([self._mins, self._maxs])
        self.ax.autoscale_view(tight=True)

        self.changed()  # set the colors

        if kwargs:
            s = ", ".join(map(repr, kwargs))
            warnings.warn('The following kwargs were not used by contour: ' +
                          s)

    def get_transform(self):
        """
        Return the :class:`~matplotlib.transforms.Transform`
        instance used by this ContourSet.
        """
        if self._transform is None:
            self._transform = self.ax.transData
        elif (not isinstance(self._transform, mtransforms.Transform)
              and hasattr(self._transform, '_as_mpl_transform')):
            self._transform = self._transform._as_mpl_transform(self.ax)
        return self._transform

    def __getstate__(self):
        state = self.__dict__.copy()
        # the C object _contour_generator cannot currently be pickled. This
        # isn't a big issue as it is not actually used once the contour has
        # been calculated.
        state['_contour_generator'] = None
        return state

    def legend_elements(self, variable_name='x', str_format=str):
        """
        Return a list of artist and labels suitable for passing through
        to :func:`plt.legend` which represent this ContourSet.

        Args:

            *variable_name*: the string used inside the inequality used
              on the labels

            *str_format*: function used to format the numbers in the labels
        """
        artists = []
        labels = []

        if self.filled:
            lowers, uppers = self._get_lowers_and_uppers()
            n_levels = len(self.collections)

            for i, (collection, lower, upper) in enumerate(
                    zip(self.collections, lowers, uppers)):
                patch = mpatches.Rectangle(
                    (0, 0), 1, 1,
                    facecolor=collection.get_facecolor()[0],
                    hatch=collection.get_hatch(),
                    alpha=collection.get_alpha())
                artists.append(patch)

                lower = str_format(lower)
                upper = str_format(upper)

                if i == 0 and self.extend in ('min', 'both'):
                    labels.append(r'$%s \leq %s$' % (variable_name,
                                                     lower))
                elif i == n_levels - 1 and self.extend in ('max', 'both'):
                    labels.append(r'$%s > %s$' % (variable_name,
                                                  upper))
                else:
                    labels.append(r'$%s < %s \leq %s$' % (lower,
                                                          variable_name,
                                                          upper))
        else:
            for collection, level in zip(self.collections, self.levels):

                patch = mcoll.LineCollection(None)
                patch.update_from(collection)

                artists.append(patch)
                # format the level for insertion into the labels
                level = str_format(level)
                labels.append(r'$%s = %s$' % (variable_name, level))

        return artists, labels

    def _process_args(self, *args, **kwargs):
        """
        Process *args* and *kwargs*; override in derived classes.

        Must set self.levels, self.zmin and self.zmax, and update axes
        limits.
        """
        self.levels = args[0]
        self.allsegs = args[1]
        self.allkinds = len(args) > 2 and args[2] or None
        self.zmax = np.max(self.levels)
        self.zmin = np.min(self.levels)
        self._auto = False

        # Check lengths of levels and allsegs.
        if self.filled:
            if len(self.allsegs) != len(self.levels) - 1:
                raise ValueError('must be one less number of segments as '
                                 'levels')
        else:
            if len(self.allsegs) != len(self.levels):
                raise ValueError('must be same number of segments as levels')

        # Check length of allkinds.
        if (self.allkinds is not None and
                len(self.allkinds) != len(self.allsegs)):
            raise ValueError('allkinds has different length to allsegs')

        # Determine x,y bounds and update axes data limits.
        flatseglist = [s for seg in self.allsegs for s in seg]
        points = np.concatenate(flatseglist, axis=0)
        self._mins = points.min(axis=0)
        self._maxs = points.max(axis=0)

        return kwargs

    def _get_allsegs_and_allkinds(self):
        """
        Override in derived classes to create and return allsegs and allkinds.
        allkinds can be None.
        """
        return self.allsegs, self.allkinds

    def _get_lowers_and_uppers(self):
        """
        Return (lowers,uppers) for filled contours.
        """
        lowers = self._levels[:-1]
        if self.zmin == lowers[0]:
            # Include minimum values in lowest interval
            lowers = lowers.copy()  # so we don't change self._levels
            if self.logscale:
                lowers[0] = 0.99 * self.zmin
            else:
                lowers[0] -= 1
        uppers = self._levels[1:]
        return (lowers, uppers)

    def _make_paths(self, segs, kinds):
        if kinds is not None:
            return [mpath.Path(seg, codes=kind)
                    for seg, kind in zip(segs, kinds)]
        else:
            return [mpath.Path(seg) for seg in segs]

    def changed(self):
        tcolors = [(tuple(rgba),)
                   for rgba in self.to_rgba(self.cvalues, alpha=self.alpha)]
        self.tcolors = tcolors
        hatches = self.hatches * len(tcolors)
        for color, hatch, collection in zip(tcolors, hatches,
                                            self.collections):
            if self.filled:
                collection.set_facecolor(color)
                # update the collection's hatch (may be None)
                collection.set_hatch(hatch)
            else:
                collection.set_color(color)
        for label, cv in zip(self.labelTexts, self.labelCValues):
            label.set_alpha(self.alpha)
            label.set_color(self.labelMappable.to_rgba(cv))
        # add label colors
        cm.ScalarMappable.changed(self)

    def _autolev(self, N):
        """
        Select contour levels to span the data.

        We need two more levels for filled contours than for
        line contours, because for the latter we need to specify
        the lower and upper boundary of each range. For example,
        a single contour boundary, say at z = 0, requires only
        one contour line, but two filled regions, and therefore
        three levels to provide boundaries for both regions.
        """
        if self.locator is None:
            if self.logscale:
                self.locator = ticker.LogLocator()
            else:
                self.locator = ticker.MaxNLocator(N + 1, min_n_ticks=1)

        lev = self.locator.tick_values(self.zmin, self.zmax)
        self._auto = True
        return lev

    def _contour_level_args(self, z, args):
        """
        Determine the contour levels and store in self.levels.
        """
        if self.filled:
            fn = 'contourf'
        else:
            fn = 'contour'
        self._auto = False
        if self.levels is None:
            if len(args) == 0:
                lev = self._autolev(7)
            else:
                level_arg = args[0]
                try:
                    if type(level_arg) == int:
                        lev = self._autolev(level_arg)
                    else:
                        lev = np.asarray(level_arg).astype(np.float64)
                except:
                    raise TypeError(
                        "Last {0} arg must give levels; see help({0})"
                        .format(fn))
            self.levels = lev
        else:
            self.levels = np.asarray(self.levels).astype(np.float64)

        if not self.filled:
            inside = (self.levels > self.zmin) & (self.levels < self.zmax)
            self.levels = self.levels[inside]
            if len(self.levels) == 0:
                self.levels = [self.zmin]
                warnings.warn("No contour levels were found"
                              " within the data range.")

        if self.filled and len(self.levels) < 2:
            raise ValueError("Filled contours require at least 2 levels.")

        if len(self.levels) > 1 and np.min(np.diff(self.levels)) <= 0.0:
            raise ValueError("Contour levels must be increasing")

    def _process_levels(self):
        """
        Assign values to :attr:`layers` based on :attr:`levels`,
        adding extended layers as needed if contours are filled.

        For line contours, layers simply coincide with levels;
        a line is a thin layer.  No extended levels are needed
        with line contours.
        """
        # Make a private _levels to include extended regions; we
        # want to leave the original levels attribute unchanged.
        # (Colorbar needs this even for line contours.)
        self._levels = list(self.levels)

        if self.extend in ('both', 'min'):
            self._levels.insert(0, min(self.levels[0], self.zmin) - 1)
        if self.extend in ('both', 'max'):
            self._levels.append(max(self.levels[-1], self.zmax) + 1)
        self._levels = np.asarray(self._levels)

        if not self.filled:
            self.layers = self.levels
            return

        # layer values are mid-way between levels
        self.layers = 0.5 * (self._levels[:-1] + self._levels[1:])
        # ...except that extended layers must be outside the
        # normed range:
        if self.extend in ('both', 'min'):
            self.layers[0] = -1e150
        if self.extend in ('both', 'max'):
            self.layers[-1] = 1e150

    def _process_colors(self):
        """
        Color argument processing for contouring.

        Note that we base the color mapping on the contour levels
        and layers, not on the actual range of the Z values.  This
        means we don't have to worry about bad values in Z, and we
        always have the full dynamic range available for the selected
        levels.

        The color is based on the midpoint of the layer, except for
        extended end layers.  By default, the norm vmin and vmax
        are the extreme values of the non-extended levels.  Hence,
        the layer color extremes are not the extreme values of
        the colormap itself, but approach those values as the number
        of levels increases.  An advantage of this scheme is that
        line contours, when added to filled contours, take on
        colors that are consistent with those of the filled regions;
        for example, a contour line on the boundary between two
        regions will have a color intermediate between those
        of the regions.

        """
        self.monochrome = self.cmap.monochrome
        if self.colors is not None:
            # Generate integers for direct indexing.
            i0, i1 = 0, len(self.levels)
            if self.filled:
                i1 -= 1
                # Out of range indices for over and under:
                if self.extend in ('both', 'min'):
                    i0 -= 1
                if self.extend in ('both', 'max'):
                    i1 += 1
            self.cvalues = list(range(i0, i1))
            self.set_norm(colors.NoNorm())
        else:
            self.cvalues = self.layers
        self.set_array(self.levels)
        self.autoscale_None()
        if self.extend in ('both', 'max', 'min'):
            self.norm.clip = False

        # self.tcolors are set by the "changed" method

    def _process_linewidths(self):
        linewidths = self.linewidths
        Nlev = len(self.levels)
        if linewidths is None:
            tlinewidths = [(mpl.rcParams['lines.linewidth'],)] * Nlev
        else:
            if not cbook.iterable(linewidths):
                linewidths = [linewidths] * Nlev
            else:
                linewidths = list(linewidths)
                if len(linewidths) < Nlev:
                    nreps = int(np.ceil(Nlev / len(linewidths)))
                    linewidths = linewidths * nreps
                if len(linewidths) > Nlev:
                    linewidths = linewidths[:Nlev]
            tlinewidths = [(w,) for w in linewidths]
        return tlinewidths

    def _process_linestyles(self):
        linestyles = self.linestyles
        Nlev = len(self.levels)
        if linestyles is None:
            tlinestyles = ['solid'] * Nlev
            if self.monochrome:
                neg_ls = mpl.rcParams['contour.negative_linestyle']
                eps = - (self.zmax - self.zmin) * 1e-15
                for i, lev in enumerate(self.levels):
                    if lev < eps:
                        tlinestyles[i] = neg_ls
        else:
            if isinstance(linestyles, six.string_types):
                tlinestyles = [linestyles] * Nlev
            elif cbook.iterable(linestyles):
                tlinestyles = list(linestyles)
                if len(tlinestyles) < Nlev:
                    nreps = int(np.ceil(Nlev / len(linestyles)))
                    tlinestyles = tlinestyles * nreps
                if len(tlinestyles) > Nlev:
                    tlinestyles = tlinestyles[:Nlev]
            else:
                raise ValueError("Unrecognized type for linestyles kwarg")
        return tlinestyles

    def get_alpha(self):
        """returns alpha to be applied to all ContourSet artists"""
        return self.alpha

    def set_alpha(self, alpha):
        """sets alpha for all ContourSet artists"""
        self.alpha = alpha
        self.changed()

    def find_nearest_contour(self, x, y, indices=None, pixel=True):
        """
        Finds contour that is closest to a point.  Defaults to
        measuring distance in pixels (screen space - useful for manual
        contour labeling), but this can be controlled via a keyword
        argument.

        Returns a tuple containing the contour, segment, index of
        segment, x & y of segment point and distance to minimum point.

        Optional keyword arguments:

          *indices*:
            Indexes of contour levels to consider when looking for
            nearest point.  Defaults to using all levels.

          *pixel*:
            If *True*, measure distance in pixel space, if not, measure
            distance in axes space.  Defaults to *True*.

        """

        # This function uses a method that is probably quite
        # inefficient based on converting each contour segment to
        # pixel coordinates and then comparing the given point to
        # those coordinates for each contour.  This will probably be
        # quite slow for complex contours, but for normal use it works
        # sufficiently well that the time is not noticeable.
        # Nonetheless, improvements could probably be made.

        if indices is None:
            indices = list(xrange(len(self.levels)))

        dmin = np.inf
        conmin = None
        segmin = None
        xmin = None
        ymin = None

        point = np.array([x, y])

        for icon in indices:
            con = self.collections[icon]
            trans = con.get_transform()
            paths = con.get_paths()

            for segNum, linepath in enumerate(paths):
                lc = linepath.vertices
                # transfer all data points to screen coordinates if desired
                if pixel:
                    lc = trans.transform(lc)

                d, xc, leg = _find_closest_point_on_path(lc, point)
                if d < dmin:
                    dmin = d
                    conmin = icon
                    segmin = segNum
                    imin = leg[1]
                    xmin = xc[0]
                    ymin = xc[1]

        return (conmin, segmin, imin, xmin, ymin, dmin)


class QuadContourSet(ContourSet):
    """
    Create and store a set of contour lines or filled regions.

    User-callable method: :meth:`clabel`

    Attributes
    ----------
    ax:
        The axes object in which the contours are drawn.

    collections:
        A silent_list of LineCollections or PolyCollections.

    levels:
        Contour levels.

    layers:
        Same as levels for line contours; half-way between
        levels for filled contours. See :meth:`_process_colors` method.
    """

    def _process_args(self, *args, **kwargs):
        """
        Process args and kwargs.
        """
        if isinstance(args[0], QuadContourSet):
            if self.levels is None:
                self.levels = args[0].levels
            self.zmin = args[0].zmin
            self.zmax = args[0].zmax
            self._corner_mask = args[0]._corner_mask
            contour_generator = args[0]._contour_generator
            self._mins = args[0]._mins
            self._maxs = args[0]._maxs
        else:
            self._corner_mask = kwargs.pop('corner_mask', None)
            if self._corner_mask is None:
                self._corner_mask = mpl.rcParams['contour.corner_mask']

            x, y, z = self._contour_args(args, kwargs)

            _mask = ma.getmask(z)
            if _mask is ma.nomask or not _mask.any():
                _mask = None

            contour_generator = _contour.QuadContourGenerator(
                x, y, z.filled(), _mask, self._corner_mask, self.nchunk)

            t = self.get_transform()

            # if the transform is not trans data, and some part of it
            # contains transData, transform the xs and ys to data coordinates
            if (t != self.ax.transData and
                    any(t.contains_branch_seperately(self.ax.transData))):
                trans_to_data = t - self.ax.transData
                pts = (np.vstack([x.flat, y.flat]).T)
                transformed_pts = trans_to_data.transform(pts)
                x = transformed_pts[..., 0]
                y = transformed_pts[..., 1]

            self._mins = [ma.min(x), ma.min(y)]
            self._maxs = [ma.max(x), ma.max(y)]

        self._contour_generator = contour_generator

        return kwargs

    def _get_allsegs_and_allkinds(self):
        """Compute ``allsegs`` and ``allkinds`` using C extension."""
        allsegs = []
        if self.filled:
            lowers, uppers = self._get_lowers_and_uppers()
            allkinds = []
            for level, level_upper in zip(lowers, uppers):
                vertices, kinds = \
                    self._contour_generator.create_filled_contour(
                        level, level_upper)
                allsegs.append(vertices)
                allkinds.append(kinds)
        else:
            allkinds = None
            for level in self.levels:
                vertices = self._contour_generator.create_contour(level)
                allsegs.append(vertices)
        return allsegs, allkinds

    def _contour_args(self, args, kwargs):
        if self.filled:
            fn = 'contourf'
        else:
            fn = 'contour'
        Nargs = len(args)
        if Nargs <= 2:
            z = ma.asarray(args[0], dtype=np.float64)
            x, y = self._initialize_x_y(z)
            args = args[1:]
        elif Nargs <= 4:
            x, y, z = self._check_xyz(args[:3], kwargs)
            args = args[3:]
        else:
            raise TypeError("Too many arguments to %s; see help(%s)" %
                            (fn, fn))
        z = ma.masked_invalid(z, copy=False)
        self.zmax = float(z.max())
        self.zmin = float(z.min())
        if self.logscale and self.zmin <= 0:
            z = ma.masked_where(z <= 0, z)
            warnings.warn('Log scale: values of z <= 0 have been masked')
            self.zmin = float(z.min())
        self._contour_level_args(z, args)
        return (x, y, z)

    def _check_xyz(self, args, kwargs):
        """
        For functions like contour, check that the dimensions
        of the input arrays match; if x and y are 1D, convert
        them to 2D using meshgrid.

        Possible change: I think we should make and use an ArgumentError
        Exception class (here and elsewhere).
        """
        x, y = args[:2]
        kwargs = self.ax._process_unit_info(xdata=x, ydata=y, kwargs=kwargs)
        x = self.ax.convert_xunits(x)
        y = self.ax.convert_yunits(y)

        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = ma.asarray(args[2], dtype=np.float64)

        if z.ndim != 2:
            raise TypeError("Input z must be a 2D array.")
        elif z.shape[0] < 2 or z.shape[1] < 2:
            raise TypeError("Input z must be at least a 2x2 array.")
        else:
            Ny, Nx = z.shape

        if x.ndim != y.ndim:
            raise TypeError("Number of dimensions of x and y should match.")

        if x.ndim == 1:

            nx, = x.shape
            ny, = y.shape

            if nx != Nx:
                raise TypeError("Length of x must be number of columns in z.")

            if ny != Ny:
                raise TypeError("Length of y must be number of rows in z.")

            x, y = np.meshgrid(x, y)

        elif x.ndim == 2:

            if x.shape != z.shape:
                raise TypeError("Shape of x does not match that of z: found "
                                "{0} instead of {1}.".format(x.shape, z.shape))

            if y.shape != z.shape:
                raise TypeError("Shape of y does not match that of z: found "
                                "{0} instead of {1}.".format(y.shape, z.shape))
        else:
            raise TypeError("Inputs x and y must be 1D or 2D.")

        return x, y, z

    def _initialize_x_y(self, z):
        """
        Return X, Y arrays such that contour(Z) will match imshow(Z)
        if origin is not None.
        The center of pixel Z[i,j] depends on origin:
        if origin is None, x = j, y = i;
        if origin is 'lower', x = j + 0.5, y = i + 0.5;
        if origin is 'upper', x = j + 0.5, y = Nrows - i - 0.5
        If extent is not None, x and y will be scaled to match,
        as in imshow.
        If origin is None and extent is not None, then extent
        will give the minimum and maximum values of x and y.
        """
        if z.ndim != 2:
            raise TypeError("Input must be a 2D array.")
        elif z.shape[0] < 2 or z.shape[1] < 2:
            raise TypeError("Input z must be at least a 2x2 array.")
        else:
            Ny, Nx = z.shape
        if self.origin is None:  # Not for image-matching.
            if self.extent is None:
                return np.meshgrid(np.arange(Nx), np.arange(Ny))
            else:
                x0, x1, y0, y1 = self.extent
                x = np.linspace(x0, x1, Nx)
                y = np.linspace(y0, y1, Ny)
                return np.meshgrid(x, y)
        # Match image behavior:
        if self.extent is None:
            x0, x1, y0, y1 = (0, Nx, 0, Ny)
        else:
            x0, x1, y0, y1 = self.extent
        dx = (x1 - x0) / Nx
        dy = (y1 - y0) / Ny
        x = x0 + (np.arange(Nx) + 0.5) * dx
        y = y0 + (np.arange(Ny) + 0.5) * dy
        if self.origin == 'upper':
            y = y[::-1]
        return np.meshgrid(x, y)

    _contour_doc = """
        Plot contours.

        :func:`~matplotlib.pyplot.contour` and
        :func:`~matplotlib.pyplot.contourf` draw contour lines and
        filled contours, respectively.  Except as noted, function
        signatures and return values are the same for both versions.

        :func:`~matplotlib.pyplot.contourf` differs from the MATLAB
        version in that it does not draw the polygon edges.
        To draw edges, add line contours with
        calls to :func:`~matplotlib.pyplot.contour`.


        Call signatures::

          contour(Z)

        make a contour plot of an array *Z*. The level values are chosen
        automatically.

        ::

          contour(X,Y,Z)

        *X*, *Y* specify the (x, y) coordinates of the surface

        ::

          contour(Z,N)
          contour(X,Y,Z,N)

        contour up to *N+1* automatically chosen contour levels
        (*N* intervals).

        ::

          contour(Z,V)
          contour(X,Y,Z,V)

        draw contour lines at the values specified in sequence *V*,
        which must be in increasing order.

        ::

          contourf(..., V)

        fill the ``len(V)-1`` regions between the values in *V*,
        which must be in increasing order.

        ::

          contour(Z, **kwargs)

        Use keyword args to control colors, linewidth, origin, cmap ... see
        below for more details.

        *X* and *Y* must both be 2-D with the same shape as *Z*, or they
        must both be 1-D such that ``len(X)`` is the number of columns in
        *Z* and ``len(Y)`` is the number of rows in *Z*.

        ``C = contour(...)`` returns a
        :class:`~matplotlib.contour.QuadContourSet` object.

        Optional keyword arguments:

          *corner_mask*: bool, optional
            Enable/disable corner masking, which only has an effect if *Z* is
            a masked array.  If ``False``, any quad touching a masked point is
            masked out.  If ``True``, only the triangular corners of quads
            nearest those points are always masked out, other triangular
            corners comprising three unmasked points are contoured as usual.

            Defaults to ``rcParams['contour.corner_mask']``, which defaults to
            ``True``.

          *colors*: [ *None* | string | (mpl_colors) ]
            If *None*, the colormap specified by cmap will be used.

            If a string, like 'r' or 'red', all levels will be plotted in this
            color.

            If a tuple of matplotlib color args (string, float, rgb, etc),
            different levels will be plotted in different colors in the order
            specified.

          *alpha*: float
            The alpha blending value

          *cmap*: [ *None* | Colormap ]
            A cm :class:`~matplotlib.colors.Colormap` instance or
            *None*. If *cmap* is *None* and *colors* is *None*, a
            default Colormap is used.

          *norm*: [ *None* | Normalize ]
            A :class:`matplotlib.colors.Normalize` instance for
            scaling data values to colors. If *norm* is *None* and
            *colors* is *None*, the default linear scaling is used.

          *vmin*, *vmax*: [ *None* | scalar ]
            If not *None*, either or both of these values will be
            supplied to the :class:`matplotlib.colors.Normalize`
            instance, overriding the default color scaling based on
            *levels*.

          *levels*: [level0, level1, ..., leveln]
            A list of floating point numbers indicating the level
            curves to draw, in increasing order; e.g., to draw just
            the zero contour pass ``levels=[0]``

          *origin*: [ *None* | 'upper' | 'lower' | 'image' ]
            If *None*, the first value of *Z* will correspond to the
            lower left corner, location (0,0). If 'image', the rc
            value for ``image.origin`` will be used.

            This keyword is not active if *X* and *Y* are specified in
            the call to contour.

          *extent*: [ *None* | (x0,x1,y0,y1) ]

            If *origin* is not *None*, then *extent* is interpreted as
            in :func:`matplotlib.pyplot.imshow`: it gives the outer
            pixel boundaries. In this case, the position of Z[0,0]
            is the center of the pixel, not a corner. If *origin* is
            *None*, then (*x0*, *y0*) is the position of Z[0,0], and
            (*x1*, *y1*) is the position of Z[-1,-1].

            This keyword is not active if *X* and *Y* are specified in
            the call to contour.

          *locator*: [ *None* | ticker.Locator subclass ]
            If *locator* is *None*, the default
            :class:`~matplotlib.ticker.MaxNLocator` is used. The
            locator is used to determine the contour levels if they
            are not given explicitly via the *V* argument.

          *extend*: [ 'neither' | 'both' | 'min' | 'max' ]
            Unless this is 'neither', contour levels are automatically
            added to one or both ends of the range so that all data
            are included. These added ranges are then mapped to the
            special colormap values which default to the ends of the
            colormap range, but can be set via
            :meth:`matplotlib.colors.Colormap.set_under` and
            :meth:`matplotlib.colors.Colormap.set_over` methods.

          *xunits*, *yunits*: [ *None* | registered units ]
            Override axis units by specifying an instance of a
            :class:`matplotlib.units.ConversionInterface`.

          *antialiased*: bool
            enable antialiasing, overriding the defaults.  For
            filled contours, the default is *True*.  For line contours,
            it is taken from rcParams['lines.antialiased'].

          *nchunk*: [ 0 | integer ]
            If 0, no subdivision of the domain.  Specify a positive integer to
            divide the domain into subdomains of *nchunk* by *nchunk* quads.
            Chunking reduces the maximum length of polygons generated by the
            contouring algorithm which reduces the rendering workload passed
            on to the backend and also requires slightly less RAM.  It can
            however introduce rendering artifacts at chunk boundaries depending
            on the backend, the *antialiased* flag and value of *alpha*.

        contour-only keyword arguments:

          *linewidths*: [ *None* | number | tuple of numbers ]
            If *linewidths* is *None*, the default width in
            ``lines.linewidth`` in ``matplotlibrc`` is used.

            If a number, all levels will be plotted with this linewidth.

            If a tuple, different levels will be plotted with different
            linewidths in the order specified.

          *linestyles*: [ *None* | 'solid' | 'dashed' | 'dashdot' | 'dotted' ]
            If *linestyles* is *None*, the default is 'solid' unless
            the lines are monochrome.  In that case, negative
            contours will take their linestyle from the ``matplotlibrc``
            ``contour.negative_linestyle`` setting.

            *linestyles* can also be an iterable of the above strings
            specifying a set of linestyles to be used. If this
            iterable is shorter than the number of contour levels
            it will be repeated as necessary.

        contourf-only keyword arguments:

          *hatches*:
            A list of cross hatch patterns to use on the filled areas.
            If None, no hatching will be added to the contour.
            Hatching is supported in the PostScript, PDF, SVG and Agg
            backends only.


        Note: contourf fills intervals that are closed at the top; that
        is, for boundaries *z1* and *z2*, the filled region is::

            z1 < z <= z2

        There is one exception: if the lowest boundary coincides with
        the minimum value of the *z* array, then that minimum value
        will be included in the lowest interval.
        """
