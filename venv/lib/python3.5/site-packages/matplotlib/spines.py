from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib

from matplotlib.artist import allow_rasterization
from matplotlib import docstring
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import numpy as np
import warnings

rcParams = matplotlib.rcParams


class Spine(mpatches.Patch):
    """an axis spine -- the line noting the data area boundaries

    Spines are the lines connecting the axis tick marks and noting the
    boundaries of the data area. They can be placed at arbitrary
    positions. See function:`~matplotlib.spines.Spine.set_position`
    for more information.

    The default position is ``('outward',0)``.

    Spines are subclasses of class:`~matplotlib.patches.Patch`, and
    inherit much of their behavior.

    Spines draw a line, a circle, or an arc depending if
    function:`~matplotlib.spines.Spine.set_patch_line`,
    function:`~matplotlib.spines.Spine.set_patch_circle`, or
    function:`~matplotlib.spines.Spine.set_patch_arc` has been called.
    Line-like is the default.

    """
    def __str__(self):
        return "Spine"

    @docstring.dedent_interpd
    def __init__(self, axes, spine_type, path, **kwargs):
        """
        - *axes* : the Axes instance containing the spine
        - *spine_type* : a string specifying the spine type
        - *path* : the path instance used to draw the spine

        Valid kwargs are:
        %(Patch)s
        """
        super(Spine, self).__init__(**kwargs)
        self.axes = axes
        self.set_figure(self.axes.figure)
        self.spine_type = spine_type
        self.set_facecolor('none')
        self.set_edgecolor(rcParams['axes.edgecolor'])
        self.set_linewidth(rcParams['axes.linewidth'])
        self.set_capstyle('projecting')
        self.axis = None

        self.set_zorder(2.5)
        self.set_transform(self.axes.transData)  # default transform

        self._bounds = None  # default bounds
        self._smart_bounds = False

        # Defer initial position determination. (Not much support for
        # non-rectangular axes is currently implemented, and this lets
        # them pass through the spines machinery without errors.)
        self._position = None
        if not isinstance(path, matplotlib.path.Path):
            raise ValueError(
                "'path' must be an instance of 'matplotlib.path.Path'")
        self._path = path

        # To support drawing both linear and circular spines, this
        # class implements Patch behavior three ways. If
        # self._patch_type == 'line', behave like a mpatches.PathPatch
        # instance. If self._patch_type == 'circle', behave like a
        # mpatches.Ellipse instance. If self._patch_type == 'arc', behave like
        # a mpatches.Arc instance.
        self._patch_type = 'line'

        # Behavior copied from mpatches.Ellipse:
        # Note: This cannot be calculated until this is added to an Axes
        self._patch_transform = mtransforms.IdentityTransform()

    def set_smart_bounds(self, value):
        """set the spine and associated axis to have smart bounds"""
        self._smart_bounds = value

        # also set the axis if possible
        if self.spine_type in ('left', 'right'):
            self.axes.yaxis.set_smart_bounds(value)
        elif self.spine_type in ('top', 'bottom'):
            self.axes.xaxis.set_smart_bounds(value)
        self.stale = True

    def get_smart_bounds(self):
        """get whether the spine has smart bounds"""
        return self._smart_bounds

    def set_patch_arc(self, center, radius, theta1, theta2):
        """set the spine to be arc-like"""
        self._patch_type = 'arc'
        self._center = center
        self._width = radius * 2
        self._height = radius * 2
        self._theta1 = theta1
        self._theta2 = theta2
        self._path = mpath.Path.arc(theta1, theta2)
        # arc drawn on axes transform
        self.set_transform(self.axes.transAxes)
        self.stale = True

    def set_patch_circle(self, center, radius):
        """set the spine to be circular"""
        self._patch_type = 'circle'
        self._center = center
        self._width = radius * 2
        self._height = radius * 2
        # circle drawn on axes transform
        self.set_transform(self.axes.transAxes)
        self.stale = True

    def set_patch_line(self):
        """set the spine to be linear"""
        self._patch_type = 'line'
        self.stale = True

    # Behavior copied from mpatches.Ellipse:
    def _recompute_transform(self):
        """NOTE: This cannot be called until after this has been added
                 to an Axes, otherwise unit conversion will fail. This
                 makes it very important to call the accessor method and
                 not directly access the transformation member variable.
        """
        assert self._patch_type in ('arc', 'circle')
        center = (self.convert_xunits(self._center[0]),
                  self.convert_yunits(self._center[1]))
        width = self.convert_xunits(self._width)
        height = self.convert_yunits(self._height)
        self._patch_transform = mtransforms.Affine2D() \
            .scale(width * 0.5, height * 0.5) \
            .translate(*center)

    def get_patch_transform(self):
        if self._patch_type in ('arc', 'circle'):
            self._recompute_transform()
            return self._patch_transform
        else:
            return super(Spine, self).get_patch_transform()

    def get_path(self):
        return self._path

    def _ensure_position_is_set(self):
        if self._position is None:
            # default position
            self._position = ('outward', 0.0)  # in points
            self.set_position(self._position)

    def register_axis(self, axis):
        """register an axis

        An axis should be registered with its corresponding spine from
        the Axes instance. This allows the spine to clear any axis
        properties when needed.
        """
        self.axis = axis
        if self.axis is not None:
            self.axis.cla()
        self.stale = True

    def cla(self):
        """Clear the current spine"""
        self._position = None  # clear position
        if self.axis is not None:
            self.axis.cla()

    def is_frame_like(self):
        """return True if directly on axes frame

        This is useful for determining if a spine is the edge of an
        old style MPL plot. If so, this function will return True.
        """
        self._ensure_position_is_set()
        position = self._position
        if isinstance(position, six.string_types):
            if position == 'center':
                position = ('axes', 0.5)
            elif position == 'zero':
                position = ('data', 0)
        if len(position) != 2:
            raise ValueError("position should be 2-tuple")
        position_type, amount = position
        if position_type == 'outward' and amount == 0:
            return True
        else:
            return False

    def _adjust_location(self):
        """automatically set spine bounds to the view interval"""

        if self.spine_type == 'circle':
            return

        if self._bounds is None:
            if self.spine_type in ('left', 'right'):
                low, high = self.axes.viewLim.intervaly
            elif self.spine_type in ('top', 'bottom'):
                low, high = self.axes.viewLim.intervalx
            else:
                raise ValueError('unknown spine spine_type: %s' %
                                 self.spine_type)

            if self._smart_bounds:
                # attempt to set bounds in sophisticated way

                # handle inverted limits
                viewlim_low, viewlim_high = sorted([low, high])

                if self.spine_type in ('left', 'right'):
                    datalim_low, datalim_high = self.axes.dataLim.intervaly
                    ticks = self.axes.get_yticks()
                elif self.spine_type in ('top', 'bottom'):
                    datalim_low, datalim_high = self.axes.dataLim.intervalx
                    ticks = self.axes.get_xticks()
                # handle inverted limits
                ticks = np.sort(ticks)
                datalim_low, datalim_high = sorted([datalim_low, datalim_high])

                if datalim_low < viewlim_low:
                    # Data extends past view. Clip line to view.
                    low = viewlim_low
                else:
                    # Data ends before view ends.
                    cond = (ticks <= datalim_low) & (ticks >= viewlim_low)
                    tickvals = ticks[cond]
                    if len(tickvals):
                        # A tick is less than or equal to lowest data point.
                        low = tickvals[-1]
                    else:
                        # No tick is available
                        low = datalim_low
                    low = max(low, viewlim_low)

                if datalim_high > viewlim_high:
                    # Data extends past view. Clip line to view.
                    high = viewlim_high
                else:
                    # Data ends before view ends.
                    cond = (ticks >= datalim_high) & (ticks <= viewlim_high)
                    tickvals = ticks[cond]
                    if len(tickvals):
                        # A tick is greater than or equal to highest data
                        # point.
                        high = tickvals[0]
                    else:
                        # No tick is available
                        high = datalim_high
                    high = min(high, viewlim_high)

        else:
            low, high = self._bounds

        if self._patch_type == 'arc':
            if self.spine_type in ('bottom', 'top'):
                try:
                    direction = self.axes.get_theta_direction()
                except AttributeError:
                    direction = 1
                try:
                    offset = self.axes.get_theta_offset()
                except AttributeError:
                    offset = 0
                low = low * direction + offset
                high = high * direction + offset
                if low > high:
                    low, high = high, low

                self._path = mpath.Path.arc(np.rad2deg(low), np.rad2deg(high))

                if self.spine_type == 'bottom':
                    rmin, rmax = self.axes.viewLim.intervaly
                    try:
                        rorigin = self.axes.get_rorigin()
                    except AttributeError:
                        rorigin = rmin
                    scaled_diameter = (rmin - rorigin) / (rmax - rorigin)
                    self._height = scaled_diameter
                    self._width = scaled_diameter

            else:
                raise ValueError('unable to set bounds for spine "%s"' %
                                 self.spine_type)
        else:
            v1 = self._path.vertices
            assert v1.shape == (2, 2), 'unexpected vertices shape'
            if self.spine_type in ['left', 'right']:
                v1[0, 1] = low
                v1[1, 1] = high
            elif self.spine_type in ['bottom', 'top']:
                v1[0, 0] = low
                v1[1, 0] = high
            else:
                raise ValueError('unable to set bounds for spine "%s"' %
                                 self.spine_type)

    @allow_rasterization
    def draw(self, renderer):
        self._adjust_location()
        ret = super(Spine, self).draw(renderer)
        self.stale = False
        return ret

    def _calc_offset_transform(self):
        """calculate the offset transform performed by the spine"""
        self._ensure_position_is_set()
        position = self._position
        if isinstance(position, six.string_types):
            if position == 'center':
                position = ('axes', 0.5)
            elif position == 'zero':
                position = ('data', 0)
        assert len(position) == 2, "position should be 2-tuple"
        position_type, amount = position
        assert position_type in ('axes', 'outward', 'data')
        if position_type == 'outward':
            if amount == 0:
                # short circuit commonest case
                self._spine_transform = ('identity',
                                         mtransforms.IdentityTransform())
            elif self.spine_type in ['left', 'right', 'top', 'bottom']:
                offset_vec = {'left': (-1, 0),
                              'right': (1, 0),
                              'bottom': (0, -1),
                              'top': (0, 1),
                              }[self.spine_type]
                # calculate x and y offset in dots
                offset_x = amount * offset_vec[0] / 72.0
                offset_y = amount * offset_vec[1] / 72.0
                self._spine_transform = ('post',
                                         mtransforms.ScaledTranslation(
                                             offset_x,
                                             offset_y,
                                             self.figure.dpi_scale_trans))
            else:
                warnings.warn('unknown spine type "%s": no spine '
                              'offset performed' % self.spine_type)
                self._spine_transform = ('identity',
                                         mtransforms.IdentityTransform())
        elif position_type == 'axes':
            if self.spine_type in ('left', 'right'):
                self._spine_transform = ('pre',
                                         mtransforms.Affine2D.from_values(
                                             # keep y unchanged, fix x at
                                             # amount
                                             0, 0, 0, 1, amount, 0))
            elif self.spine_type in ('bottom', 'top'):
                self._spine_transform = ('pre',
                                         mtransforms.Affine2D.from_values(
                                             # keep x unchanged, fix y at
                                             # amount
                                             1, 0, 0, 0, 0, amount))
            else:
                warnings.warn('unknown spine type "%s": no spine '
                              'offset performed' % self.spine_type)
                self._spine_transform = ('identity',
                                         mtransforms.IdentityTransform())
        elif position_type == 'data':
            if self.spine_type in ('right', 'top'):
                # The right and top spines have a default position of 1 in
                # axes coordinates.  When specifying the position in data
                # coordinates, we need to calculate the position relative to 0.
                amount -= 1
            if self.spine_type in ('left', 'right'):
                self._spine_transform = ('data',
                                         mtransforms.Affine2D().translate(
                                             amount, 0))
            elif self.spine_type in ('bottom', 'top'):
                self._spine_transform = ('data',
                                         mtransforms.Affine2D().translate(
                                             0, amount))
            else:
                warnings.warn('unknown spine type "%s": no spine '
                              'offset performed' % self.spine_type)
                self._spine_transform = ('identity',
                                         mtransforms.IdentityTransform())

    def set_position(self, position):
        """set the position of the spine

        Spine position is specified by a 2 tuple of (position type,
        amount). The position types are:

        * 'outward' : place the spine out from the data area by the
          specified number of points. (Negative values specify placing the
          spine inward.)

        * 'axes' : place the spine at the specified Axes coordinate (from
          0.0-1.0).

        * 'data' : place the spine at the specified data coordinate.

        Additionally, shorthand notations define a special positions:

        * 'center' -> ('axes',0.5)
        * 'zero' -> ('data', 0.0)

        """
        if position in ('center', 'zero'):
            # special positions
            pass
        else:
            if len(position) != 2:
                raise ValueError("position should be 'center' or 2-tuple")
            if position[0] not in ['outward', 'axes', 'data']:
                raise ValueError("position[0] should be one of 'outward', "
                                 "'axes', or 'data' ")
        self._position = position
        self._calc_offset_transform()

        self.set_transform(self.get_spine_transform())

        if self.axis is not None:
            self.axis.reset_ticks()
        self.stale = True

    def get_position(self):
        """get the spine position"""
        self._ensure_position_is_set()
        return self._position

    def get_spine_transform(self):
        """get the spine transform"""
        self._ensure_position_is_set()
        what, how = self._spine_transform

        if what == 'data':
            # special case data based spine locations
            data_xform = self.axes.transScale + \
                (how + self.axes.transLimits + self.axes.transAxes)
            if self.spine_type in ['left', 'right']:
                result = mtransforms.blended_transform_factory(
                    data_xform, self.axes.transData)
            elif self.spine_type in ['top', 'bottom']:
                result = mtransforms.blended_transform_factory(
                    self.axes.transData, data_xform)
            else:
                raise ValueError('unknown spine spine_type: %s' %
                                 self.spine_type)
            return result

        if self.spine_type in ['left', 'right']:
            base_transform = self.axes.get_yaxis_transform(which='grid')
        elif self.spine_type in ['top', 'bottom']:
            base_transform = self.axes.get_xaxis_transform(which='grid')
        else:
            raise ValueError('unknown spine spine_type: %s' %
                             self.spine_type)

        if what == 'identity':
            return base_transform
        elif what == 'post':
            return base_transform + how
        elif what == 'pre':
            return how + base_transform
        else:
            raise ValueError("unknown spine_transform type: %s" % what)

    def set_bounds(self, low, high):
        """Set the bounds of the spine."""
        if self.spine_type == 'circle':
            raise ValueError(
                'set_bounds() method incompatible with circular spines')
        self._bounds = (low, high)
        self.stale = True

    def get_bounds(self):
        """Get the bounds of the spine."""
        return self._bounds

    @classmethod
    def linear_spine(cls, axes, spine_type, **kwargs):
        """
        (staticmethod) Returns a linear :class:`Spine`.
        """
        # all values of 13 get replaced upon call to set_bounds()
        if spine_type == 'left':
            path = mpath.Path([(0.0, 13), (0.0, 13)])
        elif spine_type == 'right':
            path = mpath.Path([(1.0, 13), (1.0, 13)])
        elif spine_type == 'bottom':
            path = mpath.Path([(13, 0.0), (13, 0.0)])
        elif spine_type == 'top':
            path = mpath.Path([(13, 1.0), (13, 1.0)])
        else:
            raise ValueError('unable to make path for spine "%s"' % spine_type)
        result = cls(axes, spine_type, path, **kwargs)
        result.set_visible(rcParams['axes.spines.{0}'.format(spine_type)])

        return result

    @classmethod
    def arc_spine(cls, axes, spine_type, center, radius, theta1, theta2,
                  **kwargs):
        """
        (classmethod) Returns an arc :class:`Spine`.
        """
        path = mpath.Path.arc(theta1, theta2)
        result = cls(axes, spine_type, path, **kwargs)
        result.set_patch_arc(center, radius, theta1, theta2)
        return result

    @classmethod
    def circular_spine(cls, axes, center, radius, **kwargs):
        """
        (staticmethod) Returns a circular :class:`Spine`.
        """
        path = mpath.Path.unit_circle()
        spine_type = 'circle'
        result = cls(axes, spine_type, path, **kwargs)
        result.set_patch_circle(center, radius)
        return result

    def set_color(self, c):
        """
        Set the edgecolor.

        ACCEPTS: matplotlib color arg or sequence of rgba tuples

        .. seealso::

            :meth:`set_facecolor`, :meth:`set_edgecolor`
               For setting the edge or face color individually.
        """
        # The facecolor of a spine is always 'none' by default -- let
        # the user change it manually if desired.
        self.set_edgecolor(c)
        self.stale = True
