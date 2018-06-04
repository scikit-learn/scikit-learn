"""
axes3d.py, original mplot3d version by John Porter
Created: 23 Sep 2005

Parts fixed by Reinier Heeres <reinier@heeres.eu>
Minor additions by Ben Axelrod <baxelrod@coroware.com>
Significant updates and revisions by Ben Root <ben.v.root@gmail.com>

Module containing Axes3D, an object which can plot 3D objects on a
2D matplotlib figure.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import map, xrange, zip, reduce

import math
import warnings
from collections import defaultdict

import numpy as np

import matplotlib.axes as maxes
import matplotlib.cbook as cbook
import matplotlib.collections as mcoll
import matplotlib.colors as mcolors
import matplotlib.docstring as docstring
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
from matplotlib.axes import Axes, rcParams
from matplotlib.cbook import _backports
from matplotlib.colors import Normalize, LightSource
from matplotlib.transforms import Bbox
from matplotlib.tri.triangulation import Triangulation

from . import art3d
from . import proj3d
from . import axis3d


def unit_bbox():
    box = Bbox(np.array([[0, 0], [1, 1]]))
    return box


class Axes3D(Axes):
    """
    3D axes object.
    """
    name = '3d'
    _shared_z_axes = cbook.Grouper()

    def __init__(self, fig, rect=None, *args, **kwargs):
        '''
        Build an :class:`Axes3D` instance in
        :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* in
        :class:`~matplotlib.figure.Figure` coordinates

        Optional keyword arguments:

          ================   =========================================
          Keyword            Description
          ================   =========================================
          *azim*             Azimuthal viewing angle (default -60)
          *elev*             Elevation viewing angle (default 30)
          *zscale*           [%(scale)s]
          *sharez*           Other axes to share z-limits with
          *proj_type*        'persp' or 'ortho' (default 'persp')
          ================   =========================================

        .. versionadded :: 1.2.1
            *sharez*

        ''' % {'scale': ' | '.join([repr(x) for x in mscale.get_scale_names()])}

        if rect is None:
            rect = [0.0, 0.0, 1.0, 1.0]
        self._cids = []

        self.initial_azim = kwargs.pop('azim', -60)
        self.initial_elev = kwargs.pop('elev', 30)
        zscale = kwargs.pop('zscale', None)
        sharez = kwargs.pop('sharez', None)
        self.set_proj_type(kwargs.pop('proj_type', 'persp'))

        self.xy_viewLim = unit_bbox()
        self.zz_viewLim = unit_bbox()
        self.xy_dataLim = unit_bbox()
        self.zz_dataLim = unit_bbox()
        # inihibit autoscale_view until the axes are defined
        # they can't be defined until Axes.__init__ has been called
        self.view_init(self.initial_elev, self.initial_azim)
        self._ready = 0

        self._sharez = sharez
        if sharez is not None:
            self._shared_z_axes.join(self, sharez)
            self._adjustable = 'datalim'

        super(Axes3D, self).__init__(fig, rect,
                                     frameon=True,
                                     *args, **kwargs)
        # Disable drawing of axes by base class
        super(Axes3D, self).set_axis_off()
        # Enable drawing of axes by Axes3D class
        self.set_axis_on()
        self.M = None

        # func used to format z -- fall back on major formatters
        self.fmt_zdata = None

        if zscale is not None:
            self.set_zscale(zscale)

        if self.zaxis is not None:
            self._zcid = self.zaxis.callbacks.connect(
                'units finalize', lambda: self._on_units_changed(scalez=True))
        else:
            self._zcid = None

        self._ready = 1
        self.mouse_init()
        self.set_top_view()

        self.patch.set_linewidth(0)
        # Calculate the pseudo-data width and height
        pseudo_bbox = self.transLimits.inverted().transform([(0, 0), (1, 1)])
        self._pseudo_w, self._pseudo_h = pseudo_bbox[1] - pseudo_bbox[0]

        self.figure.add_axes(self)

    def set_axis_off(self):
        self._axis3don = False
        self.stale = True

    def set_axis_on(self):
        self._axis3don = True
        self.stale = True

    def have_units(self):
        """
        Return *True* if units are set on the *x*, *y*, or *z* axes

        """
        return (self.xaxis.have_units() or self.yaxis.have_units() or
                self.zaxis.have_units())

    def convert_zunits(self, z):
        """
        For artists in an axes, if the zaxis has units support,
        convert *z* using zaxis unit type

        .. versionadded :: 1.2.1

        """
        return self.zaxis.convert_units(z)

    def _process_unit_info(self, xdata=None, ydata=None, zdata=None,
                           kwargs=None):
        """
        Look for unit *kwargs* and update the axis instances as necessary

        """
        super(Axes3D, self)._process_unit_info(xdata=xdata, ydata=ydata,
                                               kwargs=kwargs)

        if self.xaxis is None or self.yaxis is None or self.zaxis is None:
            return

        if zdata is not None:
            # we only need to update if there is nothing set yet.
            if not self.zaxis.have_units():
                self.zaxis.update_units(xdata)

        # process kwargs 2nd since these will override default units
        if kwargs is not None:
            zunits = kwargs.pop('zunits', self.zaxis.units)
            if zunits != self.zaxis.units:
                self.zaxis.set_units(zunits)
                # If the units being set imply a different converter,
                # we need to update.
                if zdata is not None:
                    self.zaxis.update_units(zdata)

    def set_top_view(self):
        # this happens to be the right view for the viewing coordinates
        # moved up and to the left slightly to fit labels and axes
        xdwl = (0.95/self.dist)
        xdw = (0.9/self.dist)
        ydwl = (0.95/self.dist)
        ydw = (0.9/self.dist)

        # This is purposely using the 2D Axes's set_xlim and set_ylim,
        # because we are trying to place our viewing pane.
        super(Axes3D, self).set_xlim(-xdwl, xdw, auto=None)
        super(Axes3D, self).set_ylim(-ydwl, ydw, auto=None)

    def _init_axis(self):
        '''Init 3D axes; overrides creation of regular X/Y axes'''
        self.w_xaxis = axis3d.XAxis('x', self.xy_viewLim.intervalx,
                                    self.xy_dataLim.intervalx, self)
        self.xaxis = self.w_xaxis
        self.w_yaxis = axis3d.YAxis('y', self.xy_viewLim.intervaly,
                                    self.xy_dataLim.intervaly, self)
        self.yaxis = self.w_yaxis
        self.w_zaxis = axis3d.ZAxis('z', self.zz_viewLim.intervalx,
                                    self.zz_dataLim.intervalx, self)
        self.zaxis = self.w_zaxis

        for ax in self.xaxis, self.yaxis, self.zaxis:
            ax.init3d()

    def get_children(self):
        return [self.zaxis, ] + super(Axes3D, self).get_children()

    def _get_axis_list(self):
        return super(Axes3D, self)._get_axis_list() + (self.zaxis, )

    def unit_cube(self, vals=None):
        minx, maxx, miny, maxy, minz, maxz = vals or self.get_w_lims()
        xs, ys, zs = ([minx, maxx, maxx, minx, minx, maxx, maxx, minx],
                      [miny, miny, maxy, maxy, miny, miny, maxy, maxy],
                      [minz, minz, minz, minz, maxz, maxz, maxz, maxz])
        return list(zip(xs, ys, zs))

    def tunit_cube(self, vals=None, M=None):
        if M is None:
            M = self.M
        xyzs = self.unit_cube(vals)
        tcube = proj3d.proj_points(xyzs, M)
        return tcube

    def tunit_edges(self, vals=None, M=None):
        tc = self.tunit_cube(vals, M)
        edges = [(tc[0], tc[1]),
                 (tc[1], tc[2]),
                 (tc[2], tc[3]),
                 (tc[3], tc[0]),

                 (tc[0], tc[4]),
                 (tc[1], tc[5]),
                 (tc[2], tc[6]),
                 (tc[3], tc[7]),

                 (tc[4], tc[5]),
                 (tc[5], tc[6]),
                 (tc[6], tc[7]),
                 (tc[7], tc[4])]
        return edges

    def draw(self, renderer):
        # draw the background patch
        self.patch.draw(renderer)
        self._frameon = False

        # first, set the aspect
        # this is duplicated from `axes._base._AxesBase.draw`
        # but must be called before any of the artist are drawn as
        # it adjusts the view limits and the size of the bounding box
        # of the axes
        locator = self.get_axes_locator()
        if locator:
            pos = locator(self, renderer)
            self.apply_aspect(pos)
        else:
            self.apply_aspect()

        # add the projection matrix to the renderer
        self.M = self.get_proj()
        renderer.M = self.M
        renderer.vvec = self.vvec
        renderer.eye = self.eye
        renderer.get_axis_position = self.get_axis_position

        # Calculate projection of collections and zorder them
        for i, col in enumerate(
                sorted(self.collections,
                       key=lambda col: col.do_3d_projection(renderer),
                       reverse=True)):
            col.zorder = i

        # Calculate projection of patches and zorder them
        for i, patch in enumerate(
                sorted(self.patches,
                       key=lambda patch: patch.do_3d_projection(renderer),
                       reverse=True)):
            patch.zorder = i

        if self._axis3don:
            axes = (self.xaxis, self.yaxis, self.zaxis)
            # Draw panes first
            for ax in axes:
                ax.draw_pane(renderer)
            # Then axes
            for ax in axes:
                ax.draw(renderer)

        # Then rest
        super(Axes3D, self).draw(renderer)

    def get_axis_position(self):
        vals = self.get_w_lims()
        tc = self.tunit_cube(vals, self.M)
        xhigh = tc[1][2] > tc[2][2]
        yhigh = tc[3][2] > tc[2][2]
        zhigh = tc[0][2] > tc[2][2]
        return xhigh, yhigh, zhigh

    def _on_units_changed(self, scalex=False, scaley=False, scalez=False):
        """
        Callback for processing changes to axis units.

        Currently forces updates of data limits and view limits.
        """
        self.relim()
        self.autoscale_view(scalex=scalex, scaley=scaley, scalez=scalez)

    def update_datalim(self, xys, **kwargs):
        pass

    def get_autoscale_on(self):
        """
        Get whether autoscaling is applied for all axes on plot commands

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        return super(Axes3D, self).get_autoscale_on() and self.get_autoscalez_on()

    def get_autoscalez_on(self):
        """
        Get whether autoscaling for the z-axis is applied on plot commands

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        return self._autoscaleZon

    def set_autoscale_on(self, b):
        """
        Set whether autoscaling is applied on plot commands

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        super(Axes3D, self).set_autoscale_on(b)
        self.set_autoscalez_on(b)

    def set_autoscalez_on(self, b):
        """
        Set whether autoscaling for the z-axis is applied on plot commands

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._autoscaleZon = b

    def set_zmargin(self, m):
        """
        Set padding of Z data limits prior to autoscaling.

        *m* times the data interval will be added to each
        end of that interval before it is used in autoscaling.

        accepts: float in range 0 to 1

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        if m < 0 or m > 1 :
            raise ValueError("margin must be in range 0 to 1")
        self._zmargin = m
        self.stale = True

    def margins(self, *args, **kw):
        """
        Convenience method to set or retrieve autoscaling margins.

        signatures::
            margins()

        returns xmargin, ymargin, zmargin

        ::

            margins(margin)

            margins(xmargin, ymargin, zmargin)

            margins(x=xmargin, y=ymargin, z=zmargin)

            margins(..., tight=False)

        All forms above set the xmargin, ymargin and zmargin
        parameters. All keyword parameters are optional.  A single argument
        specifies xmargin, ymargin and zmargin.  The *tight* parameter
        is passed to :meth:`autoscale_view`, which is executed after
        a margin is changed; the default here is *True*, on the
        assumption that when margins are specified, no additional
        padding to match tick marks is usually desired.  Setting
        *tight* to *None* will preserve the previous setting.

        Specifying any margin changes only the autoscaling; for example,
        if *xmargin* is not None, then *xmargin* times the X data
        interval will be added to each end of that interval before
        it is used in autoscaling.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        if not args and not kw:
            return self._xmargin, self._ymargin, self._zmargin

        tight = kw.pop('tight', True)
        mx = kw.pop('x', None)
        my = kw.pop('y', None)
        mz = kw.pop('z', None)
        if not args:
            pass
        elif len(args) == 1:
            mx = my = mz = args[0]
        elif len(args) == 2:
            warnings.warn(
                "Passing exactly two positional arguments to Axes3D.margins "
                "is deprecated.  If needed, pass them as keyword arguments "
                "instead", cbook.mplDeprecation)
            mx, my = args
        elif len(args) == 3:
            mx, my, mz = args
        else:
            raise ValueError(
                "Axes3D.margins takes at most three positional arguments")
        if mx is not None:
            self.set_xmargin(mx)
        if my is not None:
            self.set_ymargin(my)
        if mz is not None:
            self.set_zmargin(mz)

        scalex = mx is not None
        scaley = my is not None
        scalez = mz is not None

        self.autoscale_view(tight=tight, scalex=scalex, scaley=scaley,
                                         scalez=scalez)

    def autoscale(self, enable=True, axis='both', tight=None):
        """
        Convenience method for simple axis view autoscaling.
        See :meth:`matplotlib.axes.Axes.autoscale` for full explanation.
        Note that this function behaves the same, but for all
        three axes.  Therefore, 'z' can be passed for *axis*,
        and 'both' applies to all three axes.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        if enable is None:
            scalex = True
            scaley = True
            scalez = True
        else:
            if axis in ['x', 'both']:
                self._autoscaleXon = scalex = bool(enable)
            else:
                scalex = False
            if axis in ['y', 'both']:
                self._autoscaleYon = scaley = bool(enable)
            else:
                scaley = False
            if axis in ['z', 'both']:
                self._autoscaleZon = scalez = bool(enable)
            else:
                scalez = False
        self.autoscale_view(tight=tight, scalex=scalex, scaley=scaley,
                                         scalez=scalez)

    def auto_scale_xyz(self, X, Y, Z=None, had_data=None):
        x, y, z = map(np.asarray, (X, Y, Z))
        try:
            x, y = x.flatten(), y.flatten()
            if Z is not None:
                z = z.flatten()
        except AttributeError:
            raise

        # This updates the bounding boxes as to keep a record as
        # to what the minimum sized rectangular volume holds the
        # data.
        self.xy_dataLim.update_from_data_xy(np.array([x, y]).T, not had_data)
        if z is not None:
            self.zz_dataLim.update_from_data_xy(np.array([z, z]).T, not had_data)

        # Let autoscale_view figure out how to use this data.
        self.autoscale_view()

    def autoscale_view(self, tight=None, scalex=True, scaley=True,
                             scalez=True):
        """
        Autoscale the view limits using the data limits.
        See :meth:`matplotlib.axes.Axes.autoscale_view` for documentation.
        Note that this function applies to the 3D axes, and as such
        adds the *scalez* to the function arguments.

        .. versionchanged :: 1.1.0
            Function signature was changed to better match the 2D version.
            *tight* is now explicitly a kwarg and placed first.

        .. versionchanged :: 1.2.1
            This is now fully functional.

        """
        if not self._ready:
            return

        # This method looks at the rectangular volume (see above)
        # of data and decides how to scale the view portal to fit it.
        if tight is None:
            # if image data only just use the datalim
            _tight = self._tight or (len(self.images)>0 and
                                     len(self.lines)==0 and
                                     len(self.patches)==0)
        else:
            _tight = self._tight = bool(tight)

        if scalex and self._autoscaleXon:
            xshared = self._shared_x_axes.get_siblings(self)
            dl = [ax.dataLim for ax in xshared]
            bb = mtransforms.BboxBase.union(dl)
            x0, x1 = self.xy_dataLim.intervalx
            xlocator = self.xaxis.get_major_locator()
            try:
                x0, x1 = xlocator.nonsingular(x0, x1)
            except AttributeError:
                x0, x1 = mtransforms.nonsingular(x0, x1, increasing=False,
                                                         expander=0.05)
            if self._xmargin > 0:
                delta = (x1 - x0) * self._xmargin
                x0 -= delta
                x1 += delta
            if not _tight:
                x0, x1 = xlocator.view_limits(x0, x1)
            self.set_xbound(x0, x1)

        if scaley and self._autoscaleYon:
            yshared = self._shared_y_axes.get_siblings(self)
            dl = [ax.dataLim for ax in yshared]
            bb = mtransforms.BboxBase.union(dl)
            y0, y1 = self.xy_dataLim.intervaly
            ylocator = self.yaxis.get_major_locator()
            try:
                y0, y1 = ylocator.nonsingular(y0, y1)
            except AttributeError:
                y0, y1 = mtransforms.nonsingular(y0, y1, increasing=False,
                                                         expander=0.05)
            if self._ymargin > 0:
                delta = (y1 - y0) * self._ymargin
                y0 -= delta
                y1 += delta
            if not _tight:
                y0, y1 = ylocator.view_limits(y0, y1)
            self.set_ybound(y0, y1)

        if scalez and self._autoscaleZon:
            zshared = self._shared_z_axes.get_siblings(self)
            dl = [ax.dataLim for ax in zshared]
            bb = mtransforms.BboxBase.union(dl)
            z0, z1 = self.zz_dataLim.intervalx
            zlocator = self.zaxis.get_major_locator()
            try:
                z0, z1 = zlocator.nonsingular(z0, z1)
            except AttributeError:
                z0, z1 = mtransforms.nonsingular(z0, z1, increasing=False,
                                                         expander=0.05)
            if self._zmargin > 0:
                delta = (z1 - z0) * self._zmargin
                z0 -= delta
                z1 += delta
            if not _tight:
                z0, z1 = zlocator.view_limits(z0, z1)
            self.set_zbound(z0, z1)

    def get_w_lims(self):
        '''Get 3D world limits.'''
        minx, maxx = self.get_xlim3d()
        miny, maxy = self.get_ylim3d()
        minz, maxz = self.get_zlim3d()
        return minx, maxx, miny, maxy, minz, maxz

    def _determine_lims(self, xmin=None, xmax=None, *args, **kwargs):
        if xmax is None and cbook.iterable(xmin):
            xmin, xmax = xmin
        if xmin == xmax:
            xmin -= 0.05
            xmax += 0.05
        return (xmin, xmax)

    def set_xlim3d(self, left=None, right=None, emit=True, auto=False, **kw):
        """
        Set 3D x limits.

        See :meth:`matplotlib.axes.Axes.set_xlim` for full documentation.

        """
        if 'xmin' in kw:
            left = kw.pop('xmin')
        if 'xmax' in kw:
            right = kw.pop('xmax')
        if kw:
            raise ValueError("unrecognized kwargs: %s" % list(kw))

        if right is None and cbook.iterable(left):
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
            warnings.warn(('Attempting to set identical left==right results\n'
                     'in singular transformations; automatically expanding.\n'
                     'left=%s, right=%s') % (left, right))
        left, right = mtransforms.nonsingular(left, right, increasing=False)
        left, right = self.xaxis.limit_range_for_scale(left, right)
        self.xy_viewLim.intervalx = (left, right)

        if auto is not None:
            self._autoscaleXon = bool(auto)

        if emit:
            self.callbacks.process('xlim_changed', self)
            # Call all of the other x-axes that are shared with this one
            for other in self._shared_x_axes.get_siblings(self):
                if other is not self:
                    other.set_xlim(self.xy_viewLim.intervalx,
                                            emit=False, auto=auto)
                    if (other.figure != self.figure and
                        other.figure.canvas is not None):
                        other.figure.canvas.draw_idle()
        self.stale = True
        return left, right
    set_xlim = set_xlim3d

    def set_ylim3d(self, bottom=None, top=None, emit=True, auto=False, **kw):
        """
        Set 3D y limits.

        See :meth:`matplotlib.axes.Axes.set_ylim` for full documentation.

        """
        if 'ymin' in kw:
            bottom = kw.pop('ymin')
        if 'ymax' in kw:
            top = kw.pop('ymax')
        if kw:
            raise ValueError("unrecognized kwargs: %s" % list(kw))

        if top is None and cbook.iterable(bottom):
            bottom, top = bottom

        self._process_unit_info(ydata=(bottom, top))
        bottom = self._validate_converted_limits(bottom, self.convert_yunits)
        top = self._validate_converted_limits(top, self.convert_yunits)

        old_bottom, old_top = self.get_ylim()
        if bottom is None:
            bottom = old_bottom
        if top is None:
            top = old_top

        if top == bottom:
            warnings.warn(('Attempting to set identical bottom==top results\n'
                     'in singular transformations; automatically expanding.\n'
                     'bottom=%s, top=%s') % (bottom, top))
        bottom, top = mtransforms.nonsingular(bottom, top, increasing=False)
        bottom, top = self.yaxis.limit_range_for_scale(bottom, top)
        self.xy_viewLim.intervaly = (bottom, top)

        if auto is not None:
            self._autoscaleYon = bool(auto)

        if emit:
            self.callbacks.process('ylim_changed', self)
            # Call all of the other y-axes that are shared with this one
            for other in self._shared_y_axes.get_siblings(self):
                if other is not self:
                    other.set_ylim(self.xy_viewLim.intervaly,
                                            emit=False, auto=auto)
                    if (other.figure != self.figure and
                        other.figure.canvas is not None):
                        other.figure.canvas.draw_idle()
        self.stale = True
        return bottom, top
    set_ylim = set_ylim3d

    def set_zlim3d(self, bottom=None, top=None, emit=True, auto=False, **kw):
        """
        Set 3D z limits.

        See :meth:`matplotlib.axes.Axes.set_ylim` for full documentation

        """
        if 'zmin' in kw:
            bottom = kw.pop('zmin')
        if 'zmax' in kw:
            top = kw.pop('zmax')
        if kw:
            raise ValueError("unrecognized kwargs: %s" % list(kw))

        if top is None and cbook.iterable(bottom):
            bottom, top = bottom

        self._process_unit_info(zdata=(bottom, top))
        bottom = self._validate_converted_limits(bottom, self.convert_zunits)
        top = self._validate_converted_limits(top, self.convert_zunits)

        old_bottom, old_top = self.get_zlim()
        if bottom is None:
            bottom = old_bottom
        if top is None:
            top = old_top

        if top == bottom:
            warnings.warn(('Attempting to set identical bottom==top results\n'
                     'in singular transformations; automatically expanding.\n'
                     'bottom=%s, top=%s') % (bottom, top))
        bottom, top = mtransforms.nonsingular(bottom, top, increasing=False)
        bottom, top = self.zaxis.limit_range_for_scale(bottom, top)
        self.zz_viewLim.intervalx = (bottom, top)

        if auto is not None:
            self._autoscaleZon = bool(auto)

        if emit:
            self.callbacks.process('zlim_changed', self)
            # Call all of the other y-axes that are shared with this one
            for other in self._shared_z_axes.get_siblings(self):
                if other is not self:
                    other.set_zlim(self.zz_viewLim.intervalx,
                                            emit=False, auto=auto)
                    if (other.figure != self.figure and
                        other.figure.canvas is not None):
                        other.figure.canvas.draw_idle()
        self.stale = True
        return bottom, top
    set_zlim = set_zlim3d

    def get_xlim3d(self):
        return tuple(self.xy_viewLim.intervalx)
    get_xlim3d.__doc__ = maxes.Axes.get_xlim.__doc__
    get_xlim = get_xlim3d
    if get_xlim.__doc__ is not None:
        get_xlim.__doc__ += """
            .. versionchanged :: 1.1.0
                This function now correctly refers to the 3D x-limits
            """

    def get_ylim3d(self):
        return tuple(self.xy_viewLim.intervaly)
    get_ylim3d.__doc__ = maxes.Axes.get_ylim.__doc__
    get_ylim = get_ylim3d
    if get_ylim.__doc__ is not None:
        get_ylim.__doc__ += """
            .. versionchanged :: 1.1.0
                This function now correctly refers to the 3D y-limits.
            """

    def get_zlim3d(self):
        '''Get 3D z limits.'''
        return tuple(self.zz_viewLim.intervalx)
    get_zlim = get_zlim3d

    def get_zscale(self):
        """
        Return the zaxis scale string %s

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """ % (", ".join(mscale.get_scale_names()))
        return self.zaxis.get_scale()

    # We need to slightly redefine these to pass scalez=False
    # to their calls of autoscale_view.
    def set_xscale(self, value, **kwargs):
        self.xaxis._set_scale(value, **kwargs)
        self.autoscale_view(scaley=False, scalez=False)
        self._update_transScale()
    if maxes.Axes.set_xscale.__doc__ is not None:
        set_xscale.__doc__ = maxes.Axes.set_xscale.__doc__ + """

            .. versionadded :: 1.1.0
                This function was added, but not tested. Please report any bugs.
            """

    def set_yscale(self, value, **kwargs):
        self.yaxis._set_scale(value, **kwargs)
        self.autoscale_view(scalex=False, scalez=False)
        self._update_transScale()
        self.stale = True
    if maxes.Axes.set_yscale.__doc__ is not None:
        set_yscale.__doc__ = maxes.Axes.set_yscale.__doc__ + """

            .. versionadded :: 1.1.0
                This function was added, but not tested. Please report any bugs.
            """

    @docstring.dedent_interpd
    def set_zscale(self, value, **kwargs):
        """
        Set the scaling of the z-axis: %(scale)s

        ACCEPTS: [%(scale)s]

        Different kwargs are accepted, depending on the scale:
        %(scale_docs)s

        .. note ::
            Currently, Axes3D objects only supports linear scales.
            Other scales may or may not work, and support for these
            is improving with each release.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        self.zaxis._set_scale(value, **kwargs)
        self.autoscale_view(scalex=False, scaley=False)
        self._update_transScale()
        self.stale = True

    def set_zticks(self, *args, **kwargs):
        """
        Set z-axis tick locations.
        See :meth:`matplotlib.axes.Axes.set_yticks` for more details.

        .. note::
            Minor ticks are not supported.

        .. versionadded:: 1.1.0
        """
        return self.zaxis.set_ticks(*args, **kwargs)

    def get_zticks(self, minor=False):
        """
        Return the z ticks as a list of locations
        See :meth:`matplotlib.axes.Axes.get_yticks` for more details.

        .. note::
            Minor ticks are not supported.

        .. versionadded:: 1.1.0
        """
        return self.zaxis.get_ticklocs(minor=minor)

    def get_zmajorticklabels(self):
        """
        Get the ztick labels as a list of Text instances

        .. versionadded :: 1.1.0
        """
        return cbook.silent_list('Text zticklabel',
                                 self.zaxis.get_majorticklabels())

    def get_zminorticklabels(self):
        """
        Get the ztick labels as a list of Text instances

        .. note::
            Minor ticks are not supported. This function was added
            only for completeness.

        .. versionadded :: 1.1.0
        """
        return cbook.silent_list('Text zticklabel',
                                 self.zaxis.get_minorticklabels())

    def set_zticklabels(self, *args, **kwargs):
        """
        Set z-axis tick labels.
        See :meth:`matplotlib.axes.Axes.set_yticklabels` for more details.

        .. note::
            Minor ticks are not supported by Axes3D objects.

        .. versionadded:: 1.1.0
        """
        return self.zaxis.set_ticklabels(*args, **kwargs)

    def get_zticklabels(self, minor=False):
        """
        Get ztick labels as a list of Text instances.
        See :meth:`matplotlib.axes.Axes.get_yticklabels` for more details.

        .. note::
            Minor ticks are not supported.

        .. versionadded:: 1.1.0
        """
        return cbook.silent_list('Text zticklabel',
                                 self.zaxis.get_ticklabels(minor=minor))

    def zaxis_date(self, tz=None):
        """
        Sets up z-axis ticks and labels that treat the z data as dates.

        *tz* is a timezone string or :class:`tzinfo` instance.
        Defaults to rc value.

        .. note::
            This function is merely provided for completeness.
            Axes3D objects do not officially support dates for ticks,
            and so this may or may not work as expected.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        self.zaxis.axis_date(tz)

    def get_zticklines(self):
        """
        Get ztick lines as a list of Line2D instances.
        Note that this function is provided merely for completeness.
        These lines are re-calculated as the display changes.

        .. versionadded:: 1.1.0
        """
        return self.zaxis.get_ticklines()

    def clabel(self, *args, **kwargs):
        """
        This function is currently not implemented for 3D axes.
        Returns *None*.
        """
        return None

    def view_init(self, elev=None, azim=None):
        """
        Set the elevation and azimuth of the axes.

        This can be used to rotate the axes programmatically.

        'elev' stores the elevation angle in the z plane.
        'azim' stores the azimuth angle in the x,y plane.

        if elev or azim are None (default), then the initial value
        is used which was specified in the :class:`Axes3D` constructor.
        """

        self.dist = 10

        if elev is None:
            self.elev = self.initial_elev
        else:
            self.elev = elev

        if azim is None:
            self.azim = self.initial_azim
        else:
            self.azim = azim

    def set_proj_type(self, proj_type):
        """
        Set the projection type.

        Parameters
        ----------
        proj_type : str
            Type of projection, accepts 'persp' and 'ortho'.

        """
        if proj_type == 'persp':
            self._projection = proj3d.persp_transformation
        elif proj_type == 'ortho':
            self._projection = proj3d.ortho_transformation
        else:
            raise ValueError("unrecognized projection: %s" % proj_type)

    def get_proj(self):
        """
        Create the projection matrix from the current viewing position.

        elev stores the elevation angle in the z plane
        azim stores the azimuth angle in the x,y plane

        dist is the distance of the eye viewing point from the object
        point.

        """
        relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

        xmin, xmax = self.get_xlim3d()
        ymin, ymax = self.get_ylim3d()
        zmin, zmax = self.get_zlim3d()

        # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
        worldM = proj3d.world_transformation(xmin, xmax,
                                             ymin, ymax,
                                             zmin, zmax)

        # look into the middle of the new coordinates
        R = np.array([0.5, 0.5, 0.5])

        xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
        yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
        zp = R[2] + np.sin(relev) * self.dist
        E = np.array((xp, yp, zp))

        self.eye = E
        self.vvec = R - E
        self.vvec = self.vvec / proj3d.mod(self.vvec)

        if abs(relev) > np.pi/2:
            # upside down
            V = np.array((0, 0, -1))
        else:
            V = np.array((0, 0, 1))
        zfront, zback = -self.dist, self.dist

        viewM = proj3d.view_transformation(E, R, V)
        projM = self._projection(zfront, zback)
        M0 = np.dot(viewM, worldM)
        M = np.dot(projM, M0)
        return M

    def mouse_init(self, rotate_btn=1, zoom_btn=3):
        """Initializes mouse button callbacks to enable 3D rotation of
        the axes.  Also optionally sets the mouse buttons for 3D rotation
        and zooming.

        ============  =======================================================
        Argument      Description
        ============  =======================================================
        *rotate_btn*  The integer or list of integers specifying which mouse
                      button or buttons to use for 3D rotation of the axes.
                      Default = 1.

        *zoom_btn*    The integer or list of integers specifying which mouse
                      button or buttons to use to zoom the 3D axes.
                      Default = 3.
        ============  =======================================================

        """
        self.button_pressed = None
        canv = self.figure.canvas
        if canv is not None:
            c1 = canv.mpl_connect('motion_notify_event', self._on_move)
            c2 = canv.mpl_connect('button_press_event', self._button_press)
            c3 = canv.mpl_connect('button_release_event', self._button_release)
            self._cids = [c1, c2, c3]
        else:
            warnings.warn(
                "Axes3D.figure.canvas is 'None', mouse rotation disabled.  "
                "Set canvas then call Axes3D.mouse_init().")

        # coerce scalars into array-like, then convert into
        # a regular list to avoid comparisons against None
        # which breaks in recent versions of numpy.
        self._rotate_btn = np.atleast_1d(rotate_btn).tolist()
        self._zoom_btn = np.atleast_1d(zoom_btn).tolist()

    def can_zoom(self):
        """
        Return *True* if this axes supports the zoom box button functionality.

        3D axes objects do not use the zoom box button.
        """
        return False

    def can_pan(self):
        """
        Return *True* if this axes supports the pan/zoom button functionality.

        3D axes objects do not use the pan/zoom button.
        """
        return False

    def cla(self):
        """
        Clear axes
        """
        # Disabling mouse interaction might have been needed a long
        # time ago, but I can't find a reason for it now - BVR (2012-03)
        #self.disable_mouse_rotation()
        super(Axes3D, self).cla()
        self.zaxis.cla()

        if self._sharez is not None:
            self.zaxis.major = self._sharez.zaxis.major
            self.zaxis.minor = self._sharez.zaxis.minor
            z0, z1 = self._sharez.get_zlim()
            self.set_zlim(z0, z1, emit=False, auto=None)
            self.zaxis._set_scale(self._sharez.zaxis.get_scale())
        else:
            self.zaxis._set_scale('linear')
            try:
                self.set_zlim(0, 1)
            except TypeError:
                pass

        self._autoscaleZon = True
        self._zmargin = 0

        self.grid(rcParams['axes3d.grid'])

    def disable_mouse_rotation(self):
        """Disable mouse button callbacks.
        """
        # Disconnect the various events we set.
        for cid in self._cids:
            self.figure.canvas.mpl_disconnect(cid)

        self._cids = []

    def _button_press(self, event):
        if event.inaxes == self:
            self.button_pressed = event.button
            self.sx, self.sy = event.xdata, event.ydata

    def _button_release(self, event):
        self.button_pressed = None

    def format_zdata(self, z):
        """
        Return *z* string formatted.  This function will use the
        :attr:`fmt_zdata` attribute if it is callable, else will fall
        back on the zaxis major formatter
        """
        try: return self.fmt_zdata(z)
        except (AttributeError, TypeError):
            func = self.zaxis.get_major_formatter().format_data_short
            val = func(z)
            return val

    def format_coord(self, xd, yd):
        """
        Given the 2D view coordinates attempt to guess a 3D coordinate.
        Looks for the nearest edge to the point and then assumes that
        the point is at the same z location as the nearest point on the edge.
        """

        if self.M is None:
            return ''

        if self.button_pressed in self._rotate_btn:
            return 'azimuth=%d deg, elevation=%d deg ' % (self.azim, self.elev)
            # ignore xd and yd and display angles instead

        # nearest edge
        p0, p1 = min(self.tunit_edges(),
                     key=lambda edge: proj3d.line2d_seg_dist(
                         edge[0], edge[1], (xd, yd)))

        # scale the z value to match
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        d0 = np.hypot(x0-xd, y0-yd)
        d1 = np.hypot(x1-xd, y1-yd)
        dt = d0+d1
        z = d1/dt * z0 + d0/dt * z1

        x, y, z = proj3d.inv_transform(xd, yd, z, self.M)

        xs = self.format_xdata(x)
        ys = self.format_ydata(y)
        zs = self.format_zdata(z)
        return 'x=%s, y=%s, z=%s' % (xs, ys, zs)

    def _on_move(self, event):
        """Mouse moving

        button-1 rotates by default.  Can be set explicitly in mouse_init().
        button-3 zooms by default.  Can be set explicitly in mouse_init().
        """

        if not self.button_pressed:
            return

        if self.M is None:
            return

        x, y = event.xdata, event.ydata
        # In case the mouse is out of bounds.
        if x is None:
            return

        dx, dy = x - self.sx, y - self.sy
        w = self._pseudo_w
        h = self._pseudo_h
        self.sx, self.sy = x, y

        # Rotation
        if self.button_pressed in self._rotate_btn:
            # rotate viewing point
            # get the x and y pixel coords
            if dx == 0 and dy == 0:
                return
            self.elev = art3d.norm_angle(self.elev - (dy/h)*180)
            self.azim = art3d.norm_angle(self.azim - (dx/w)*180)
            self.get_proj()
            self.stale = True
            self.figure.canvas.draw_idle()

#        elif self.button_pressed == 2:
            # pan view
            # project xv,yv,zv -> xw,yw,zw
            # pan
#            pass

        # Zoom
        elif self.button_pressed in self._zoom_btn:
            # zoom view
            # hmmm..this needs some help from clipping....
            minx, maxx, miny, maxy, minz, maxz = self.get_w_lims()
            df = 1-((h - dy)/h)
            dx = (maxx-minx)*df
            dy = (maxy-miny)*df
            dz = (maxz-minz)*df
            self.set_xlim3d(minx - dx, maxx + dx)
            self.set_ylim3d(miny - dy, maxy + dy)
            self.set_zlim3d(minz - dz, maxz + dz)
            self.get_proj()
            self.figure.canvas.draw_idle()

    def set_zlabel(self, zlabel, fontdict=None, labelpad=None, **kwargs):
        '''
        Set zlabel.  See doc for :meth:`set_ylabel` for description.

        '''
        if labelpad is not None : self.zaxis.labelpad = labelpad
        return self.zaxis.set_label_text(zlabel, fontdict, **kwargs)

    def get_zlabel(self):
        """
        Get the z-label text string.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        label = self.zaxis.get_label()
        return label.get_text()

    #### Axes rectangle characteristics

    def get_frame_on(self):
        """
        Get whether the 3D axes panels are drawn.

        .. versionadded :: 1.1.0
        """
        return self._frameon

    def set_frame_on(self, b):
        """
        Set whether the 3D axes panels are drawn.

        .. versionadded :: 1.1.0

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._frameon = bool(b)
        self.stale = True

    def get_axisbelow(self):
        """
        Get whether axis below is true or not.

        For axes3d objects, this will always be *True*

        .. versionadded :: 1.1.0
            This function was added for completeness.
        """
        return True

    def set_axisbelow(self, b):
        """
        Set whether axis ticks and gridlines are above or below most artists.

        For axes3d objects, this will ignore any settings and just use *True*

        .. versionadded :: 1.1.0
            This function was added for completeness.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._axisbelow = True
        self.stale = True

    def grid(self, b=True, **kwargs):
        '''
        Set / unset 3D grid.

        .. note::

            Currently, this function does not behave the same as
            :meth:`matplotlib.axes.Axes.grid`, but it is intended to
            eventually support that behavior.

        .. versionchanged :: 1.1.0
            This function was changed, but not tested. Please report any bugs.
        '''
        # TODO: Operate on each axes separately
        if len(kwargs):
            b = True
        self._draw_grid = cbook._string_to_bool(b)
        self.stale = True

    def ticklabel_format(self, **kwargs):
        """
        Convenience method for manipulating the ScalarFormatter
        used by default for linear axes in Axed3D objects.

        See :meth:`matplotlib.axes.Axes.ticklabel_format` for full
        documentation.  Note that this version applies to all three
        axes of the Axes3D object.  Therefore, the *axis* argument
        will also accept a value of 'z' and the value of 'both' will
        apply to all three axes.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        style = kwargs.pop('style', '').lower()
        scilimits = kwargs.pop('scilimits', None)
        useOffset = kwargs.pop('useOffset', None)
        axis = kwargs.pop('axis', 'both').lower()
        if scilimits is not None:
            try:
                m, n = scilimits
                m+n+1  # check that both are numbers
            except (ValueError, TypeError):
                raise ValueError("scilimits must be a sequence of 2 integers")
        if style[:3] == 'sci':
            sb = True
        elif style in ['plain', 'comma']:
            sb = False
            if style == 'plain':
                cb = False
            else:
                cb = True
                raise NotImplementedError("comma style remains to be added")
        elif style == '':
            sb = None
        else:
            raise ValueError("%s is not a valid style value")
        try:
            if sb is not None:
                if axis in ['both', 'z']:
                    self.xaxis.major.formatter.set_scientific(sb)
                if axis in ['both', 'y']:
                    self.yaxis.major.formatter.set_scientific(sb)
                if axis in ['both', 'z'] :
                    self.zaxis.major.formatter.set_scientific(sb)
            if scilimits is not None:
                if axis in ['both', 'x']:
                    self.xaxis.major.formatter.set_powerlimits(scilimits)
                if axis in ['both', 'y']:
                    self.yaxis.major.formatter.set_powerlimits(scilimits)
                if axis in ['both', 'z']:
                    self.zaxis.major.formatter.set_powerlimits(scilimits)
            if useOffset is not None:
                if axis in ['both', 'x']:
                    self.xaxis.major.formatter.set_useOffset(useOffset)
                if axis in ['both', 'y']:
                    self.yaxis.major.formatter.set_useOffset(useOffset)
                if axis in ['both', 'z']:
                    self.zaxis.major.formatter.set_useOffset(useOffset)
        except AttributeError:
            raise AttributeError(
                "This method only works with the ScalarFormatter.")

    def locator_params(self, axis='both', tight=None, **kwargs):
        """
        Convenience method for controlling tick locators.

        See :meth:`matplotlib.axes.Axes.locator_params` for full
        documentation  Note that this is for Axes3D objects,
        therefore, setting *axis* to 'both' will result in the
        parameters being set for all three axes.  Also, *axis*
        can also take a value of 'z' to apply parameters to the
        z axis.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        _x = axis in ['x', 'both']
        _y = axis in ['y', 'both']
        _z = axis in ['z', 'both']
        if _x:
            self.xaxis.get_major_locator().set_params(**kwargs)
        if _y:
            self.yaxis.get_major_locator().set_params(**kwargs)
        if _z:
            self.zaxis.get_major_locator().set_params(**kwargs)
        self.autoscale_view(tight=tight, scalex=_x, scaley=_y, scalez=_z)

    def tick_params(self, axis='both', **kwargs):
        """
        Convenience method for changing the appearance of ticks and
        tick labels.

        See :meth:`matplotlib.axes.Axes.tick_params` for more complete
        documentation.

        The only difference is that setting *axis* to 'both' will
        mean that the settings are applied to all three axes. Also,
        the *axis* parameter also accepts a value of 'z', which
        would mean to apply to only the z-axis.

        Also, because of how Axes3D objects are drawn very differently
        from regular 2D axes, some of these settings may have
        ambiguous meaning.  For simplicity, the 'z' axis will
        accept settings as if it was like the 'y' axis.

        .. note::
            While this function is currently implemented, the core part
            of the Axes3D object may ignore some of these settings.
            Future releases will fix this. Priority will be given to
            those who file bugs.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        super(Axes3D, self).tick_params(axis, **kwargs)
        if axis in ['z', 'both'] :
            zkw = dict(kwargs)
            zkw.pop('top', None)
            zkw.pop('bottom', None)
            zkw.pop('labeltop', None)
            zkw.pop('labelbottom', None)
            self.zaxis.set_tick_params(**zkw)

    ### data limits, ticks, tick labels, and formatting

    def invert_zaxis(self):
        """
        Invert the z-axis.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        bottom, top = self.get_zlim()
        self.set_zlim(top, bottom, auto=None)

    def zaxis_inverted(self):
        '''
        Returns True if the z-axis is inverted.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        '''
        bottom, top = self.get_zlim()
        return top < bottom

    def get_zbound(self):
        """
        Returns the z-axis numerical bounds where::

          lowerBound < upperBound

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        bottom, top = self.get_zlim()
        if bottom < top:
            return bottom, top
        else:
            return top, bottom

    def set_zbound(self, lower=None, upper=None):
        """
        Set the lower and upper numerical bounds of the z-axis.
        This method will honor axes inversion regardless of parameter order.
        It will not change the :attr:`_autoscaleZon` attribute.

        .. versionadded :: 1.1.0
            This function was added, but not tested. Please report any bugs.
        """
        if upper is None and cbook.iterable(lower):
            lower,upper = lower

        old_lower,old_upper = self.get_zbound()

        if lower is None: lower = old_lower
        if upper is None: upper = old_upper

        if self.zaxis_inverted():
            if lower < upper:
                self.set_zlim(upper, lower, auto=None)
            else:
                self.set_zlim(lower, upper, auto=None)
        else :
            if lower < upper:
                self.set_zlim(lower, upper, auto=None)
            else :
                self.set_zlim(upper, lower, auto=None)

    def text(self, x, y, z, s, zdir=None, **kwargs):
        '''
        Add text to the plot. kwargs will be passed on to Axes.text,
        except for the `zdir` keyword, which sets the direction to be
        used as the z direction.
        '''
        text = super(Axes3D, self).text(x, y, s, **kwargs)
        art3d.text_2d_to_3d(text, z, zdir)
        return text

    text3D = text
    text2D = Axes.text

    def plot(self, xs, ys, *args, **kwargs):
        '''
        Plot 2D or 3D data.

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *xs*, *ys*  x, y coordinates of vertices

        *zs*        z value(s), either one for all points or one for
                    each point.
        *zdir*      Which direction to use as z ('x', 'y' or 'z')
                    when plotting a 2D set.
        ==========  ================================================

        Other arguments are passed on to
        :func:`~matplotlib.axes.Axes.plot`
        '''
        had_data = self.has_data()

        # `zs` can be passed positionally or as keyword; checking whether
        # args[0] is a string matches the behavior of 2D `plot` (via
        # `_process_plot_var_args`).
        if args and not isinstance(args[0], six.string_types):
            zs = args[0]
            args = args[1:]
            if 'zs' in kwargs:
                raise TypeError("plot() for multiple values for argument 'z'")
        else:
            zs = kwargs.pop('zs', 0)
        zdir = kwargs.pop('zdir', 'z')

        # Match length
        zs = _backports.broadcast_to(zs, len(xs))

        lines = super(Axes3D, self).plot(xs, ys, *args, **kwargs)
        for line in lines:
            art3d.line_2d_to_3d(line, zs=zs, zdir=zdir)

        xs, ys, zs = art3d.juggle_axes(xs, ys, zs, zdir)
        self.auto_scale_xyz(xs, ys, zs, had_data)
        return lines

    plot3D = plot

    def plot_surface(self, X, Y, Z, *args, **kwargs):
        """
        Create a surface plot.

        By default it will be colored in shades of a solid color, but it also
        supports color mapping by supplying the *cmap* argument.

        .. note::

           The *rcount* and *ccount* kwargs, which both default to 50,
           determine the maximum number of samples used in each direction.  If
           the input data is larger, it will be downsampled (by slicing) to
           these numbers of points.

        Parameters
        ----------
        X, Y, Z : 2d arrays
            Data values.

        rcount, ccount : int
            Maximum number of samples used in each direction.  If the input
            data is larger, it will be downsampled (by slicing) to these
            numbers of points.  Defaults to 50.

            .. versionadded:: 2.0

        rstride, cstride : int
            Downsampling stride in each direction.  These arguments are
            mutually exclusive with *rcount* and *ccount*.  If only one of
            *rstride* or *cstride* is set, the other defaults to 10.

            'classic' mode uses a default of ``rstride = cstride = 10`` instead
            of the new default of ``rcount = ccount = 50``.

        color : color-like
            Color of the surface patches.

        cmap : Colormap
            Colormap of the surface patches.

        facecolors : array-like of colors.
            Colors of each individual patch.

        norm : Normalize
            Normalization for the colormap.

        vmin, vmax : float
            Bounds for the normalization.

        shade : bool
            Whether to shade the face colors.

        **kwargs :
            Other arguments are forwarded to `.Poly3DCollection`.
        """

        had_data = self.has_data()

        if Z.ndim != 2:
            raise ValueError("Argument Z must be 2-dimensional.")
        # TODO: Support masked arrays
        X, Y, Z = np.broadcast_arrays(X, Y, Z)
        rows, cols = Z.shape

        has_stride = 'rstride' in kwargs or 'cstride' in kwargs
        has_count = 'rcount' in kwargs or 'ccount' in kwargs

        if has_stride and has_count:
            raise ValueError("Cannot specify both stride and count arguments")

        rstride = kwargs.pop('rstride', 10)
        cstride = kwargs.pop('cstride', 10)
        rcount = kwargs.pop('rcount', 50)
        ccount = kwargs.pop('ccount', 50)

        if rcParams['_internal.classic_mode']:
            # Strides have priority over counts in classic mode.
            # So, only compute strides from counts
            # if counts were explicitly given
            if has_count:
                rstride = int(max(np.ceil(rows / rcount), 1))
                cstride = int(max(np.ceil(cols / ccount), 1))
        else:
            # If the strides are provided then it has priority.
            # Otherwise, compute the strides from the counts.
            if not has_stride:
                rstride = int(max(np.ceil(rows / rcount), 1))
                cstride = int(max(np.ceil(cols / ccount), 1))

        if 'facecolors' in kwargs:
            fcolors = kwargs.pop('facecolors')
        else:
            color = kwargs.pop('color', None)
            if color is None:
                color = self._get_lines.get_next_color()
            color = np.array(mcolors.to_rgba(color))
            fcolors = None

        cmap = kwargs.get('cmap', None)
        norm = kwargs.pop('norm', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        linewidth = kwargs.get('linewidth', None)
        shade = kwargs.pop('shade', cmap is None)
        lightsource = kwargs.pop('lightsource', None)

        # Shade the data
        if shade and cmap is not None and fcolors is not None:
            fcolors = self._shade_colors_lightsource(Z, cmap, lightsource)

        polys = []
        # Only need these vectors to shade if there is no cmap
        if cmap is None and shade :
            totpts = int(np.ceil((rows - 1) / rstride) *
                         np.ceil((cols - 1) / cstride))
            v1 = np.empty((totpts, 3))
            v2 = np.empty((totpts, 3))
            # This indexes the vertex points
            which_pt = 0


        #colset contains the data for coloring: either average z or the facecolor
        colset = []
        for rs in xrange(0, rows-1, rstride):
            for cs in xrange(0, cols-1, cstride):
                ps = []
                for a in (X, Y, Z):
                    ztop = a[rs,cs:min(cols, cs+cstride+1)]
                    zleft = a[rs+1:min(rows, rs+rstride+1),
                              min(cols-1, cs+cstride)]
                    zbase = a[min(rows-1, rs+rstride), cs:min(cols, cs+cstride+1):][::-1]
                    zright = a[rs:min(rows-1, rs+rstride):, cs][::-1]
                    z = np.concatenate((ztop, zleft, zbase, zright))
                    ps.append(z)

                # The construction leaves the array with duplicate points, which
                # are removed here.
                ps = list(zip(*ps))
                lastp = np.array([])
                ps2 = [ps[0]] + [ps[i] for i in xrange(1, len(ps)) if ps[i] != ps[i-1]]
                avgzsum = sum(p[2] for p in ps2)
                polys.append(ps2)

                if fcolors is not None:
                    colset.append(fcolors[rs][cs])
                else:
                    colset.append(avgzsum / len(ps2))

                # Only need vectors to shade if no cmap
                if cmap is None and shade:
                    i1, i2, i3 = 0, int(len(ps2)/3), int(2*len(ps2)/3)
                    v1[which_pt] = np.array(ps2[i1]) - np.array(ps2[i2])
                    v2[which_pt] = np.array(ps2[i2]) - np.array(ps2[i3])
                    which_pt += 1
        if cmap is None and shade:
            normals = np.cross(v1, v2)
        else :
            normals = []

        polyc = art3d.Poly3DCollection(polys, *args, **kwargs)

        if fcolors is not None:
            if shade:
                colset = self._shade_colors(colset, normals)
            polyc.set_facecolors(colset)
            polyc.set_edgecolors(colset)
        elif cmap:
            colset = np.array(colset)
            polyc.set_array(colset)
            if vmin is not None or vmax is not None:
                polyc.set_clim(vmin, vmax)
            if norm is not None:
                polyc.set_norm(norm)
        else:
            if shade:
                colset = self._shade_colors(color, normals)
            else:
                colset = color
            polyc.set_facecolors(colset)

        self.add_collection(polyc)
        self.auto_scale_xyz(X, Y, Z, had_data)

        return polyc

    def _generate_normals(self, polygons):
        '''
        Generate normals for polygons by using the first three points.
        This normal of course might not make sense for polygons with
        more than three points not lying in a plane.
        '''

        normals = []
        for verts in polygons:
            v1 = np.array(verts[0]) - np.array(verts[1])
            v2 = np.array(verts[2]) - np.array(verts[0])
            normals.append(np.cross(v1, v2))
        return normals

    def _shade_colors(self, color, normals):
        '''
        Shade *color* using normal vectors given by *normals*.
        *color* can also be an array of the same length as *normals*.
        '''

        shade = np.array([np.dot(n / proj3d.mod(n), [-1, -1, 0.5])
                          if proj3d.mod(n) else np.nan
                          for n in normals])
        mask = ~np.isnan(shade)

        if len(shade[mask]) > 0:
            norm = Normalize(min(shade[mask]), max(shade[mask]))
            shade[~mask] = min(shade[mask])
            color = mcolors.to_rgba_array(color)
            # shape of color should be (M, 4) (where M is number of faces)
            # shape of shade should be (M,)
            # colors should have final shape of (M, 4)
            alpha = color[:, 3]
            colors = (0.5 + norm(shade)[:, np.newaxis] * 0.5) * color
            colors[:, 3] = alpha
        else:
            colors = np.asanyarray(color).copy()

        return colors

    def _shade_colors_lightsource(self, data, cmap, lightsource):
        if lightsource is None:
            lightsource = LightSource(azdeg=135, altdeg=55)
        return lightsource.shade(data, cmap)

    def plot_wireframe(self, X, Y, Z, *args, **kwargs):
        """
        Plot a 3D wireframe.

        .. note::

           The *rcount* and *ccount* kwargs, which both default to 50,
           determine the maximum number of samples used in each direction.  If
           the input data is larger, it will be downsampled (by slicing) to
           these numbers of points.

        Parameters
        ----------
        X, Y, Z : 2d arrays
            Data values.

        rcount, ccount : int
            Maximum number of samples used in each direction.  If the input
            data is larger, it will be downsampled (by slicing) to these
            numbers of points.  Setting a count to zero causes the data to be
            not sampled in the corresponding direction, producing a 3D line
            plot rather than a wireframe plot.  Defaults to 50.

            .. versionadded:: 2.0

        rstride, cstride : int
            Downsampling stride in each direction.  These arguments are
            mutually exclusive with *rcount* and *ccount*.  If only one of
            *rstride* or *cstride* is set, the other defaults to 1.  Setting a
            stride to zero causes the data to be not sampled in the
            corresponding direction, producing a 3D line plot rather than a
            wireframe plot.

            'classic' mode uses a default of ``rstride = cstride = 1`` instead
            of the new default of ``rcount = ccount = 50``.

        **kwargs :
            Other arguments are forwarded to `.Line3DCollection`.
        """

        had_data = self.has_data()
        if Z.ndim != 2:
            raise ValueError("Argument Z must be 2-dimensional.")
        # FIXME: Support masked arrays
        X, Y, Z = np.broadcast_arrays(X, Y, Z)
        rows, cols = Z.shape

        has_stride = 'rstride' in kwargs or 'cstride' in kwargs
        has_count = 'rcount' in kwargs or 'ccount' in kwargs

        if has_stride and has_count:
            raise ValueError("Cannot specify both stride and count arguments")

        rstride = kwargs.pop('rstride', 1)
        cstride = kwargs.pop('cstride', 1)
        rcount = kwargs.pop('rcount', 50)
        ccount = kwargs.pop('ccount', 50)

        if rcParams['_internal.classic_mode']:
            # Strides have priority over counts in classic mode.
            # So, only compute strides from counts
            # if counts were explicitly given
            if has_count:
                rstride = int(max(np.ceil(rows / rcount), 1)) if rcount else 0
                cstride = int(max(np.ceil(cols / ccount), 1)) if ccount else 0
        else:
            # If the strides are provided then it has priority.
            # Otherwise, compute the strides from the counts.
            if not has_stride:
                rstride = int(max(np.ceil(rows / rcount), 1)) if rcount else 0
                cstride = int(max(np.ceil(cols / ccount), 1)) if ccount else 0

        # We want two sets of lines, one running along the "rows" of
        # Z and another set of lines running along the "columns" of Z.
        # This transpose will make it easy to obtain the columns.
        tX, tY, tZ = np.transpose(X), np.transpose(Y), np.transpose(Z)

        if rstride:
            rii = list(xrange(0, rows, rstride))
            # Add the last index only if needed
            if rows > 0 and rii[-1] != (rows - 1):
                rii += [rows-1]
        else:
            rii = []
        if cstride:
            cii = list(xrange(0, cols, cstride))
            # Add the last index only if needed
            if cols > 0 and cii[-1] != (cols - 1):
                cii += [cols-1]
        else:
            cii = []

        if rstride == 0 and cstride == 0:
            raise ValueError("Either rstride or cstride must be non zero")

        # If the inputs were empty, then just
        # reset everything.
        if Z.size == 0:
            rii = []
            cii = []

        xlines = [X[i] for i in rii]
        ylines = [Y[i] for i in rii]
        zlines = [Z[i] for i in rii]

        txlines = [tX[i] for i in cii]
        tylines = [tY[i] for i in cii]
        tzlines = [tZ[i] for i in cii]

        lines = ([list(zip(xl, yl, zl))
                  for xl, yl, zl in zip(xlines, ylines, zlines)]
                + [list(zip(xl, yl, zl))
                   for xl, yl, zl in zip(txlines, tylines, tzlines)])

        linec = art3d.Line3DCollection(lines, *args, **kwargs)
        self.add_collection(linec)
        self.auto_scale_xyz(X, Y, Z, had_data)

        return linec

    def plot_trisurf(self, *args, **kwargs):
        """
        ============= ================================================
        Argument      Description
        ============= ================================================
        *X*, *Y*, *Z* Data values as 1D arrays
        *color*       Color of the surface patches
        *cmap*        A colormap for the surface patches.
        *norm*        An instance of Normalize to map values to colors
        *vmin*        Minimum value to map
        *vmax*        Maximum value to map
        *shade*       Whether to shade the facecolors
        ============= ================================================

        The (optional) triangulation can be specified in one of two ways;
        either::

          plot_trisurf(triangulation, ...)

        where triangulation is a :class:`~matplotlib.tri.Triangulation`
        object, or::

          plot_trisurf(X, Y, ...)
          plot_trisurf(X, Y, triangles, ...)
          plot_trisurf(X, Y, triangles=triangles, ...)

        in which case a Triangulation object will be created.  See
        :class:`~matplotlib.tri.Triangulation` for a explanation of
        these possibilities.

        The remaining arguments are::

          plot_trisurf(..., Z)

        where *Z* is the array of values to contour, one per point
        in the triangulation.

        Other arguments are passed on to
        :class:`~mpl_toolkits.mplot3d.art3d.Poly3DCollection`

        **Examples:**

        .. plot:: gallery/mplot3d/trisurf3d.py
        .. plot:: gallery/mplot3d/trisurf3d_2.py

        .. versionadded:: 1.2.0
            This plotting function was added for the v1.2.0 release.
        """

        had_data = self.has_data()

        # TODO: Support custom face colours
        color = kwargs.pop('color', None)
        if color is None:
            color = self._get_lines.get_next_color()
        color = np.array(mcolors.to_rgba(color))

        cmap = kwargs.get('cmap', None)
        norm = kwargs.pop('norm', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        linewidth = kwargs.get('linewidth', None)
        shade = kwargs.pop('shade', cmap is None)
        lightsource = kwargs.pop('lightsource', None)

        tri, args, kwargs = Triangulation.get_from_args_and_kwargs(*args, **kwargs)
        if 'Z' in kwargs:
            z = np.asarray(kwargs.pop('Z'))
        else:
            z = np.asarray(args[0])
            # We do this so Z doesn't get passed as an arg to PolyCollection
            args = args[1:]

        triangles = tri.get_masked_triangles()
        xt = tri.x[triangles]
        yt = tri.y[triangles]
        zt = z[triangles]

        # verts = np.stack((xt, yt, zt), axis=-1)
        verts = np.concatenate((
            xt[..., np.newaxis], yt[..., np.newaxis], zt[..., np.newaxis]
        ), axis=-1)

        polyc = art3d.Poly3DCollection(verts, *args, **kwargs)

        if cmap:
            # average over the three points of each triangle
            avg_z = verts[:, :, 2].mean(axis=1)
            polyc.set_array(avg_z)
            if vmin is not None or vmax is not None:
                polyc.set_clim(vmin, vmax)
            if norm is not None:
                polyc.set_norm(norm)
        else:
            if shade:
                v1 = verts[:, 0, :] - verts[:, 1, :]
                v2 = verts[:, 1, :] - verts[:, 2, :]
                normals = np.cross(v1, v2)
                colset = self._shade_colors(color, normals)
            else:
                colset = color
            polyc.set_facecolors(colset)

        self.add_collection(polyc)
        self.auto_scale_xyz(tri.x, tri.y, z, had_data)

        return polyc

    def _3d_extend_contour(self, cset, stride=5):
        '''
        Extend a contour in 3D by creating
        '''

        levels = cset.levels
        colls = cset.collections
        dz = (levels[1] - levels[0]) / 2

        for z, linec in zip(levels, colls):
            topverts = art3d.paths_to_3d_segments(linec.get_paths(), z - dz)
            botverts = art3d.paths_to_3d_segments(linec.get_paths(), z + dz)

            color = linec.get_color()[0]

            polyverts = []
            normals = []
            nsteps = np.round(len(topverts[0]) / stride)
            if nsteps <= 1:
                if len(topverts[0]) > 1:
                    nsteps = 2
                else:
                    continue

            stepsize = (len(topverts[0]) - 1) / (nsteps - 1)
            for i in range(int(np.round(nsteps)) - 1):
                i1 = int(np.round(i * stepsize))
                i2 = int(np.round((i + 1) * stepsize))
                polyverts.append([topverts[0][i1],
                    topverts[0][i2],
                    botverts[0][i2],
                    botverts[0][i1]])

                v1 = np.array(topverts[0][i1]) - np.array(topverts[0][i2])
                v2 = np.array(topverts[0][i1]) - np.array(botverts[0][i1])
                normals.append(np.cross(v1, v2))

            colors = self._shade_colors(color, normals)
            colors2 = self._shade_colors(color, normals)
            polycol = art3d.Poly3DCollection(polyverts,
                                             facecolors=colors,
                                             edgecolors=colors2)
            polycol.set_sort_zpos(z)
            self.add_collection3d(polycol)

        for col in colls:
            self.collections.remove(col)

    def add_contour_set(self, cset, extend3d=False, stride=5, zdir='z', offset=None):
        zdir = '-' + zdir
        if extend3d:
            self._3d_extend_contour(cset, stride)
        else:
            for z, linec in zip(cset.levels, cset.collections):
                if offset is not None:
                    z = offset
                art3d.line_collection_2d_to_3d(linec, z, zdir=zdir)

    def add_contourf_set(self, cset, zdir='z', offset=None):
        zdir = '-' + zdir
        for z, linec in zip(cset.levels, cset.collections):
            if offset is not None :
                z = offset
            art3d.poly_collection_2d_to_3d(linec, z, zdir=zdir)
            linec.set_sort_zpos(z)

    def contour(self, X, Y, Z, *args, **kwargs):
        '''
        Create a 3D contour plot.

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *X*, *Y*,   Data values as numpy.arrays
        *Z*
        *extend3d*  Whether to extend contour in 3D (default: False)
        *stride*    Stride (step size) for extending contour
        *zdir*      The direction to use: x, y or z (default)
        *offset*    If specified plot a projection of the contour
                    lines on this position in plane normal to zdir
        ==========  ================================================

        The positional and other keyword arguments are passed on to
        :func:`~matplotlib.axes.Axes.contour`

        Returns a :class:`~matplotlib.axes.Axes.contour`
        '''

        extend3d = kwargs.pop('extend3d', False)
        stride = kwargs.pop('stride', 5)
        zdir = kwargs.pop('zdir', 'z')
        offset = kwargs.pop('offset', None)

        had_data = self.has_data()

        jX, jY, jZ = art3d.rotate_axes(X, Y, Z, zdir)
        cset = super(Axes3D, self).contour(jX, jY, jZ, *args, **kwargs)
        self.add_contour_set(cset, extend3d, stride, zdir, offset)

        self.auto_scale_xyz(X, Y, Z, had_data)
        return cset

    contour3D = contour

    def tricontour(self, *args, **kwargs):
        """
        Create a 3D contour plot.

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *X*, *Y*,   Data values as numpy.arrays
        *Z*
        *extend3d*  Whether to extend contour in 3D (default: False)
        *stride*    Stride (step size) for extending contour
        *zdir*      The direction to use: x, y or z (default)
        *offset*    If specified plot a projection of the contour
                    lines on this position in plane normal to zdir
        ==========  ================================================

        Other keyword arguments are passed on to
        :func:`~matplotlib.axes.Axes.tricontour`

        Returns a :class:`~matplotlib.axes.Axes.contour`

        .. versionchanged:: 1.3.0
            Added support for custom triangulations

        EXPERIMENTAL:  This method currently produces incorrect output due to a
        longstanding bug in 3D PolyCollection rendering.
        """

        extend3d = kwargs.pop('extend3d', False)
        stride = kwargs.pop('stride', 5)
        zdir = kwargs.pop('zdir', 'z')
        offset = kwargs.pop('offset', None)

        had_data = self.has_data()

        tri, args, kwargs = Triangulation.get_from_args_and_kwargs(
                *args, **kwargs)
        X = tri.x
        Y = tri.y
        if 'Z' in kwargs:
            Z = kwargs.pop('Z')
        else:
            Z = args[0]
            # We do this so Z doesn't get passed as an arg to Axes.tricontour
            args = args[1:]

        jX, jY, jZ = art3d.rotate_axes(X, Y, Z, zdir)
        tri = Triangulation(jX, jY, tri.triangles, tri.mask)

        cset = super(Axes3D, self).tricontour(tri, jZ, *args, **kwargs)
        self.add_contour_set(cset, extend3d, stride, zdir, offset)

        self.auto_scale_xyz(X, Y, Z, had_data)
        return cset

    def contourf(self, X, Y, Z, *args, **kwargs):
        '''
        Create a 3D contourf plot.

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *X*, *Y*,   Data values as numpy.arrays
        *Z*
        *zdir*      The direction to use: x, y or z (default)
        *offset*    If specified plot a projection of the filled contour
                    on this position in plane normal to zdir
        ==========  ================================================

        The positional and keyword arguments are passed on to
        :func:`~matplotlib.axes.Axes.contourf`

        Returns a :class:`~matplotlib.axes.Axes.contourf`

        .. versionchanged :: 1.1.0
            The *zdir* and *offset* kwargs were added.
        '''

        zdir = kwargs.pop('zdir', 'z')
        offset = kwargs.pop('offset', None)

        had_data = self.has_data()

        jX, jY, jZ = art3d.rotate_axes(X, Y, Z, zdir)
        cset = super(Axes3D, self).contourf(jX, jY, jZ, *args, **kwargs)
        self.add_contourf_set(cset, zdir, offset)

        self.auto_scale_xyz(X, Y, Z, had_data)
        return cset

    contourf3D = contourf

    def tricontourf(self, *args, **kwargs):
        """
        Create a 3D contourf plot.

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *X*, *Y*,   Data values as numpy.arrays
        *Z*
        *zdir*      The direction to use: x, y or z (default)
        *offset*    If specified plot a projection of the contour
                    lines on this position in plane normal to zdir
        ==========  ================================================

        Other keyword arguments are passed on to
        :func:`~matplotlib.axes.Axes.tricontour`

        Returns a :class:`~matplotlib.axes.Axes.contour`

        .. versionchanged :: 1.3.0
            Added support for custom triangulations

        EXPERIMENTAL:  This method currently produces incorrect output due to a
        longstanding bug in 3D PolyCollection rendering.
        """
        zdir = kwargs.pop('zdir', 'z')
        offset = kwargs.pop('offset', None)

        had_data = self.has_data()

        tri, args, kwargs = Triangulation.get_from_args_and_kwargs(
                *args, **kwargs)
        X = tri.x
        Y = tri.y
        if 'Z' in kwargs:
            Z = kwargs.pop('Z')
        else:
            Z = args[0]
            # We do this so Z doesn't get passed as an arg to Axes.tricontourf
            args = args[1:]

        jX, jY, jZ = art3d.rotate_axes(X, Y, Z, zdir)
        tri = Triangulation(jX, jY, tri.triangles, tri.mask)

        cset = super(Axes3D, self).tricontourf(tri, jZ, *args, **kwargs)
        self.add_contourf_set(cset, zdir, offset)

        self.auto_scale_xyz(X, Y, Z, had_data)
        return cset

    def add_collection3d(self, col, zs=0, zdir='z'):
        '''
        Add a 3D collection object to the plot.

        2D collection types are converted to a 3D version by
        modifying the object and adding z coordinate information.

        Supported are:
            - PolyCollection
            - LineCollection
            - PatchCollection
        '''
        zvals = np.atleast_1d(zs)
        if len(zvals) > 0 :
            zsortval = min(zvals)
        else :
            zsortval = 0   # FIXME: Fairly arbitrary. Is there a better value?

        # FIXME: use issubclass() (although, then a 3D collection
        #       object would also pass.)  Maybe have a collection3d
        #       abstract class to test for and exclude?
        if type(col) is mcoll.PolyCollection:
            art3d.poly_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)
        elif type(col) is mcoll.LineCollection:
            art3d.line_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)
        elif type(col) is mcoll.PatchCollection:
            art3d.patch_collection_2d_to_3d(col, zs=zs, zdir=zdir)
            col.set_sort_zpos(zsortval)

        super(Axes3D, self).add_collection(col)

    def scatter(self, xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True,
                *args, **kwargs):
        '''
        Create a scatter plot.

        ============  ========================================================
        Argument      Description
        ============  ========================================================
        *xs*, *ys*    Positions of data points.
        *zs*          Either an array of the same length as *xs* and
                      *ys* or a single value to place all points in
                      the same plane. Default is 0.
        *zdir*        Which direction to use as z ('x', 'y' or 'z')
                      when plotting a 2D set.
        *s*           Size in points^2.  It is a scalar or an array of the
                      same length as *x* and *y*.

        *c*           A color. *c* can be a single color format string, or a
                      sequence of color specifications of length *N*, or a
                      sequence of *N* numbers to be mapped to colors using the
                      *cmap* and *norm* specified via kwargs (see below). Note
                      that *c* should not be a single numeric RGB or RGBA
                      sequence because that is indistinguishable from an array
                      of values to be colormapped.  *c* can be a 2-D array in
                      which the rows are RGB or RGBA, however, including the
                      case of a single row to specify the same color for
                      all points.

        *depthshade*
                      Whether or not to shade the scatter markers to give
                      the appearance of depth. Default is *True*.
        ============  ========================================================

        Keyword arguments are passed on to
        :func:`~matplotlib.axes.Axes.scatter`.

        Returns a :class:`~mpl_toolkits.mplot3d.art3d.Patch3DCollection`
        '''

        had_data = self.has_data()

        xs, ys, zs = np.broadcast_arrays(
            *[np.ravel(np.ma.filled(t, np.nan)) for t in [xs, ys, zs]])
        s = np.ma.ravel(s)  # This doesn't have to match x, y in size.

        xs, ys, zs, s, c = cbook.delete_masked_points(xs, ys, zs, s, c)

        patches = super(Axes3D, self).scatter(
            xs, ys, s=s, c=c, *args, **kwargs)
        is_2d = not cbook.iterable(zs)
        zs = _backports.broadcast_to(zs, len(xs))
        art3d.patch_collection_2d_to_3d(patches, zs=zs, zdir=zdir,
                                        depthshade=depthshade)

        if self._zmargin < 0.05 and xs.size > 0:
            self.set_zmargin(0.05)

        #FIXME: why is this necessary?
        if not is_2d:
            self.auto_scale_xyz(xs, ys, zs, had_data)

        return patches

    scatter3D = scatter

    def bar(self, left, height, zs=0, zdir='z', *args, **kwargs):
        '''
        Add 2D bar(s).

        ==========  ================================================
        Argument    Description
        ==========  ================================================
        *left*      The x coordinates of the left sides of the bars.
        *height*    The height of the bars.
        *zs*        Z coordinate of bars, if one value is specified
                    they will all be placed at the same z.
        *zdir*      Which direction to use as z ('x', 'y' or 'z')
                    when plotting a 2D set.
        ==========  ================================================

        Keyword arguments are passed onto :func:`~matplotlib.axes.Axes.bar`.

        Returns a :class:`~mpl_toolkits.mplot3d.art3d.Patch3DCollection`
        '''

        had_data = self.has_data()

        patches = super(Axes3D, self).bar(left, height, *args, **kwargs)

        zs = _backports.broadcast_to(zs, len(left))

        verts = []
        verts_zs = []
        for p, z in zip(patches, zs):
            vs = art3d.get_patch_verts(p)
            verts += vs.tolist()
            verts_zs += [z] * len(vs)
            art3d.patch_2d_to_3d(p, z, zdir)
            if 'alpha' in kwargs:
                p.set_alpha(kwargs['alpha'])

        if len(verts) > 0 :
            # the following has to be skipped if verts is empty
            # NOTE: Bugs could still occur if len(verts) > 0,
            #       but the "2nd dimension" is empty.
            xs, ys = list(zip(*verts))
        else :
            xs, ys = [], []

        xs, ys, verts_zs = art3d.juggle_axes(xs, ys, verts_zs, zdir)
        self.auto_scale_xyz(xs, ys, verts_zs, had_data)

        return patches

    def bar3d(self, x, y, z, dx, dy, dz, color=None,
              zsort='average', shade=True, *args, **kwargs):
        """Generate a 3D barplot.

        This method creates three dimensional barplot where the width,
        depth, height, and color of the bars can all be uniquely set.

        Parameters
        ----------
        x, y, z : array-like
            The coordinates of the anchor point of the bars.

        dx, dy, dz : scalar or array-like
            The width, depth, and height of the bars, respectively.

        color : sequence of valid color specifications, optional
            The color of the bars can be specified globally or
            individually. This parameter can be:

              - A single color value, to color all bars the same color.
              - An array of colors of length N bars, to color each bar
                independently.
              - An array of colors of length 6, to color the faces of the
                bars similarly.
              - An array of colors of length 6 * N bars, to color each face
                independently.

            When coloring the faces of the boxes specifically, this is
            the order of the coloring:

              1. -Z (bottom of box)
              2. +Z (top of box)
              3. -Y
              4. +Y
              5. -X
              6. +X

        zsort : str, optional
            The z-axis sorting scheme passed onto
            :func:`~mpl_toolkits.mplot3d.art3d.Poly3DCollection`

        shade : bool, optional (default = True)
            When true, this shades the dark sides of the bars (relative
            to the plot's source of light).

        Any additional keyword arguments are passed onto
        :func:`~mpl_toolkits.mplot3d.art3d.Poly3DCollection`

        Returns
        -------
        collection : Poly3DCollection
            A collection of three dimensional polygons representing
            the bars.
        """

        had_data = self.has_data()

        x, y, z, dx, dy, dz = np.broadcast_arrays(
            np.atleast_1d(x), y, z, dx, dy, dz)
        minx = np.min(x)
        maxx = np.max(x + dx)
        miny = np.min(y)
        maxy = np.max(y + dy)
        minz = np.min(z)
        maxz = np.max(z + dz)

        polys = []
        for xi, yi, zi, dxi, dyi, dzi in zip(x, y, z, dx, dy, dz):
            polys.extend([
                ((xi, yi, zi), (xi + dxi, yi, zi),
                    (xi + dxi, yi + dyi, zi), (xi, yi + dyi, zi)),
                ((xi, yi, zi + dzi), (xi + dxi, yi, zi + dzi),
                    (xi + dxi, yi + dyi, zi + dzi), (xi, yi + dyi, zi + dzi)),

                ((xi, yi, zi), (xi + dxi, yi, zi),
                    (xi + dxi, yi, zi + dzi), (xi, yi, zi + dzi)),
                ((xi, yi + dyi, zi), (xi + dxi, yi + dyi, zi),
                    (xi + dxi, yi + dyi, zi + dzi), (xi, yi + dyi, zi + dzi)),

                ((xi, yi, zi), (xi, yi + dyi, zi),
                    (xi, yi + dyi, zi + dzi), (xi, yi, zi + dzi)),
                ((xi + dxi, yi, zi), (xi + dxi, yi + dyi, zi),
                    (xi + dxi, yi + dyi, zi + dzi), (xi + dxi, yi, zi + dzi)),
            ])

        facecolors = []
        if color is None:
            color = [self._get_patches_for_fill.get_next_color()]

        if len(color) == len(x):
            # bar colors specified, need to expand to number of faces
            for c in color:
                facecolors.extend([c] * 6)
        else:
            # a single color specified, or face colors specified explicitly
            facecolors = list(mcolors.to_rgba_array(color))
            if len(facecolors) < len(x):
                facecolors *= (6 * len(x))

        if shade:
            normals = self._generate_normals(polys)
            sfacecolors = self._shade_colors(facecolors, normals)
        else:
            sfacecolors = facecolors

        col = art3d.Poly3DCollection(polys,
                                     zsort=zsort,
                                     facecolor=sfacecolors,
                                     *args, **kwargs)
        self.add_collection(col)

        self.auto_scale_xyz((minx, maxx), (miny, maxy), (minz, maxz), had_data)

        return col

    def set_title(self, label, fontdict=None, loc='center', **kwargs):
        ret = super(Axes3D, self).set_title(label, fontdict=fontdict, loc=loc,
                                            **kwargs)
        (x, y) = self.title.get_position()
        self.title.set_y(0.92 * y)
        return ret
    set_title.__doc__ = maxes.Axes.set_title.__doc__

    def quiver(self, *args, **kwargs):
        """
        Plot a 3D field of arrows.

        call signatures::

            quiver(X, Y, Z, U, V, W, **kwargs)

        Arguments:

            *X*, *Y*, *Z*:
                The x, y and z coordinates of the arrow locations (default is
                tail of arrow; see *pivot* kwarg)

            *U*, *V*, *W*:
                The x, y and z components of the arrow vectors

        The arguments could be array-like or scalars, so long as they
        they can be broadcast together. The arguments can also be
        masked arrays. If an element in any of argument is masked, then
        that corresponding quiver element will not be plotted.

        Keyword arguments:

            *length*: [1.0 | float]
                The length of each quiver, default to 1.0, the unit is
                the same with the axes

            *arrow_length_ratio*: [0.3 | float]
                The ratio of the arrow head with respect to the quiver,
                default to 0.3

            *pivot*: [ 'tail' | 'middle' | 'tip' ]
                The part of the arrow that is at the grid point; the arrow
                rotates about this point, hence the name *pivot*.
                Default is 'tail'

            *normalize*: bool
                When True, all of the arrows will be the same length. This
                defaults to False, where the arrows will be different lengths
                depending on the values of u,v,w.

        Any additional keyword arguments are delegated to
        :class:`~matplotlib.collections.LineCollection`

        """
        def calc_arrow(uvw, angle=15):
            """
            To calculate the arrow head. uvw should be a unit vector.
            We normalize it here:
            """
            # get unit direction vector perpendicular to (u,v,w)
            norm = np.linalg.norm(uvw[:2])
            if norm > 0:
                x = uvw[1] / norm
                y = -uvw[0] / norm
            else:
                x, y = 0, 1

            # compute the two arrowhead direction unit vectors
            ra = math.radians(angle)
            c = math.cos(ra)
            s = math.sin(ra)

            # construct the rotation matrices
            Rpos = np.array([[c+(x**2)*(1-c), x*y*(1-c), y*s],
                             [y*x*(1-c), c+(y**2)*(1-c), -x*s],
                             [-y*s, x*s, c]])
            # opposite rotation negates all the sin terms
            Rneg = Rpos.copy()
            Rneg[[0,1,2,2],[2,2,0,1]] = -Rneg[[0,1,2,2],[2,2,0,1]]

            # multiply them to get the rotated vector
            return Rpos.dot(uvw), Rneg.dot(uvw)

        had_data = self.has_data()

        # handle kwargs
        # shaft length
        length = kwargs.pop('length', 1)
        # arrow length ratio to the shaft length
        arrow_length_ratio = kwargs.pop('arrow_length_ratio', 0.3)
        # pivot point
        pivot = kwargs.pop('pivot', 'tail')
        # normalize
        normalize = kwargs.pop('normalize', False)

        # handle args
        argi = 6
        if len(args) < argi:
            raise ValueError('Wrong number of arguments. Expected %d got %d' %
                             (argi, len(args)))

        # first 6 arguments are X, Y, Z, U, V, W
        input_args = args[:argi]
        # if any of the args are scalar, convert into list
        input_args = [[k] if isinstance(k, (int, float)) else k
                      for k in input_args]

        # extract the masks, if any
        masks = [k.mask for k in input_args if isinstance(k, np.ma.MaskedArray)]
        # broadcast to match the shape
        bcast = np.broadcast_arrays(*(input_args + masks))
        input_args = bcast[:argi]
        masks = bcast[argi:]
        if masks:
            # combine the masks into one
            mask = reduce(np.logical_or, masks)
            # put mask on and compress
            input_args = [np.ma.array(k, mask=mask).compressed()
                          for k in input_args]
        else:
            input_args = [k.flatten() for k in input_args]

        if any(len(v) == 0 for v in input_args):
            # No quivers, so just make an empty collection and return early
            linec = art3d.Line3DCollection([], *args[argi:], **kwargs)
            self.add_collection(linec)
            return linec

        # Following assertions must be true before proceeding
        # must all be ndarray
        assert all(isinstance(k, np.ndarray) for k in input_args)
        # must all in same shape
        assert len({k.shape for k in input_args}) == 1

        shaft_dt = np.linspace(0, length, num=2)
        arrow_dt = shaft_dt * arrow_length_ratio

        if pivot == 'tail':
            shaft_dt -= length
        elif pivot == 'middle':
            shaft_dt -= length/2.
        elif pivot != 'tip':
            raise ValueError('Invalid pivot argument: ' + str(pivot))

        XYZ = np.column_stack(input_args[:3])
        UVW = np.column_stack(input_args[3:argi]).astype(float)

        # Normalize rows of UVW
        # Note: with numpy 1.9+, could use np.linalg.norm(UVW, axis=1)
        norm = np.sqrt(np.sum(UVW**2, axis=1))

        # If any row of UVW is all zeros, don't make a quiver for it
        mask = norm > 0
        XYZ = XYZ[mask]
        if normalize:
            UVW = UVW[mask] / norm[mask].reshape((-1, 1))
        else:
            UVW = UVW[mask]

        if len(XYZ) > 0:
            # compute the shaft lines all at once with an outer product
            shafts = (XYZ - np.multiply.outer(shaft_dt, UVW)).swapaxes(0, 1)
            # compute head direction vectors, n heads by 2 sides by 3 dimensions
            head_dirs = np.array([calc_arrow(d) for d in UVW])
            # compute all head lines at once, starting from where the shaft ends
            heads = shafts[:, :1] - np.multiply.outer(arrow_dt, head_dirs)
            # stack left and right head lines together
            heads.shape = (len(arrow_dt), -1, 3)
            # transpose to get a list of lines
            heads = heads.swapaxes(0, 1)

            lines = list(shafts) + list(heads)
        else:
            lines = []

        linec = art3d.Line3DCollection(lines, *args[argi:], **kwargs)
        self.add_collection(linec)

        self.auto_scale_xyz(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], had_data)

        return linec

    quiver3D = quiver

    def voxels(self, *args, **kwargs):
        """
        ax.voxels([x, y, z,] /, filled, **kwargs)

        Plot a set of filled voxels

        All voxels are plotted as 1x1x1 cubes on the axis, with filled[0,0,0]
        placed with its lower corner at the origin. Occluded faces are not
        plotted.

        Call signatures::

            voxels(filled, facecolors=fc, edgecolors=ec, **kwargs)
            voxels(x, y, z, filled, facecolors=fc, edgecolors=ec, **kwargs)

        .. versionadded:: 2.1

        Parameters
        ----------
        filled : 3D np.array of bool
            A 3d array of values, with truthy values indicating which voxels
            to fill

        x, y, z : 3D np.array, optional
            The coordinates of the corners of the voxels. This should broadcast
            to a shape one larger in every dimension than the shape of `filled`.
            These can be used to plot non-cubic voxels.

            If not specified, defaults to increasing integers along each axis,
            like those returned by :func:`~numpy.indices`.
            As indicated by the ``/`` in the function signature, these arguments
            can only be passed positionally.

        facecolors, edgecolors : array_like, optional
            The color to draw the faces and edges of the voxels. Can only be
            passed as keyword arguments.
            This parameter can be:

              - A single color value, to color all voxels the same color. This
                can be either a string, or a 1D rgb/rgba array
              - ``None``, the default, to use a single color for the faces, and
                the style default for the edges.
              - A 3D ndarray of color names, with each item the color for the
                corresponding voxel. The size must match the voxels.
              - A 4D ndarray of rgb/rgba data, with the components along the
                last axis.

        **kwargs
            Additional keyword arguments to pass onto
            :func:`~mpl_toolkits.mplot3d.art3d.Poly3DCollection`

        Returns
        -------
        faces : dict
            A dictionary indexed by coordinate, where ``faces[i,j,k]`` is a
            `Poly3DCollection` of the faces drawn for the voxel
            ``filled[i,j,k]``. If no faces were drawn for a given voxel, either
            because it was not asked to be drawn, or it is fully occluded, then
            ``(i,j,k) not in faces``.

        Examples
        --------
        .. plot:: gallery/mplot3d/voxels.py
        .. plot:: gallery/mplot3d/voxels_rgb.py
        .. plot:: gallery/mplot3d/voxels_torus.py
        .. plot:: gallery/mplot3d/voxels_numpy_logo.py
        """

        # work out which signature we should be using, and use it to parse
        # the arguments. Name must be voxels for the correct error message
        if len(args) >= 3:
            # underscores indicate position only
            def voxels(__x, __y, __z, filled, **kwargs):
                return (__x, __y, __z), filled, kwargs
        else:
            def voxels(filled, **kwargs):
                return None, filled, kwargs

        xyz, filled, kwargs = voxels(*args, **kwargs)

        # check dimensions
        if filled.ndim != 3:
            raise ValueError("Argument filled must be 3-dimensional")
        size = np.array(filled.shape, dtype=np.intp)

        # check xyz coordinates, which are one larger than the filled shape
        coord_shape = tuple(size + 1)
        if xyz is None:
            x, y, z = np.indices(coord_shape)
        else:
            x, y, z = (_backports.broadcast_to(c, coord_shape) for c in xyz)

        def _broadcast_color_arg(color, name):
            if np.ndim(color) in (0, 1):
                # single color, like "red" or [1, 0, 0]
                return _backports.broadcast_to(
                    color, filled.shape + np.shape(color))
            elif np.ndim(color) in (3, 4):
                # 3D array of strings, or 4D array with last axis rgb
                if np.shape(color)[:3] != filled.shape:
                    raise ValueError(
                        "When multidimensional, {} must match the shape of "
                        "filled".format(name))
                return color
            else:
                raise ValueError("Invalid {} argument".format(name))

        # intercept the facecolors, handling defaults and broacasting
        facecolors = kwargs.pop('facecolors', None)
        if facecolors is None:
            facecolors = self._get_patches_for_fill.get_next_color()
        facecolors = _broadcast_color_arg(facecolors, 'facecolors')

        # broadcast but no default on edgecolors
        edgecolors = kwargs.pop('edgecolors', None)
        edgecolors = _broadcast_color_arg(edgecolors, 'edgecolors')

        # always scale to the full array, even if the data is only in the center
        self.auto_scale_xyz(x, y, z)

        # points lying on corners of a square
        square = np.array([
            [0, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [1, 0, 0]
        ], dtype=np.intp)

        voxel_faces = defaultdict(list)

        def permutation_matrices(n):
            """ Generator of cyclic permutation matices """
            mat = np.eye(n, dtype=np.intp)
            for i in range(n):
                yield mat
                mat = np.roll(mat, 1, axis=0)

        # iterate over each of the YZ, ZX, and XY orientations, finding faces to
        # render
        for permute in permutation_matrices(3):
            # find the set of ranges to iterate over
            pc, qc, rc = permute.T.dot(size)
            pinds = np.arange(pc)
            qinds = np.arange(qc)
            rinds = np.arange(rc)

            square_rot = square.dot(permute.T)

            # iterate within the current plane
            for p in pinds:
                for q in qinds:
                    # iterate perpendicularly to the current plane, handling
                    # boundaries. We only draw faces between a voxel and an
                    # empty space, to avoid drawing internal faces.

                    # draw lower faces
                    p0 = permute.dot([p, q, 0])
                    i0 = tuple(p0)
                    if filled[i0]:
                        voxel_faces[i0].append(p0 + square_rot)

                    # draw middle faces
                    for r1, r2 in zip(rinds[:-1], rinds[1:]):
                        p1 = permute.dot([p, q, r1])
                        p2 = permute.dot([p, q, r2])

                        i1 = tuple(p1)
                        i2 = tuple(p2)

                        if filled[i1] and not filled[i2]:
                            voxel_faces[i1].append(p2 + square_rot)
                        elif not filled[i1] and filled[i2]:
                            voxel_faces[i2].append(p2 + square_rot)

                    # draw upper faces
                    pk = permute.dot([p, q, rc-1])
                    pk2 = permute.dot([p, q, rc])
                    ik = tuple(pk)
                    if filled[ik]:
                        voxel_faces[ik].append(pk2 + square_rot)

        # iterate over the faces, and generate a Poly3DCollection for each voxel
        polygons = {}
        for coord, faces_inds in voxel_faces.items():
            # convert indices into 3D positions
            if xyz is None:
                faces = faces_inds
            else:
                faces = []
                for face_inds in faces_inds:
                    ind = face_inds[:, 0], face_inds[:, 1], face_inds[:, 2]
                    face = np.empty(face_inds.shape)
                    face[:, 0] = x[ind]
                    face[:, 1] = y[ind]
                    face[:, 2] = z[ind]
                    faces.append(face)

            poly = art3d.Poly3DCollection(faces,
                facecolors=facecolors[coord],
                edgecolors=edgecolors[coord],
                **kwargs
            )
            self.add_collection3d(poly)
            polygons[coord] = poly

        return polygons


def get_test_data(delta=0.05):
    '''
    Return a tuple X, Y, Z with a test data set.
    '''
    x = y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)

    Z1 = np.exp(-(X**2 + Y**2) / 2) / (2 * np.pi)
    Z2 = (np.exp(-(((X - 1) / 1.5)**2 + ((Y - 1) / 0.5)**2) / 2) /
          (2 * np.pi * 0.5 * 1.5))
    Z = Z2 - Z1

    X = X * 10
    Y = Y * 10
    Z = Z * 500
    return X, Y, Z


########################################################
# Register Axes3D as a 'projection' object available
# for use just like any other axes
########################################################
import matplotlib.projections as proj
proj.projection_registry.register(Axes3D)
