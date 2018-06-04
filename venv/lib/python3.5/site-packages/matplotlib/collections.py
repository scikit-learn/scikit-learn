"""
Classes for the efficient drawing of large collections of objects that
share most properties, e.g., a large number of line segments or
polygons.

The classes are not meant to be as flexible as their single element
counterparts (e.g., you may not be able to select all line styles) but
they are meant to be fast for common use cases (e.g., a large set of solid
line segemnts)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import warnings

import six
from six.moves import zip
try:
    from math import gcd
except ImportError:
    # LPy workaround
    from fractions import gcd

import numpy as np
import matplotlib as mpl
from . import (_path, artist, cbook, cm, colors as mcolors, docstring,
               lines as mlines, path as mpath, transforms)

CIRCLE_AREA_FACTOR = 1.0 / np.sqrt(np.pi)


_color_aliases = {'facecolors': ['facecolor'],
                  'edgecolors': ['edgecolor']}


class Collection(artist.Artist, cm.ScalarMappable):
    """
    Base class for Collections.  Must be subclassed to be usable.

    All properties in a collection must be sequences or scalars;
    if scalars, they will be converted to sequences.  The
    property of the ith element of the collection is::

      prop[i % len(props)]

    Exceptions are *capstyle* and *joinstyle* properties, these can
    only be set globally for the whole collection.

    Keyword arguments and default values:

        * *edgecolors*: None
        * *facecolors*: None
        * *linewidths*: None
        * *capstyle*:   None
        * *joinstyle*:  None
        * *antialiaseds*: None
        * *offsets*: None
        * *transOffset*: transforms.IdentityTransform()
        * *offset_position*: 'screen' (default) or 'data'
        * *norm*: None (optional for
          :class:`matplotlib.cm.ScalarMappable`)
        * *cmap*: None (optional for
          :class:`matplotlib.cm.ScalarMappable`)
        * *hatch*: None
        * *zorder*: 1


    *offsets* and *transOffset* are used to translate the patch after
    rendering (default no offsets).  If offset_position is 'screen'
    (default) the offset is applied after the master transform has
    been applied, that is, the offsets are in screen coordinates.  If
    offset_position is 'data', the offset is applied before the master
    transform, i.e., the offsets are in data coordinates.

    If any of *edgecolors*, *facecolors*, *linewidths*, *antialiaseds*
    are None, they default to their :data:`matplotlib.rcParams` patch
    setting, in sequence form.

    The use of :class:`~matplotlib.cm.ScalarMappable` is optional.  If
    the :class:`~matplotlib.cm.ScalarMappable` matrix _A is not None
    (i.e., a call to set_array has been made), at draw time a call to
    scalar mappable will be made to set the face colors.
    """
    _offsets = np.zeros((0, 2))
    _transOffset = transforms.IdentityTransform()
    #: Either a list of 3x3 arrays or an Nx3x3 array of transforms, suitable
    #: for the `all_transforms` argument to
    #: :meth:`~matplotlib.backend_bases.RendererBase.draw_path_collection`;
    #: each 3x3 array is used to initialize an
    #: :class:`~matplotlib.transforms.Affine2D` object.
    #: Each kind of collection defines this based on its arguments.
    _transforms = np.empty((0, 3, 3))

    # Whether to draw an edge by default.  Set on a
    # subclass-by-subclass basis.
    _edge_default = False

    def __init__(self,
                 edgecolors=None,
                 facecolors=None,
                 linewidths=None,
                 linestyles='solid',
                 capstyle=None,
                 joinstyle=None,
                 antialiaseds=None,
                 offsets=None,
                 transOffset=None,
                 norm=None,  # optional for ScalarMappable
                 cmap=None,  # ditto
                 pickradius=5.0,
                 hatch=None,
                 urls=None,
                 offset_position='screen',
                 zorder=1,
                 **kwargs
                 ):
        """
        Create a Collection

        %(Collection)s
        """
        artist.Artist.__init__(self)
        cm.ScalarMappable.__init__(self, norm, cmap)
        # list of un-scaled dash patterns
        # this is needed scaling the dash pattern by linewidth
        self._us_linestyles = [(None, None)]
        # list of dash patterns
        self._linestyles = [(None, None)]
        # list of unbroadcast/scaled linewidths
        self._us_lw = [0]
        self._linewidths = [0]
        self._is_filled = True  # May be modified by set_facecolor().

        self._hatch_color = mcolors.to_rgba(mpl.rcParams['hatch.color'])
        self.set_facecolor(facecolors)
        self.set_edgecolor(edgecolors)
        self.set_linewidth(linewidths)
        self.set_linestyle(linestyles)
        self.set_antialiased(antialiaseds)
        self.set_pickradius(pickradius)
        self.set_urls(urls)
        self.set_hatch(hatch)
        self.set_offset_position(offset_position)
        self.set_zorder(zorder)

        if capstyle:
            self.set_capstyle(capstyle)
        else:
            self._capstyle = None

        if joinstyle:
            self.set_joinstyle(joinstyle)
        else:
            self._joinstyle = None

        self._offsets = np.zeros((1, 2))
        self._uniform_offsets = None
        if offsets is not None:
            offsets = np.asanyarray(offsets, float)
            # Broadcast (2,) -> (1, 2) but nothing else.
            if offsets.shape == (2,):
                offsets = offsets[None, :]
            if transOffset is not None:
                self._offsets = offsets
                self._transOffset = transOffset
            else:
                self._uniform_offsets = offsets

        self._path_effects = None
        self.update(kwargs)
        self._paths = None

    def get_paths(self):
        return self._paths

    def set_paths(self):
        raise NotImplementedError

    def get_transforms(self):
        return self._transforms

    def get_offset_transform(self):
        t = self._transOffset
        if (not isinstance(t, transforms.Transform)
                and hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)
        return t

    def get_datalim(self, transData):
        transform = self.get_transform()
        transOffset = self.get_offset_transform()
        offsets = self._offsets
        paths = self.get_paths()

        if not transform.is_affine:
            paths = [transform.transform_path_non_affine(p) for p in paths]
            transform = transform.get_affine()
        if not transOffset.is_affine:
            offsets = transOffset.transform_non_affine(offsets)
            transOffset = transOffset.get_affine()

        if isinstance(offsets, np.ma.MaskedArray):
            offsets = offsets.filled(np.nan)
            # get_path_collection_extents handles nan but not masked arrays

        if len(paths) and len(offsets):
            result = mpath.get_path_collection_extents(
                transform.frozen(), paths, self.get_transforms(),
                offsets, transOffset.frozen())
            result = result.inverse_transformed(transData)
        else:
            result = transforms.Bbox.null()
        return result

    def get_window_extent(self, renderer):
        # TODO:check to ensure that this does not fail for
        # cases other than scatter plot legend
        return self.get_datalim(transforms.IdentityTransform())

    def _prepare_points(self):
        """Point prep for drawing and hit testing"""

        transform = self.get_transform()
        transOffset = self.get_offset_transform()
        offsets = self._offsets
        paths = self.get_paths()

        if self.have_units():
            paths = []
            for path in self.get_paths():
                vertices = path.vertices
                xs, ys = vertices[:, 0], vertices[:, 1]
                xs = self.convert_xunits(xs)
                ys = self.convert_yunits(ys)
                paths.append(mpath.Path(np.column_stack([xs, ys]), path.codes))

            if offsets.size > 0:
                xs = self.convert_xunits(offsets[:, 0])
                ys = self.convert_yunits(offsets[:, 1])
                offsets = np.column_stack([xs, ys])

        if not transform.is_affine:
            paths = [transform.transform_path_non_affine(path)
                     for path in paths]
            transform = transform.get_affine()
        if not transOffset.is_affine:
            offsets = transOffset.transform_non_affine(offsets)
            # This might have changed an ndarray into a masked array.
            transOffset = transOffset.get_affine()

        if isinstance(offsets, np.ma.MaskedArray):
            offsets = offsets.filled(np.nan)
            # Changing from a masked array to nan-filled ndarray
            # is probably most efficient at this point.

        return transform, transOffset, offsets, paths

    @artist.allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        renderer.open_group(self.__class__.__name__, self.get_gid())

        self.update_scalarmappable()

        transform, transOffset, offsets, paths = self._prepare_points()

        gc = renderer.new_gc()
        self._set_gc_clip(gc)
        gc.set_snap(self.get_snap())

        if self._hatch:
            gc.set_hatch(self._hatch)
            try:
                gc.set_hatch_color(self._hatch_color)
            except AttributeError:
                # if we end up with a GC that does not have this method
                warnings.warn("Your backend does not support setting the "
                              "hatch color.")

        if self.get_sketch_params() is not None:
            gc.set_sketch_params(*self.get_sketch_params())

        if self.get_path_effects():
            from matplotlib.patheffects import PathEffectRenderer
            renderer = PathEffectRenderer(self.get_path_effects(), renderer)

        # If the collection is made up of a single shape/color/stroke,
        # it can be rendered once and blitted multiple times, using
        # `draw_markers` rather than `draw_path_collection`.  This is
        # *much* faster for Agg, and results in smaller file sizes in
        # PDF/SVG/PS.

        trans = self.get_transforms()
        facecolors = self.get_facecolor()
        edgecolors = self.get_edgecolor()
        do_single_path_optimization = False
        if (len(paths) == 1 and len(trans) <= 1 and
            len(facecolors) == 1 and len(edgecolors) == 1 and
            len(self._linewidths) == 1 and
            self._linestyles == [(None, None)] and
            len(self._antialiaseds) == 1 and len(self._urls) == 1 and
            self.get_hatch() is None):
            if len(trans):
                combined_transform = (transforms.Affine2D(trans[0]) +
                                      transform)
            else:
                combined_transform = transform
            extents = paths[0].get_extents(combined_transform)
            width, height = renderer.get_canvas_width_height()
            if (extents.width < width and
                extents.height < height):
                do_single_path_optimization = True

        if self._joinstyle:
            gc.set_joinstyle(self._joinstyle)

        if self._capstyle:
            gc.set_capstyle(self._capstyle)

        if do_single_path_optimization:
            gc.set_foreground(tuple(edgecolors[0]))
            gc.set_linewidth(self._linewidths[0])
            gc.set_dashes(*self._linestyles[0])
            gc.set_antialiased(self._antialiaseds[0])
            gc.set_url(self._urls[0])
            renderer.draw_markers(
                gc, paths[0], combined_transform.frozen(),
                mpath.Path(offsets), transOffset, tuple(facecolors[0]))
        else:
            renderer.draw_path_collection(
                gc, transform.frozen(), paths,
                self.get_transforms(), offsets, transOffset,
                self.get_facecolor(), self.get_edgecolor(),
                self._linewidths, self._linestyles,
                self._antialiaseds, self._urls,
                self._offset_position)

        gc.restore()
        renderer.close_group(self.__class__.__name__)
        self.stale = False

    def set_pickradius(self, pr):
        """Set the pick radius used for containment tests.

        .. ACCEPTS: float distance in points

        Parameters
        ----------
        d : float
            Pick radius, in points.
        """
        self._pickradius = pr

    def get_pickradius(self):
        return self._pickradius

    def contains(self, mouseevent):
        """
        Test whether the mouse event occurred in the collection.

        Returns True | False, ``dict(ind=itemlist)``, where every
        item in itemlist contains the event.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)

        if not self.get_visible():
            return False, {}

        pickradius = (
            float(self._picker)
            if cbook.is_numlike(self._picker) and
               self._picker is not True  # the bool, not just nonzero or 1
            else self._pickradius)

        transform, transOffset, offsets, paths = self._prepare_points()

        ind = _path.point_in_path_collection(
            mouseevent.x, mouseevent.y, pickradius,
            transform.frozen(), paths, self.get_transforms(),
            offsets, transOffset, pickradius <= 0,
            self.get_offset_position())

        return len(ind) > 0, dict(ind=ind)

    def set_urls(self, urls):
        """
        Parameters
        ----------
        urls : List[str] or None
            .. ACCEPTS: List[str] or None
        """
        self._urls = urls if urls is not None else [None]
        self.stale = True

    def get_urls(self):
        return self._urls

    def set_hatch(self, hatch):
        r"""
        Set the hatching pattern

        *hatch* can be one of::

          /   - diagonal hatching
          \   - back diagonal
          |   - vertical
          -   - horizontal
          +   - crossed
          x   - crossed diagonal
          o   - small circle
          O   - large circle
          .   - dots
          *   - stars

        Letters can be combined, in which case all the specified
        hatchings are done.  If same letter repeats, it increases the
        density of hatching of that pattern.

        Hatching is supported in the PostScript, PDF, SVG and Agg
        backends only.

        Unlike other properties such as linewidth and colors, hatching
        can only be specified for the collection as a whole, not separately
        for each member.

        ACCEPTS: [ '/' | '\\' | '|' | '-' | '+' | 'x' | 'o' | 'O' | '.' | '*' ]
        """
        self._hatch = hatch
        self.stale = True

    def get_hatch(self):
        """Return the current hatching pattern."""
        return self._hatch

    def set_offsets(self, offsets):
        """
        Set the offsets for the collection.  *offsets* can be a scalar
        or a sequence.

        ACCEPTS: float or sequence of floats
        """
        offsets = np.asanyarray(offsets, float)
        if offsets.shape == (2,):  # Broadcast (2,) -> (1, 2) but nothing else.
            offsets = offsets[None, :]
        # This decision is based on how they are initialized above in __init__.
        if self._uniform_offsets is None:
            self._offsets = offsets
        else:
            self._uniform_offsets = offsets
        self.stale = True

    def get_offsets(self):
        """Return the offsets for the collection."""
        # This decision is based on how they are initialized above in __init__.
        if self._uniform_offsets is None:
            return self._offsets
        else:
            return self._uniform_offsets

    def set_offset_position(self, offset_position):
        """
        Set how offsets are applied.  If *offset_position* is 'screen'
        (default) the offset is applied after the master transform has
        been applied, that is, the offsets are in screen coordinates.
        If offset_position is 'data', the offset is applied before the
        master transform, i.e., the offsets are in data coordinates.

        .. ACCEPTS: [ 'screen' | 'data' ]
        """
        if offset_position not in ('screen', 'data'):
            raise ValueError("offset_position must be 'screen' or 'data'")
        self._offset_position = offset_position
        self.stale = True

    def get_offset_position(self):
        """
        Returns how offsets are applied for the collection.  If
        *offset_position* is 'screen', the offset is applied after the
        master transform has been applied, that is, the offsets are in
        screen coordinates.  If offset_position is 'data', the offset
        is applied before the master transform, i.e., the offsets are
        in data coordinates.
        """
        return self._offset_position

    def set_linewidth(self, lw):
        """
        Set the linewidth(s) for the collection.  *lw* can be a scalar
        or a sequence; if it is a sequence the patches will cycle
        through the sequence

        ACCEPTS: float or sequence of floats
        """
        if lw is None:
            lw = mpl.rcParams['patch.linewidth']
            if lw is None:
                lw = mpl.rcParams['lines.linewidth']
        # get the un-scaled/broadcast lw
        self._us_lw = np.atleast_1d(np.asarray(lw))

        # scale all of the dash patterns.
        self._linewidths, self._linestyles = self._bcast_lwls(
            self._us_lw, self._us_linestyles)
        self.stale = True

    def set_linewidths(self, lw):
        """alias for set_linewidth"""
        return self.set_linewidth(lw)

    def set_lw(self, lw):
        """alias for set_linewidth"""
        return self.set_linewidth(lw)

    def set_linestyle(self, ls):
        """
        Set the linestyle(s) for the collection.

        ===========================   =================
        linestyle                     description
        ===========================   =================
        ``'-'`` or ``'solid'``        solid line
        ``'--'`` or  ``'dashed'``     dashed line
        ``'-.'`` or  ``'dashdot'``    dash-dotted line
        ``':'`` or ``'dotted'``       dotted line
        ===========================   =================

        Alternatively a dash tuple of the following form can be provided::

            (offset, onoffseq),

        where ``onoffseq`` is an even length tuple of on and off ink
        in points.

        ACCEPTS: ['solid' | 'dashed', 'dashdot', 'dotted' |
                   (offset, on-off-dash-seq) |
                   ``'-'`` | ``'--'`` | ``'-.'`` | ``':'`` | ``'None'`` |
                   ``' '`` | ``''``]

        Parameters
        ----------
        ls : { '-',  '--', '-.', ':'} and more see description
            The line style.
        """
        try:
            if isinstance(ls, six.string_types):
                ls = cbook.ls_mapper.get(ls, ls)
                dashes = [mlines._get_dash_pattern(ls)]
            else:
                try:
                    dashes = [mlines._get_dash_pattern(ls)]
                except ValueError:
                    dashes = [mlines._get_dash_pattern(x) for x in ls]

        except ValueError:
            raise ValueError(
                'Do not know how to convert {!r} to dashes'.format(ls))

        # get the list of raw 'unscaled' dash patterns
        self._us_linestyles = dashes

        # broadcast and scale the lw and dash patterns
        self._linewidths, self._linestyles = self._bcast_lwls(
            self._us_lw, self._us_linestyles)

    def set_capstyle(self, cs):
        """
        Set the capstyle for the collection. The capstyle can
        only be set globally for all elements in the collection

        Parameters
        ----------
        cs : ['butt' | 'round' | 'projecting']
            The capstyle
        """
        if cs in ('butt', 'round', 'projecting'):
            self._capstyle = cs
        else:
            raise ValueError('Unrecognized cap style.  Found %s' % cs)

    def get_capstyle(self):
        return self._capstyle

    def set_joinstyle(self, js):
        """
        Set the joinstyle for the collection. The joinstyle can only be
        set globally for all elements in the collection.

        Parameters
        ----------
        js : ['miter' | 'round' | 'bevel']
            The joinstyle
        """
        if js in ('miter', 'round', 'bevel'):
            self._joinstyle = js
        else:
            raise ValueError('Unrecognized join style.  Found %s' % js)

    def get_joinstyle(self):
        return self._joinstyle

    @staticmethod
    def _bcast_lwls(linewidths, dashes):
        '''Internal helper function to broadcast + scale ls/lw

        In the collection drawing code the linewidth and linestyle are
        cycled through as circular buffers (via v[i % len(v)]).  Thus,
        if we are going to scale the dash pattern at set time (not
        draw time) we need to do the broadcasting now and expand both
        lists to be the same length.

        Parameters
        ----------
        linewidths : list
            line widths of collection

        dashes : list
            dash specification (offset, (dash pattern tuple))

        Returns
        -------
        linewidths, dashes : list
             Will be the same length, dashes are scaled by paired linewidth

        '''
        if mpl.rcParams['_internal.classic_mode']:
            return linewidths, dashes
        # make sure they are the same length so we can zip them
        if len(dashes) != len(linewidths):
            l_dashes = len(dashes)
            l_lw = len(linewidths)
            GCD = gcd(l_dashes, l_lw)
            dashes = list(dashes) * (l_lw // GCD)
            linewidths = list(linewidths) * (l_dashes // GCD)

        # scale the dash patters
        dashes = [mlines._scale_dashes(o, d, lw)
                  for (o, d), lw in zip(dashes, linewidths)]

        return linewidths, dashes

    def set_linestyles(self, ls):
        """alias for set_linestyle"""
        return self.set_linestyle(ls)

    def set_dashes(self, ls):
        """alias for set_linestyle"""
        return self.set_linestyle(ls)

    def set_antialiased(self, aa):
        """
        Set the antialiasing state for rendering.

        ACCEPTS: Boolean or sequence of booleans
        """
        if aa is None:
            aa = mpl.rcParams['patch.antialiased']
        self._antialiaseds = np.atleast_1d(np.asarray(aa, bool))
        self.stale = True

    def set_antialiaseds(self, aa):
        """alias for set_antialiased"""
        return self.set_antialiased(aa)

    def set_color(self, c):
        """
        Set both the edgecolor and the facecolor.

        ACCEPTS: matplotlib color arg or sequence of rgba tuples

        .. seealso::

            :meth:`set_facecolor`, :meth:`set_edgecolor`
               For setting the edge or face color individually.
        """
        self.set_facecolor(c)
        self.set_edgecolor(c)

    def _set_facecolor(self, c):
        if c is None:
            c = mpl.rcParams['patch.facecolor']

        self._is_filled = True
        try:
            if c.lower() == 'none':
                self._is_filled = False
        except AttributeError:
            pass
        self._facecolors = mcolors.to_rgba_array(c, self._alpha)
        self.stale = True

    def set_facecolor(self, c):
        """
        Set the facecolor(s) of the collection.  *c* can be a
        matplotlib color spec (all patches have same color), or a
        sequence of specs; if it is a sequence the patches will
        cycle through the sequence.

        If *c* is 'none', the patch will not be filled.

        ACCEPTS: matplotlib color spec or sequence of specs
        """
        self._original_facecolor = c
        self._set_facecolor(c)

    def set_facecolors(self, c):
        """alias for set_facecolor"""
        return self.set_facecolor(c)

    def get_facecolor(self):
        return self._facecolors
    get_facecolors = get_facecolor

    def get_edgecolor(self):
        if (isinstance(self._edgecolors, six.string_types)
                   and self._edgecolors == str('face')):
            return self.get_facecolors()
        else:
            return self._edgecolors
    get_edgecolors = get_edgecolor

    def _set_edgecolor(self, c):
        set_hatch_color = True
        if c is None:
            if (mpl.rcParams['patch.force_edgecolor'] or
                    not self._is_filled or self._edge_default):
                c = mpl.rcParams['patch.edgecolor']
            else:
                c = 'none'
                set_hatch_color = False

        self._is_stroked = True
        try:
            if c.lower() == 'none':
                self._is_stroked = False
        except AttributeError:
            pass

        try:
            if c.lower() == 'face':   # Special case: lookup in "get" method.
                self._edgecolors = 'face'
                return
        except AttributeError:
            pass
        self._edgecolors = mcolors.to_rgba_array(c, self._alpha)
        if set_hatch_color and len(self._edgecolors):
            self._hatch_color = tuple(self._edgecolors[0])
        self.stale = True

    def set_edgecolor(self, c):
        """
        Set the edgecolor(s) of the collection. *c* can be a
        matplotlib color spec (all patches have same color), or a
        sequence of specs; if it is a sequence the patches will
        cycle through the sequence.

        If *c* is 'face', the edge color will always be the same as
        the face color.  If it is 'none', the patch boundary will not
        be drawn.

        ACCEPTS: matplotlib color spec or sequence of specs
        """
        self._original_edgecolor = c
        self._set_edgecolor(c)

    def set_edgecolors(self, c):
        """alias for set_edgecolor"""
        return self.set_edgecolor(c)

    def set_alpha(self, alpha):
        """
        Set the alpha tranparencies of the collection.  *alpha* must be
        a float or *None*.

        ACCEPTS: float or None
        """
        if alpha is not None:
            try:
                float(alpha)
            except TypeError:
                raise TypeError('alpha must be a float or None')
        self.update_dict['array'] = True
        artist.Artist.set_alpha(self, alpha)
        self._set_facecolor(self._original_facecolor)
        self._set_edgecolor(self._original_edgecolor)

    def get_linewidths(self):
        return self._linewidths
    get_linewidth = get_linewidths

    def get_linestyles(self):
        return self._linestyles
    get_dashes = get_linestyle = get_linestyles

    def update_scalarmappable(self):
        """
        If the scalar mappable array is not none, update colors
        from scalar data
        """
        if self._A is None:
            return
        if self._A.ndim > 1:
            raise ValueError('Collections can only map rank 1 arrays')
        if not self.check_update("array"):
            return
        if self._is_filled:
            self._facecolors = self.to_rgba(self._A, self._alpha)
        elif self._is_stroked:
            self._edgecolors = self.to_rgba(self._A, self._alpha)
        self.stale = True

    def get_fill(self):
        'return whether fill is set'
        return self._is_filled

    def update_from(self, other):
        'copy properties from other to self'

        artist.Artist.update_from(self, other)
        self._antialiaseds = other._antialiaseds
        self._original_edgecolor = other._original_edgecolor
        self._edgecolors = other._edgecolors
        self._original_facecolor = other._original_facecolor
        self._facecolors = other._facecolors
        self._linewidths = other._linewidths
        self._linestyles = other._linestyles
        self._us_linestyles = other._us_linestyles
        self._pickradius = other._pickradius
        self._hatch = other._hatch

        # update_from for scalarmappable
        self._A = other._A
        self.norm = other.norm
        self.cmap = other.cmap
        # self.update_dict = other.update_dict # do we need to copy this? -JJL
        self.stale = True

# these are not available for the object inspector until after the
# class is built so we define an initial set here for the init
# function and they will be overridden after object defn
docstring.interpd.update(Collection="""\
    Valid Collection keyword arguments:

        * *edgecolors*: None
        * *facecolors*: None
        * *linewidths*: None
        * *antialiaseds*: None
        * *offsets*: None
        * *transOffset*: transforms.IdentityTransform()
        * *norm*: None (optional for
          :class:`matplotlib.cm.ScalarMappable`)
        * *cmap*: None (optional for
          :class:`matplotlib.cm.ScalarMappable`)

    *offsets* and *transOffset* are used to translate the patch after
    rendering (default no offsets)

    If any of *edgecolors*, *facecolors*, *linewidths*, *antialiaseds*
    are None, they default to their :data:`matplotlib.rcParams` patch
    setting, in sequence form.
""")


class _CollectionWithSizes(Collection):
    """
    Base class for collections that have an array of sizes.
    """
    _factor = 1.0

    def get_sizes(self):
        """
        Returns the sizes of the elements in the collection.  The
        value represents the 'area' of the element.

        Returns
        -------
        sizes : array
            The 'area' of each element.
        """
        return self._sizes

    def set_sizes(self, sizes, dpi=72.0):
        """
        Set the sizes of each member of the collection.

        Parameters
        ----------
        sizes : ndarray or None
            The size to set for each element of the collection.  The
            value is the 'area' of the element.

        dpi : float
            The dpi of the canvas. Defaults to 72.0.
        """
        if sizes is None:
            self._sizes = np.array([])
            self._transforms = np.empty((0, 3, 3))
        else:
            self._sizes = np.asarray(sizes)
            self._transforms = np.zeros((len(self._sizes), 3, 3))
            scale = np.sqrt(self._sizes) * dpi / 72.0 * self._factor
            self._transforms[:, 0, 0] = scale
            self._transforms[:, 1, 1] = scale
            self._transforms[:, 2, 2] = 1.0
        self.stale = True

    @artist.allow_rasterization
    def draw(self, renderer):
        self.set_sizes(self._sizes, self.figure.dpi)
        Collection.draw(self, renderer)


class PathCollection(_CollectionWithSizes):
    """
    This is the most basic :class:`Collection` subclass.
    """
    @docstring.dedent_interpd
    def __init__(self, paths, sizes=None, **kwargs):
        """
        *paths* is a sequence of :class:`matplotlib.path.Path`
        instances.

        %(Collection)s
        """

        Collection.__init__(self, **kwargs)
        self.set_paths(paths)
        self.set_sizes(sizes)
        self.stale = True

    def set_paths(self, paths):
        self._paths = paths
        self.stale = True

    def get_paths(self):
        return self._paths


class PolyCollection(_CollectionWithSizes):
    @docstring.dedent_interpd
    def __init__(self, verts, sizes=None, closed=True, **kwargs):
        """
        *verts* is a sequence of ( *verts0*, *verts1*, ...) where
        *verts_i* is a sequence of *xy* tuples of vertices, or an
        equivalent :mod:`numpy` array of shape (*nv*, 2).

        *sizes* is *None* (default) or a sequence of floats that
        scale the corresponding *verts_i*.  The scaling is applied
        before the Artist master transform; if the latter is an identity
        transform, then the overall scaling is such that if
        *verts_i* specify a unit square, then *sizes_i* is the area
        of that square in points^2.
        If len(*sizes*) < *nv*, the additional values will be
        taken cyclically from the array.

        *closed*, when *True*, will explicitly close the polygon.

        %(Collection)s
        """
        Collection.__init__(self, **kwargs)
        self.set_sizes(sizes)
        self.set_verts(verts, closed)
        self.stale = True

    def set_verts(self, verts, closed=True):
        '''This allows one to delay initialization of the vertices.'''
        if isinstance(verts, np.ma.MaskedArray):
            verts = verts.astype(float).filled(np.nan)
            # This is much faster than having Path do it one at a time.
        if closed:
            self._paths = []
            for xy in verts:
                if len(xy):
                    if isinstance(xy, np.ma.MaskedArray):
                        xy = np.ma.concatenate([xy, xy[0:1]])
                    else:
                        xy = np.asarray(xy)
                        xy = np.concatenate([xy, xy[0:1]])
                    codes = np.empty(xy.shape[0], dtype=mpath.Path.code_type)
                    codes[:] = mpath.Path.LINETO
                    codes[0] = mpath.Path.MOVETO
                    codes[-1] = mpath.Path.CLOSEPOLY
                    self._paths.append(mpath.Path(xy, codes))
                else:
                    self._paths.append(mpath.Path(xy))
        else:
            self._paths = [mpath.Path(xy) for xy in verts]
        self.stale = True

    set_paths = set_verts

    def set_verts_and_codes(self, verts, codes):
        '''This allows one to initialize vertices with path codes.'''
        if (len(verts) != len(codes)):
            raise ValueError("'codes' must be a 1D list or array "
                             "with the same length of 'verts'")
        self._paths = []
        for xy, cds in zip(verts, codes):
            if len(xy):
                self._paths.append(mpath.Path(xy, cds))
            else:
                self._paths.append(mpath.Path(xy))
        self.stale = True


class BrokenBarHCollection(PolyCollection):
    """
    A collection of horizontal bars spanning *yrange* with a sequence of
    *xranges*.
    """
    @docstring.dedent_interpd
    def __init__(self, xranges, yrange, **kwargs):
        """
        *xranges*
            sequence of (*xmin*, *xwidth*)

        *yrange*
            *ymin*, *ywidth*

        %(Collection)s
        """
        ymin, ywidth = yrange
        ymax = ymin + ywidth
        verts = [[(xmin, ymin),
                  (xmin, ymax),
                  (xmin + xwidth, ymax),
                  (xmin + xwidth, ymin),
                  (xmin, ymin)] for xmin, xwidth in xranges]
        PolyCollection.__init__(self, verts, **kwargs)

    @staticmethod
    def span_where(x, ymin, ymax, where, **kwargs):
        """
        Create a BrokenBarHCollection to plot horizontal bars from
        over the regions in *x* where *where* is True.  The bars range
        on the y-axis from *ymin* to *ymax*

        A :class:`BrokenBarHCollection` is returned.  *kwargs* are
        passed on to the collection.
        """
        xranges = []
        for ind0, ind1 in cbook.contiguous_regions(where):
            xslice = x[ind0:ind1]
            if not len(xslice):
                continue
            xranges.append((xslice[0], xslice[-1] - xslice[0]))

        collection = BrokenBarHCollection(
            xranges, [ymin, ymax - ymin], **kwargs)
        return collection


class RegularPolyCollection(_CollectionWithSizes):
    """Draw a collection of regular polygons with *numsides*."""
    _path_generator = mpath.Path.unit_regular_polygon

    _factor = CIRCLE_AREA_FACTOR

    @docstring.dedent_interpd
    def __init__(self,
                 numsides,
                 rotation=0,
                 sizes=(1,),
                 **kwargs):
        """
        *numsides*
            the number of sides of the polygon

        *rotation*
            the rotation of the polygon in radians

        *sizes*
            gives the area of the circle circumscribing the
            regular polygon in points^2

        %(Collection)s

        Example: see :file:`examples/dynamic_collection.py` for
        complete example::

            offsets = np.random.rand(20,2)
            facecolors = [cm.jet(x) for x in np.random.rand(20)]
            black = (0,0,0,1)

            collection = RegularPolyCollection(
                numsides=5, # a pentagon
                rotation=0, sizes=(50,),
                facecolors = facecolors,
                edgecolors = (black,),
                linewidths = (1,),
                offsets = offsets,
                transOffset = ax.transData,
                )
        """
        Collection.__init__(self, **kwargs)
        self.set_sizes(sizes)
        self._numsides = numsides
        self._paths = [self._path_generator(numsides)]
        self._rotation = rotation
        self.set_transform(transforms.IdentityTransform())

    def get_numsides(self):
        return self._numsides

    def get_rotation(self):
        return self._rotation

    @artist.allow_rasterization
    def draw(self, renderer):
        self.set_sizes(self._sizes, self.figure.dpi)
        self._transforms = [
            transforms.Affine2D(x).rotate(-self._rotation).get_matrix()
            for x in self._transforms
        ]
        Collection.draw(self, renderer)


class StarPolygonCollection(RegularPolyCollection):
    """
    Draw a collection of regular stars with *numsides* points."""

    _path_generator = mpath.Path.unit_regular_star


class AsteriskPolygonCollection(RegularPolyCollection):
    """
    Draw a collection of regular asterisks with *numsides* points."""

    _path_generator = mpath.Path.unit_regular_asterisk


class LineCollection(Collection):
    """
    All parameters must be sequences or scalars; if scalars, they will
    be converted to sequences.  The property of the ith line
    segment is::

       prop[i % len(props)]

    i.e., the properties cycle if the ``len`` of props is less than the
    number of segments.
    """

    _edge_default = True

    def __init__(self, segments,     # Can be None.
                 linewidths=None,
                 colors=None,
                 antialiaseds=None,
                 linestyles='solid',
                 offsets=None,
                 transOffset=None,
                 norm=None,
                 cmap=None,
                 pickradius=5,
                 zorder=2,
                 facecolors='none',
                 **kwargs
                 ):
        """
        Parameters
        ----------
        segments :
            A sequence of (*line0*, *line1*, *line2*), where::

                linen = (x0, y0), (x1, y1), ... (xm, ym)

            or the equivalent numpy array with two columns. Each line
            can be a different length.

        colors : sequence, optional
            A sequence of RGBA tuples (e.g., arbitrary color
            strings, etc, not allowed).

        antialiaseds : sequence, optional
            A sequence of ones or zeros.

        linestyles : string, tuple, optional
            Either one of [ 'solid' | 'dashed' | 'dashdot' | 'dotted' ], or
            a dash tuple. The dash tuple is::

                (offset, onoffseq)

            where ``onoffseq`` is an even length tuple of on and off ink
            in points.

        norm : Normalize, optional
            `~.colors.Normalize` instance.

        cmap : string or Colormap, optional
            Colormap name or `~.colors.Colormap` instance.

        pickradius : float, optional
            The tolerance in points for mouse clicks picking a line.
            Default is 5 pt.

        zorder : int, optional
           zorder of the LineCollection. Default is 2.

        facecolors : optional
           The facecolors of the LineCollection. Default is 'none'.
           Setting to a value other than 'none' will lead to a filled
           polygon being drawn between points on each line.

        Notes
        -----
        If *linewidths*, *colors*, or *antialiaseds* is None, they
        default to their rcParams setting, in sequence form.

        If *offsets* and *transOffset* are not None, then
        *offsets* are transformed by *transOffset* and applied after
        the segments have been transformed to display coordinates.

        If *offsets* is not None but *transOffset* is None, then the
        *offsets* are added to the segments before any transformation.
        In this case, a single offset can be specified as::

            offsets=(xo,yo)

        and this value will be added cumulatively to each successive
        segment, so as to produce a set of successively offset curves.

        The use of :class:`~matplotlib.cm.ScalarMappable` is optional.
        If the :class:`~matplotlib.cm.ScalarMappable` array
        :attr:`~matplotlib.cm.ScalarMappable._A` is not None (i.e., a call to
        :meth:`~matplotlib.cm.ScalarMappable.set_array` has been made), at
        draw time a call to scalar mappable will be made to set the colors.
        """
        if colors is None:
            colors = mpl.rcParams['lines.color']
        if linewidths is None:
            linewidths = (mpl.rcParams['lines.linewidth'],)
        if antialiaseds is None:
            antialiaseds = (mpl.rcParams['lines.antialiased'],)

        colors = mcolors.to_rgba_array(colors)

        Collection.__init__(
            self,
            edgecolors=colors,
            facecolors=facecolors,
            linewidths=linewidths,
            linestyles=linestyles,
            antialiaseds=antialiaseds,
            offsets=offsets,
            transOffset=transOffset,
            norm=norm,
            cmap=cmap,
            pickradius=pickradius,
            zorder=zorder,
            **kwargs)

        self.set_segments(segments)

    def set_segments(self, segments):
        if segments is None:
            return
        _segments = []

        for seg in segments:
            if not isinstance(seg, np.ma.MaskedArray):
                seg = np.asarray(seg, float)
            _segments.append(seg)

        if self._uniform_offsets is not None:
            _segments = self._add_offsets(_segments)

        self._paths = [mpath.Path(_seg) for _seg in _segments]
        self.stale = True

    set_verts = set_segments  # for compatibility with PolyCollection
    set_paths = set_segments

    def get_segments(self):
        """
        Returns
        -------
        segments : list
            List of segments in the LineCollection. Each list item contains an
            array of vertices.
        """
        segments = []

        for path in self._paths:
            vertices = [vertex for vertex, _ in path.iter_segments()]
            vertices = np.asarray(vertices)
            segments.append(vertices)

        return segments

    def _add_offsets(self, segs):
        offsets = self._uniform_offsets
        Nsegs = len(segs)
        Noffs = offsets.shape[0]
        if Noffs == 1:
            for i in range(Nsegs):
                segs[i] = segs[i] + i * offsets
        else:
            for i in range(Nsegs):
                io = i % Noffs
                segs[i] = segs[i] + offsets[io:io + 1]
        return segs

    def set_color(self, c):
        """
        Set the color(s) of the LineCollection.

        Parameters
        ----------
        c :
            Matplotlib color argument (all patches have same color), or a
            sequence or rgba tuples; if it is a sequence the patches will
            cycle through the sequence.
        """
        self.set_edgecolor(c)
        self.stale = True

    def get_color(self):
        return self._edgecolors

    get_colors = get_color  # for compatibility with old versions


class EventCollection(LineCollection):
    '''
    A collection of discrete events.

    The events are given by a 1-dimensional array, usually the position of
    something along an axis, such as time or length.  They do not have an
    amplitude and are displayed as vertical or horizontal parallel bars.
    '''

    _edge_default = True

    def __init__(self,
                 positions,     # Cannot be None.
                 orientation=None,
                 lineoffset=0,
                 linelength=1,
                 linewidth=None,
                 color=None,
                 linestyle='solid',
                 antialiased=None,
                 **kwargs
                 ):
        """
        Parameters
        ----------
        positions : 1D array-like object
            Each value is an event.

        orientation : {None, 'horizontal', 'vertical'}, optional
            The orientation of the **collection** (the event bars are along
            the orthogonal direction). Defaults to 'horizontal' if not
            specified or None.

        lineoffset : scalar, optional, default: 0
            The offset of the center of the markers from the origin, in the
            direction orthogonal to *orientation*.

        linelength : scalar, optional, default: 1
            The total height of the marker (i.e. the marker stretches from
            ``lineoffset - linelength/2`` to ``lineoffset + linelength/2``).

        linewidth : scalar or None, optional, default: None
            If it is None, defaults to its rcParams setting, in sequence form.

        color : color, sequence of colors or None, optional, default: None
            If it is None, defaults to its rcParams setting, in sequence form.

        linestyle : str or tuple, optional, default: 'solid'
            Valid strings are ['solid', 'dashed', 'dashdot', 'dotted',
            '-', '--', '-.', ':']. Dash tuples should be of the form::

                (offset, onoffseq),

            where *onoffseq* is an even length tuple of on and off ink
            in points.

        antialiased : {None, 1, 2}, optional
            If it is None, defaults to its rcParams setting, in sequence form.

        **kwargs : optional
            Other keyword arguments are line collection properties.  See
            :class:`~matplotlib.collections.LineCollection` for a list of
            the valid properties.

        Examples
        --------

        .. plot:: gallery/lines_bars_and_markers/eventcollection_demo.py
        """

        segment = (lineoffset + linelength / 2.,
                   lineoffset - linelength / 2.)
        if positions is None or len(positions) == 0:
            segments = []
        elif hasattr(positions, 'ndim') and positions.ndim > 1:
            raise ValueError('positions cannot be an array with more than '
                             'one dimension.')
        elif (orientation is None or orientation.lower() == 'none' or
              orientation.lower() == 'horizontal'):
            positions.sort()
            segments = [[(coord1, coord2) for coord2 in segment] for
                        coord1 in positions]
            self._is_horizontal = True
        elif orientation.lower() == 'vertical':
            positions.sort()
            segments = [[(coord2, coord1) for coord2 in segment] for
                        coord1 in positions]
            self._is_horizontal = False
        else:
            raise ValueError("orientation must be 'horizontal' or 'vertical'")

        LineCollection.__init__(self,
                                segments,
                                linewidths=linewidth,
                                colors=color,
                                antialiaseds=antialiased,
                                linestyles=linestyle,
                                **kwargs)

        self._linelength = linelength
        self._lineoffset = lineoffset

    def get_positions(self):
        '''
        return an array containing the floating-point values of the positions
        '''
        segments = self.get_segments()
        pos = 0 if self.is_horizontal() else 1
        positions = []
        for segment in segments:
            positions.append(segment[0, pos])
        return positions

    def set_positions(self, positions):
        '''
        set the positions of the events to the specified value
        '''
        if positions is None or (hasattr(positions, 'len') and
                                 len(positions) == 0):
            self.set_segments([])
            return

        lineoffset = self.get_lineoffset()
        linelength = self.get_linelength()
        segment = (lineoffset + linelength / 2.,
                   lineoffset - linelength / 2.)
        positions = np.asanyarray(positions)
        positions.sort()
        if self.is_horizontal():
            segments = [[(coord1, coord2) for coord2 in segment] for
                        coord1 in positions]
        else:
            segments = [[(coord2, coord1) for coord2 in segment] for
                        coord1 in positions]
        self.set_segments(segments)

    def add_positions(self, position):
        '''
        add one or more events at the specified positions
        '''
        if position is None or (hasattr(position, 'len') and
                                len(position) == 0):
            return
        positions = self.get_positions()
        positions = np.hstack([positions, np.asanyarray(position)])
        self.set_positions(positions)
    extend_positions = append_positions = add_positions

    def is_horizontal(self):
        '''
        True if the eventcollection is horizontal, False if vertical
        '''
        return self._is_horizontal

    def get_orientation(self):
        '''
        get the orientation of the event line, may be:
        [ 'horizontal' | 'vertical' ]
        '''
        return 'horizontal' if self.is_horizontal() else 'vertical'

    def switch_orientation(self):
        '''
        switch the orientation of the event line, either from vertical to
        horizontal or vice versus
        '''
        segments = self.get_segments()
        for i, segment in enumerate(segments):
            segments[i] = np.fliplr(segment)
        self.set_segments(segments)
        self._is_horizontal = not self.is_horizontal()
        self.stale = True

    def set_orientation(self, orientation=None):
        '''
        set the orientation of the event line
        [ 'horizontal' | 'vertical' | None ]
        defaults to 'horizontal' if not specified or None
        '''
        if (orientation is None or orientation.lower() == 'none' or
                orientation.lower() == 'horizontal'):
            is_horizontal = True
        elif orientation.lower() == 'vertical':
            is_horizontal = False
        else:
            raise ValueError("orientation must be 'horizontal' or 'vertical'")

        if is_horizontal == self.is_horizontal():
            return
        self.switch_orientation()

    def get_linelength(self):
        '''
        get the length of the lines used to mark each event
        '''
        return self._linelength

    def set_linelength(self, linelength):
        '''
        set the length of the lines used to mark each event
        '''
        if linelength == self.get_linelength():
            return
        lineoffset = self.get_lineoffset()
        segments = self.get_segments()
        pos = 1 if self.is_horizontal() else 0
        for segment in segments:
            segment[0, pos] = lineoffset + linelength / 2.
            segment[1, pos] = lineoffset - linelength / 2.
        self.set_segments(segments)
        self._linelength = linelength

    def get_lineoffset(self):
        '''
        get the offset of the lines used to mark each event
        '''
        return self._lineoffset

    def set_lineoffset(self, lineoffset):
        '''
        set the offset of the lines used to mark each event
        '''
        if lineoffset == self.get_lineoffset():
            return
        linelength = self.get_linelength()
        segments = self.get_segments()
        pos = 1 if self.is_horizontal() else 0
        for segment in segments:
            segment[0, pos] = lineoffset + linelength / 2.
            segment[1, pos] = lineoffset - linelength / 2.
        self.set_segments(segments)
        self._lineoffset = lineoffset

    def get_linewidth(self):
        '''
        get the width of the lines used to mark each event
        '''
        return self.get_linewidths()[0]

    def get_linestyle(self):
        '''
        get the style of the lines used to mark each event
        [ 'solid' | 'dashed' | 'dashdot' | 'dotted' ]
        '''
        return self.get_linestyles()

    def get_color(self):
        '''
        get the color of the lines used to mark each event
        '''
        return self.get_colors()[0]


class CircleCollection(_CollectionWithSizes):
    """
    A collection of circles, drawn using splines.
    """
    _factor = CIRCLE_AREA_FACTOR

    @docstring.dedent_interpd
    def __init__(self, sizes, **kwargs):
        """
        *sizes*
            Gives the area of the circle in points^2

        %(Collection)s
        """
        Collection.__init__(self, **kwargs)
        self.set_sizes(sizes)
        self.set_transform(transforms.IdentityTransform())
        self._paths = [mpath.Path.unit_circle()]


class EllipseCollection(Collection):
    """
    A collection of ellipses, drawn using splines.
    """
    @docstring.dedent_interpd
    def __init__(self, widths, heights, angles, units='points', **kwargs):
        """
        *widths*: sequence
            lengths of first axes (e.g., major axis lengths)

        *heights*: sequence
            lengths of second axes

        *angles*: sequence
            angles of first axes, degrees CCW from the X-axis

        *units*: ['points' | 'inches' | 'dots' | 'width' | 'height'
        | 'x' | 'y' | 'xy']

            units in which majors and minors are given; 'width' and
            'height' refer to the dimensions of the axes, while 'x'
            and 'y' refer to the *offsets* data units. 'xy' differs
            from all others in that the angle as plotted varies with
            the aspect ratio, and equals the specified angle only when
            the aspect ratio is unity.  Hence it behaves the same as
            the :class:`~matplotlib.patches.Ellipse` with
            axes.transData as its transform.

        Additional kwargs inherited from the base :class:`Collection`:

        %(Collection)s
        """
        Collection.__init__(self, **kwargs)
        self._widths = 0.5 * np.asarray(widths).ravel()
        self._heights = 0.5 * np.asarray(heights).ravel()
        self._angles = np.deg2rad(angles).ravel()
        self._units = units
        self.set_transform(transforms.IdentityTransform())
        self._transforms = np.empty((0, 3, 3))
        self._paths = [mpath.Path.unit_circle()]

    def _set_transforms(self):
        """
        Calculate transforms immediately before drawing.
        """
        ax = self.axes
        fig = self.figure

        if self._units == 'xy':
            sc = 1
        elif self._units == 'x':
            sc = ax.bbox.width / ax.viewLim.width
        elif self._units == 'y':
            sc = ax.bbox.height / ax.viewLim.height
        elif self._units == 'inches':
            sc = fig.dpi
        elif self._units == 'points':
            sc = fig.dpi / 72.0
        elif self._units == 'width':
            sc = ax.bbox.width
        elif self._units == 'height':
            sc = ax.bbox.height
        elif self._units == 'dots':
            sc = 1.0
        else:
            raise ValueError('unrecognized units: %s' % self._units)

        self._transforms = np.zeros((len(self._widths), 3, 3))
        widths = self._widths * sc
        heights = self._heights * sc
        sin_angle = np.sin(self._angles)
        cos_angle = np.cos(self._angles)
        self._transforms[:, 0, 0] = widths * cos_angle
        self._transforms[:, 0, 1] = heights * -sin_angle
        self._transforms[:, 1, 0] = widths * sin_angle
        self._transforms[:, 1, 1] = heights * cos_angle
        self._transforms[:, 2, 2] = 1.0

        _affine = transforms.Affine2D
        if self._units == 'xy':
            m = ax.transData.get_affine().get_matrix().copy()
            m[:2, 2:] = 0
            self.set_transform(_affine(m))

    @artist.allow_rasterization
    def draw(self, renderer):
        self._set_transforms()
        Collection.draw(self, renderer)


class PatchCollection(Collection):
    """
    A generic collection of patches.

    This makes it easier to assign a color map to a heterogeneous
    collection of patches.

    This also may improve plotting speed, since PatchCollection will
    draw faster than a large number of patches.
    """

    def __init__(self, patches, match_original=False, **kwargs):
        """
        *patches*
            a sequence of Patch objects.  This list may include
            a heterogeneous assortment of different patch types.

        *match_original*
            If True, use the colors and linewidths of the original
            patches.  If False, new colors may be assigned by
            providing the standard collection arguments, facecolor,
            edgecolor, linewidths, norm or cmap.

        If any of *edgecolors*, *facecolors*, *linewidths*,
        *antialiaseds* are None, they default to their
        :data:`matplotlib.rcParams` patch setting, in sequence form.

        The use of :class:`~matplotlib.cm.ScalarMappable` is optional.
        If the :class:`~matplotlib.cm.ScalarMappable` matrix _A is not
        None (i.e., a call to set_array has been made), at draw time a
        call to scalar mappable will be made to set the face colors.
        """

        if match_original:
            def determine_facecolor(patch):
                if patch.get_fill():
                    return patch.get_facecolor()
                return [0, 0, 0, 0]

            kwargs['facecolors'] = [determine_facecolor(p) for p in patches]
            kwargs['edgecolors'] = [p.get_edgecolor() for p in patches]
            kwargs['linewidths'] = [p.get_linewidth() for p in patches]
            kwargs['linestyles'] = [p.get_linestyle() for p in patches]
            kwargs['antialiaseds'] = [p.get_antialiased() for p in patches]

        Collection.__init__(self, **kwargs)

        self.set_paths(patches)

    def set_paths(self, patches):
        paths = [p.get_transform().transform_path(p.get_path())
                 for p in patches]
        self._paths = paths


class TriMesh(Collection):
    """
    Class for the efficient drawing of a triangular mesh using
    Gouraud shading.

    A triangular mesh is a :class:`~matplotlib.tri.Triangulation`
    object.
    """
    def __init__(self, triangulation, **kwargs):
        Collection.__init__(self, **kwargs)
        self._triangulation = triangulation
        self._shading = 'gouraud'
        self._is_filled = True

        self._bbox = transforms.Bbox.unit()

        # Unfortunately this requires a copy, unless Triangulation
        # was rewritten.
        xy = np.hstack((triangulation.x.reshape(-1, 1),
                        triangulation.y.reshape(-1, 1)))
        self._bbox.update_from_data_xy(xy)

    def get_paths(self):
        if self._paths is None:
            self.set_paths()
        return self._paths

    def set_paths(self):
        self._paths = self.convert_mesh_to_paths(self._triangulation)

    @staticmethod
    def convert_mesh_to_paths(tri):
        """
        Converts a given mesh into a sequence of
        :class:`matplotlib.path.Path` objects for easier rendering by
        backends that do not directly support meshes.

        This function is primarily of use to backend implementers.
        """
        Path = mpath.Path
        triangles = tri.get_masked_triangles()
        verts = np.concatenate((tri.x[triangles][..., np.newaxis],
                                tri.y[triangles][..., np.newaxis]), axis=2)
        return [Path(x) for x in verts]

    @artist.allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        renderer.open_group(self.__class__.__name__)
        transform = self.get_transform()

        # Get a list of triangles and the color at each vertex.
        tri = self._triangulation
        triangles = tri.get_masked_triangles()

        verts = np.concatenate((tri.x[triangles][..., np.newaxis],
                                tri.y[triangles][..., np.newaxis]), axis=2)

        self.update_scalarmappable()
        colors = self._facecolors[triangles]

        gc = renderer.new_gc()
        self._set_gc_clip(gc)
        gc.set_linewidth(self.get_linewidth()[0])
        renderer.draw_gouraud_triangles(gc, verts, colors, transform.frozen())
        gc.restore()
        renderer.close_group(self.__class__.__name__)


class QuadMesh(Collection):
    """
    Class for the efficient drawing of a quadrilateral mesh.

    A quadrilateral mesh consists of a grid of vertices. The
    dimensions of this array are (*meshWidth* + 1, *meshHeight* +
    1). Each vertex in the mesh has a different set of "mesh
    coordinates" representing its position in the topology of the
    mesh. For any values (*m*, *n*) such that 0 <= *m* <= *meshWidth*
    and 0 <= *n* <= *meshHeight*, the vertices at mesh coordinates
    (*m*, *n*), (*m*, *n* + 1), (*m* + 1, *n* + 1), and (*m* + 1, *n*)
    form one of the quadrilaterals in the mesh. There are thus
    (*meshWidth* * *meshHeight*) quadrilaterals in the mesh.  The mesh
    need not be regular and the polygons need not be convex.

    A quadrilateral mesh is represented by a (2 x ((*meshWidth* + 1) *
    (*meshHeight* + 1))) numpy array *coordinates*, where each row is
    the *x* and *y* coordinates of one of the vertices.  To define the
    function that maps from a data point to its corresponding color,
    use the :meth:`set_cmap` method.  Each of these arrays is indexed in
    row-major order by the mesh coordinates of the vertex (or the mesh
    coordinates of the lower left vertex, in the case of the
    colors).

    For example, the first entry in *coordinates* is the
    coordinates of the vertex at mesh coordinates (0, 0), then the one
    at (0, 1), then at (0, 2) .. (0, meshWidth), (1, 0), (1, 1), and
    so on.

    *shading* may be 'flat', or 'gouraud'
    """
    def __init__(self, meshWidth, meshHeight, coordinates,
                 antialiased=True, shading='flat', **kwargs):
        Collection.__init__(self, **kwargs)
        self._meshWidth = meshWidth
        self._meshHeight = meshHeight
        # By converting to floats now, we can avoid that on every draw.
        self._coordinates = np.asarray(coordinates, float).reshape(
            (meshHeight + 1, meshWidth + 1, 2))
        self._antialiased = antialiased
        self._shading = shading

        self._bbox = transforms.Bbox.unit()
        self._bbox.update_from_data_xy(coordinates.reshape(
            ((meshWidth + 1) * (meshHeight + 1), 2)))

    def get_paths(self):
        if self._paths is None:
            self.set_paths()
        return self._paths

    def set_paths(self):
        self._paths = self.convert_mesh_to_paths(
            self._meshWidth, self._meshHeight, self._coordinates)
        self.stale = True

    def get_datalim(self, transData):
        return (self.get_transform() - transData).transform_bbox(self._bbox)

    @staticmethod
    def convert_mesh_to_paths(meshWidth, meshHeight, coordinates):
        """
        Converts a given mesh into a sequence of
        :class:`matplotlib.path.Path` objects for easier rendering by
        backends that do not directly support quadmeshes.

        This function is primarily of use to backend implementers.
        """
        Path = mpath.Path

        if isinstance(coordinates, np.ma.MaskedArray):
            c = coordinates.data
        else:
            c = coordinates

        points = np.concatenate((
                    c[0:-1, 0:-1],
                    c[0:-1, 1:],
                    c[1:, 1:],
                    c[1:, 0:-1],
                    c[0:-1, 0:-1]
                ), axis=2)
        points = points.reshape((meshWidth * meshHeight, 5, 2))
        return [Path(x) for x in points]

    def convert_mesh_to_triangles(self, meshWidth, meshHeight, coordinates):
        """
        Converts a given mesh into a sequence of triangles, each point
        with its own color.  This is useful for experiments using
        `draw_qouraud_triangle`.
        """
        if isinstance(coordinates, np.ma.MaskedArray):
            p = coordinates.data
        else:
            p = coordinates

        p_a = p[:-1, :-1]
        p_b = p[:-1, 1:]
        p_c = p[1:, 1:]
        p_d = p[1:, :-1]
        p_center = (p_a + p_b + p_c + p_d) / 4.0

        triangles = np.concatenate((
                p_a, p_b, p_center,
                p_b, p_c, p_center,
                p_c, p_d, p_center,
                p_d, p_a, p_center,
            ), axis=2)
        triangles = triangles.reshape((meshWidth * meshHeight * 4, 3, 2))

        c = self.get_facecolor().reshape((meshHeight + 1, meshWidth + 1, 4))
        c_a = c[:-1, :-1]
        c_b = c[:-1, 1:]
        c_c = c[1:, 1:]
        c_d = c[1:, :-1]
        c_center = (c_a + c_b + c_c + c_d) / 4.0

        colors = np.concatenate((
                        c_a, c_b, c_center,
                        c_b, c_c, c_center,
                        c_c, c_d, c_center,
                        c_d, c_a, c_center,
                    ), axis=2)
        colors = colors.reshape((meshWidth * meshHeight * 4, 3, 4))

        return triangles, colors

    @artist.allow_rasterization
    def draw(self, renderer):
        if not self.get_visible():
            return
        renderer.open_group(self.__class__.__name__, self.get_gid())
        transform = self.get_transform()
        transOffset = self.get_offset_transform()
        offsets = self._offsets

        if self.have_units():
            if len(self._offsets):
                xs = self.convert_xunits(self._offsets[:, 0])
                ys = self.convert_yunits(self._offsets[:, 1])
                offsets = np.column_stack([xs, ys])

        self.update_scalarmappable()

        if not transform.is_affine:
            coordinates = self._coordinates.reshape((-1, 2))
            coordinates = transform.transform(coordinates)
            coordinates = coordinates.reshape(self._coordinates.shape)
            transform = transforms.IdentityTransform()
        else:
            coordinates = self._coordinates

        if not transOffset.is_affine:
            offsets = transOffset.transform_non_affine(offsets)
            transOffset = transOffset.get_affine()

        gc = renderer.new_gc()
        self._set_gc_clip(gc)
        gc.set_linewidth(self.get_linewidth()[0])

        if self._shading == 'gouraud':
            triangles, colors = self.convert_mesh_to_triangles(
                self._meshWidth, self._meshHeight, coordinates)
            renderer.draw_gouraud_triangles(
                gc, triangles, colors, transform.frozen())
        else:
            renderer.draw_quad_mesh(
                gc, transform.frozen(), self._meshWidth, self._meshHeight,
                coordinates, offsets, transOffset, self.get_facecolor(),
                self._antialiased, self.get_edgecolors())
        gc.restore()
        renderer.close_group(self.__class__.__name__)
        self.stale = False


patchstr = artist.kwdoc(Collection)
for k in ('QuadMesh', 'TriMesh', 'PolyCollection', 'BrokenBarHCollection',
          'RegularPolyCollection', 'PathCollection',
          'StarPolygonCollection', 'PatchCollection',
          'CircleCollection', 'Collection',):
    docstring.interpd.update({k: patchstr})
docstring.interpd.update(LineCollection=artist.kwdoc(LineCollection))
