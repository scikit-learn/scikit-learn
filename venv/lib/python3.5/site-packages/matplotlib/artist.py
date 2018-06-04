from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from collections import OrderedDict, namedtuple
from functools import wraps
import inspect
import re
import warnings

import numpy as np

import matplotlib
from . import cbook, docstring, rcParams
from .path import Path
from .transforms import (Bbox, IdentityTransform, Transform, TransformedBbox,
                         TransformedPatchPath, TransformedPath)
# Note, matplotlib artists use the doc strings for set and get
# methods to enable the introspection methods of setp and getp.  Every
# set_* method should have a docstring containing the line
#
# ACCEPTS: [ legal | values ]
#
# and aliases for setters and getters should have a docstring that
# starts with 'alias for ', as in 'alias for set_somemethod'
#
# You may wonder why we use so much boiler-plate manually defining the
# set_alias and get_alias functions, rather than using some clever
# python trick.  The answer is that I need to be able to manipulate
# the docstring, and there is no clever way to do that in python 2.2,
# as far as I can see - see
#
# https://mail.python.org/pipermail/python-list/2004-October/242925.html


def allow_rasterization(draw):
    """
    Decorator for Artist.draw method. Provides routines
    that run before and after the draw call. The before and after functions
    are useful for changing artist-dependent renderer attributes or making
    other setup function calls, such as starting and flushing a mixed-mode
    renderer.
    """

    # the axes class has a second argument inframe for its draw method.
    @wraps(draw)
    def draw_wrapper(artist, renderer, *args, **kwargs):
        try:
            if artist.get_rasterized():
                renderer.start_rasterizing()
            if artist.get_agg_filter() is not None:
                renderer.start_filter()

            return draw(artist, renderer, *args, **kwargs)
        finally:
            if artist.get_agg_filter() is not None:
                renderer.stop_filter(artist.get_agg_filter())
            if artist.get_rasterized():
                renderer.stop_rasterizing()

    draw_wrapper._supports_rasterization = True
    return draw_wrapper


def _stale_axes_callback(self, val):
    if self.axes:
        self.axes.stale = val


_XYPair = namedtuple("_XYPair", "x y")


class Artist(object):
    """
    Abstract base class for someone who renders into a
    :class:`FigureCanvas`.
    """

    aname = 'Artist'
    zorder = 0
    # order of precedence when bulk setting/updating properties
    # via update.  The keys should be property names and the values
    # integers
    _prop_order = dict(color=-1)

    def __init__(self):
        self._stale = True
        self.stale_callback = None
        self._axes = None
        self.figure = None

        self._transform = None
        self._transformSet = False
        self._visible = True
        self._animated = False
        self._alpha = None
        self.clipbox = None
        self._clippath = None
        self._clipon = True
        self._label = ''
        self._picker = None
        self._contains = None
        self._rasterized = None
        self._agg_filter = None
        self._mouseover = False
        self.eventson = False  # fire events only if eventson
        self._oid = 0  # an observer id
        self._propobservers = {}  # a dict from oids to funcs
        try:
            self.axes = None
        except AttributeError:
            # Handle self.axes as a read-only property, as in Figure.
            pass
        self._remove_method = None
        self._url = None
        self._gid = None
        self._snap = None
        self._sketch = rcParams['path.sketch']
        self._path_effects = rcParams['path.effects']
        self._sticky_edges = _XYPair([], [])

    def __getstate__(self):
        d = self.__dict__.copy()
        # remove the unpicklable remove method, this will get re-added on load
        # (by the axes) if the artist lives on an axes.
        d['_remove_method'] = None
        d['stale_callback'] = None
        return d

    def remove(self):
        """
        Remove the artist from the figure if possible.  The effect
        will not be visible until the figure is redrawn, e.g., with
        :meth:`matplotlib.axes.Axes.draw_idle`.  Call
        :meth:`matplotlib.axes.Axes.relim` to update the axes limits
        if desired.

        Note: :meth:`~matplotlib.axes.Axes.relim` will not see
        collections even if the collection was added to axes with
        *autolim* = True.

        Note: there is no support for removing the artist's legend entry.
        """

        # There is no method to set the callback.  Instead the parent should
        # set the _remove_method attribute directly.  This would be a
        # protected attribute if Python supported that sort of thing.  The
        # callback has one parameter, which is the child to be removed.
        if self._remove_method is not None:
            self._remove_method(self)
            # clear stale callback
            self.stale_callback = None
            _ax_flag = False
            if hasattr(self, 'axes') and self.axes:
                # remove from the mouse hit list
                self.axes.mouseover_set.discard(self)
                # mark the axes as stale
                self.axes.stale = True
                # decouple the artist from the axes
                self.axes = None
                _ax_flag = True

            if self.figure:
                self.figure = None
                if not _ax_flag:
                    self.figure = True

        else:
            raise NotImplementedError('cannot remove artist')
        # TODO: the fix for the collections relim problem is to move the
        # limits calculation into the artist itself, including the property of
        # whether or not the artist should affect the limits.  Then there will
        # be no distinction between axes.add_line, axes.add_patch, etc.
        # TODO: add legend support

    def have_units(self):
        'Return *True* if units are set on the *x* or *y* axes'
        ax = self.axes
        if ax is None or ax.xaxis is None:
            return False
        return ax.xaxis.have_units() or ax.yaxis.have_units()

    def convert_xunits(self, x):
        """For artists in an axes, if the xaxis has units support,
        convert *x* using xaxis unit type
        """
        ax = getattr(self, 'axes', None)
        if ax is None or ax.xaxis is None:
            return x
        return ax.xaxis.convert_units(x)

    def convert_yunits(self, y):
        """For artists in an axes, if the yaxis has units support,
        convert *y* using yaxis unit type
        """
        ax = getattr(self, 'axes', None)
        if ax is None or ax.yaxis is None:
            return y
        return ax.yaxis.convert_units(y)

    @property
    def axes(self):
        """
        The :class:`~matplotlib.axes.Axes` instance the artist
        resides in, or *None*.
        """
        return self._axes

    @axes.setter
    def axes(self, new_axes):
        if (new_axes is not None and self._axes is not None
                and new_axes != self._axes):
            raise ValueError("Can not reset the axes.  You are probably "
                             "trying to re-use an artist in more than one "
                             "Axes which is not supported")
        self._axes = new_axes
        if new_axes is not None and new_axes is not self:
            self.stale_callback = _stale_axes_callback
        return new_axes

    @property
    def stale(self):
        """
        If the artist is 'stale' and needs to be re-drawn for the output to
        match the internal state of the artist.
        """
        return self._stale

    @stale.setter
    def stale(self, val):
        self._stale = val

        # if the artist is animated it does not take normal part in the
        # draw stack and is not expected to be drawn as part of the normal
        # draw loop (when not saving) so do not propagate this change
        if self.get_animated():
            return

        if val and self.stale_callback is not None:
            self.stale_callback(self, val)

    def get_window_extent(self, renderer):
        """
        Get the axes bounding box in display space.
        Subclasses should override for inclusion in the bounding box
        "tight" calculation. Default is to return an empty bounding
        box at 0, 0.

        Be careful when using this function, the results will not update
        if the artist window extent of the artist changes.  The extent
        can change due to any changes in the transform stack, such as
        changing the axes limits, the figure size, or the canvas used
        (as is done when saving a figure).  This can lead to unexpected
        behavior where interactive figures will look fine on the screen,
        but will save incorrectly.
        """
        return Bbox([[0, 0], [0, 0]])

    def add_callback(self, func):
        """
        Adds a callback function that will be called whenever one of
        the :class:`Artist`'s properties changes.

        Returns an *id* that is useful for removing the callback with
        :meth:`remove_callback` later.
        """
        oid = self._oid
        self._propobservers[oid] = func
        self._oid += 1
        return oid

    def remove_callback(self, oid):
        """
        Remove a callback based on its *id*.

        .. seealso::

            :meth:`add_callback`
               For adding callbacks

        """
        try:
            del self._propobservers[oid]
        except KeyError:
            pass

    def pchanged(self):
        """
        Fire an event when property changed, calling all of the
        registered callbacks.
        """
        for oid, func in six.iteritems(self._propobservers):
            func(self)

    def is_transform_set(self):
        """
        Returns *True* if :class:`Artist` has a transform explicitly
        set.
        """
        return self._transformSet

    def set_transform(self, t):
        """
        Set the artist transform.

        Parameters
        ----------
        t : `.Transform`
            .. ACCEPTS: `.Transform`
        """
        self._transform = t
        self._transformSet = True
        self.pchanged()
        self.stale = True

    def get_transform(self):
        """
        Return the :class:`~matplotlib.transforms.Transform`
        instance used by this artist.
        """
        if self._transform is None:
            self._transform = IdentityTransform()
        elif (not isinstance(self._transform, Transform)
              and hasattr(self._transform, '_as_mpl_transform')):
            self._transform = self._transform._as_mpl_transform(self.axes)
        return self._transform

    @cbook.deprecated("2.2")
    def hitlist(self, event):
        """
        List the children of the artist which contain the mouse event *event*.
        """
        L = []
        try:
            hascursor, info = self.contains(event)
            if hascursor:
                L.append(self)
        except:
            import traceback
            traceback.print_exc()
            print("while checking", self.__class__)

        for a in self.get_children():
            L.extend(a.hitlist(event))
        return L

    def get_children(self):
        """
        Return a list of the child :class:`Artist`s this
        :class:`Artist` contains.
        """
        return []

    def contains(self, mouseevent):
        """Test whether the artist contains the mouse event.

        Returns the truth value and a dictionary of artist specific details of
        selection, such as which points are contained in the pick radius.  See
        individual artists for details.
        """
        if callable(self._contains):
            return self._contains(self, mouseevent)
        warnings.warn("'%s' needs 'contains' method" % self.__class__.__name__)
        return False, {}

    def set_contains(self, picker):
        """
        Replace the contains test used by this artist. The new picker
        should be a callable function which determines whether the
        artist is hit by the mouse event::

            hit, props = picker(artist, mouseevent)

        If the mouse event is over the artist, return *hit* = *True*
        and *props* is a dictionary of properties you want returned
        with the contains test.

        Parameters
        ----------
        picker : callable
            .. ACCEPTS: a callable function
        """
        self._contains = picker

    def get_contains(self):
        """
        Return the _contains test used by the artist, or *None* for default.
        """
        return self._contains

    def pickable(self):
        'Return *True* if :class:`Artist` is pickable.'
        return (self.figure is not None and
                self.figure.canvas is not None and
                self._picker is not None)

    def pick(self, mouseevent):
        """
        Process pick event

        each child artist will fire a pick event if *mouseevent* is over
        the artist and the artist has picker set
        """
        # Pick self
        if self.pickable():
            picker = self.get_picker()
            if callable(picker):
                inside, prop = picker(self, mouseevent)
            else:
                inside, prop = self.contains(mouseevent)
            if inside:
                self.figure.canvas.pick_event(mouseevent, self, **prop)

        # Pick children
        for a in self.get_children():
            # make sure the event happened in the same axes
            ax = getattr(a, 'axes', None)
            if (mouseevent.inaxes is None or ax is None
                    or mouseevent.inaxes == ax):
                # we need to check if mouseevent.inaxes is None
                # because some objects associated with an axes (e.g., a
                # tick label) can be outside the bounding box of the
                # axes and inaxes will be None
                # also check that ax is None so that it traverse objects
                # which do no have an axes property but children might
                a.pick(mouseevent)

    def set_picker(self, picker):
        """
        Set the epsilon for picking used by this artist

        *picker* can be one of the following:

          * *None*: picking is disabled for this artist (default)

          * A boolean: if *True* then picking will be enabled and the
            artist will fire a pick event if the mouse event is over
            the artist

          * A float: if picker is a number it is interpreted as an
            epsilon tolerance in points and the artist will fire
            off an event if it's data is within epsilon of the mouse
            event.  For some artists like lines and patch collections,
            the artist may provide additional data to the pick event
            that is generated, e.g., the indices of the data within
            epsilon of the pick event

          * A function: if picker is callable, it is a user supplied
            function which determines whether the artist is hit by the
            mouse event::

              hit, props = picker(artist, mouseevent)

            to determine the hit test.  if the mouse event is over the
            artist, return *hit=True* and props is a dictionary of
            properties you want added to the PickEvent attributes.

        Parameters
        ----------
        picker : None or bool or float or callable
            .. ACCEPTS: [None | bool | float | callable]
        """
        self._picker = picker

    def get_picker(self):
        """Return the picker object used by this artist."""
        return self._picker

    @cbook.deprecated("2.2", "artist.figure is not None")
    def is_figure_set(self):
        """Returns whether the artist is assigned to a `.Figure`."""
        return self.figure is not None

    def get_url(self):
        """Returns the url."""
        return self._url

    def set_url(self, url):
        """
        Sets the url for the artist.

        Parameters
        ----------
        url : str
            .. ACCEPTS: a url string
        """
        self._url = url

    def get_gid(self):
        """Returns the group id."""
        return self._gid

    def set_gid(self, gid):
        """
        Sets the (group) id for the artist.

        Parameters
        ----------
        gid : str
            .. ACCEPTS: an id string
        """
        self._gid = gid

    def get_snap(self):
        """
        Returns the snap setting which may be:

          * True: snap vertices to the nearest pixel center

          * False: leave vertices as-is

          * None: (auto) If the path contains only rectilinear line
            segments, round to the nearest pixel center

        Only supported by the Agg and MacOSX backends.
        """
        if rcParams['path.snap']:
            return self._snap
        else:
            return False

    def set_snap(self, snap):
        """
        Sets the snap setting which may be:

          * True: snap vertices to the nearest pixel center

          * False: leave vertices as-is

          * None: (auto) If the path contains only rectilinear line
            segments, round to the nearest pixel center

        Only supported by the Agg and MacOSX backends.

        Parameters
        ----------
        snap : bool or None
            .. ACCEPTS: bool or None
        """
        self._snap = snap
        self.stale = True

    def get_sketch_params(self):
        """
        Returns the sketch parameters for the artist.

        Returns
        -------
        sketch_params : tuple or `None`

        A 3-tuple with the following elements:

          * `scale`: The amplitude of the wiggle perpendicular to the
            source line.

          * `length`: The length of the wiggle along the line.

          * `randomness`: The scale factor by which the length is
            shrunken or expanded.

        May return `None` if no sketch parameters were set.
        """
        return self._sketch

    def set_sketch_params(self, scale=None, length=None, randomness=None):
        """
        Sets the sketch parameters.

        Parameters
        ----------

        scale : float, optional
            The amplitude of the wiggle perpendicular to the source
            line, in pixels.  If scale is `None`, or not provided, no
            sketch filter will be provided.

        length : float, optional
             The length of the wiggle along the line, in pixels
             (default 128.0)

        randomness : float, optional
            The scale factor by which the length is shrunken or
            expanded (default 16.0)

            .. ACCEPTS: (scale: float, length: float, randomness: float)
        """
        if scale is None:
            self._sketch = None
        else:
            self._sketch = (scale, length or 128.0, randomness or 16.0)
        self.stale = True

    def set_path_effects(self, path_effects):
        """Set the path effects.

        Parameters
        ----------
        path_effects : `.AbstractPathEffect`
            .. ACCEPTS: `.AbstractPathEffect`
        """
        self._path_effects = path_effects
        self.stale = True

    def get_path_effects(self):
        return self._path_effects

    def get_figure(self):
        """Return the `.Figure` instance the artist belongs to."""
        return self.figure

    def set_figure(self, fig):
        """
        Set the `.Figure` instance the artist belongs to.

        Parameters
        ----------
        fig : `.Figure`
            .. ACCEPTS: a `.Figure` instance
        """
        # if this is a no-op just return
        if self.figure is fig:
            return
        # if we currently have a figure (the case of both `self.figure`
        # and `fig` being none is taken care of above) we then user is
        # trying to change the figure an artist is associated with which
        # is not allowed for the same reason as adding the same instance
        # to more than one Axes
        if self.figure is not None:
            raise RuntimeError("Can not put single artist in "
                               "more than one figure")
        self.figure = fig
        if self.figure and self.figure is not self:
            self.pchanged()
        self.stale = True

    def set_clip_box(self, clipbox):
        """
        Set the artist's clip `.Bbox`.

        Parameters
        ----------
        clipbox : `.Bbox`
            .. ACCEPTS: a `.Bbox` instance
        """
        self.clipbox = clipbox
        self.pchanged()
        self.stale = True

    def set_clip_path(self, path, transform=None):
        """
        Set the artist's clip path, which may be:

        - a :class:`~matplotlib.patches.Patch` (or subclass) instance; or
        - a :class:`~matplotlib.path.Path` instance, in which case a
          :class:`~matplotlib.transforms.Transform` instance, which will be
          applied to the path before using it for clipping, must be provided;
          or
        - ``None``, to remove a previously set clipping path.

        For efficiency, if the path happens to be an axis-aligned rectangle,
        this method will set the clipping box to the corresponding rectangle
        and set the clipping path to ``None``.

        ACCEPTS: [(`~matplotlib.path.Path`, `.Transform`) | `.Patch` | None]
        """
        from matplotlib.patches import Patch, Rectangle

        success = False
        if transform is None:
            if isinstance(path, Rectangle):
                self.clipbox = TransformedBbox(Bbox.unit(),
                                               path.get_transform())
                self._clippath = None
                success = True
            elif isinstance(path, Patch):
                self._clippath = TransformedPatchPath(path)
                success = True
            elif isinstance(path, tuple):
                path, transform = path

        if path is None:
            self._clippath = None
            success = True
        elif isinstance(path, Path):
            self._clippath = TransformedPath(path, transform)
            success = True
        elif isinstance(path, TransformedPatchPath):
            self._clippath = path
            success = True
        elif isinstance(path, TransformedPath):
            self._clippath = path
            success = True

        if not success:
            raise TypeError(
                "Invalid arguments to set_clip_path, of type {} and {}"
                .format(type(path).__name__, type(transform).__name__))
        # This may result in the callbacks being hit twice, but guarantees they
        # will be hit at least once.
        self.pchanged()
        self.stale = True

    def get_alpha(self):
        """
        Return the alpha value used for blending - not supported on all
        backends
        """
        return self._alpha

    def get_visible(self):
        "Return the artist's visiblity"
        return self._visible

    def get_animated(self):
        "Return the artist's animated state"
        return self._animated

    def get_clip_on(self):
        'Return whether artist uses clipping'
        return self._clipon

    def get_clip_box(self):
        'Return artist clipbox'
        return self.clipbox

    def get_clip_path(self):
        'Return artist clip path'
        return self._clippath

    def get_transformed_clip_path_and_affine(self):
        '''
        Return the clip path with the non-affine part of its
        transformation applied, and the remaining affine part of its
        transformation.
        '''
        if self._clippath is not None:
            return self._clippath.get_transformed_path_and_affine()
        return None, None

    def set_clip_on(self, b):
        """
        Set whether artist uses clipping.

        When False artists will be visible out side of the axes which
        can lead to unexpected results.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._clipon = b
        # This may result in the callbacks being hit twice, but ensures they
        # are hit at least once
        self.pchanged()
        self.stale = True

    def _set_gc_clip(self, gc):
        'Set the clip properly for the gc'
        if self._clipon:
            if self.clipbox is not None:
                gc.set_clip_rectangle(self.clipbox)
            gc.set_clip_path(self._clippath)
        else:
            gc.set_clip_rectangle(None)
            gc.set_clip_path(None)

    def get_rasterized(self):
        """Return whether the artist is to be rasterized."""
        return self._rasterized

    def set_rasterized(self, rasterized):
        """
        Force rasterized (bitmap) drawing in vector backend output.

        Defaults to None, which implies the backend's default behavior.

        Parameters
        ----------
        rasterized : bool or None
            .. ACCEPTS: bool or None
        """
        if rasterized and not hasattr(self.draw, "_supports_rasterization"):
            warnings.warn("Rasterization of '%s' will be ignored" % self)

        self._rasterized = rasterized

    def get_agg_filter(self):
        """Return filter function to be used for agg filter."""
        return self._agg_filter

    def set_agg_filter(self, filter_func):
        """Set the agg filter.

        Parameters
        ----------
        filter_func : callable
            A filter function, which takes a (m, n, 3) float array and a dpi
            value, and returns a (m, n, 3) array.

            .. ACCEPTS: a filter function, which takes a (m, n, 3) float array
                and a dpi value, and returns a (m, n, 3) array
        """
        self._agg_filter = filter_func
        self.stale = True

    def draw(self, renderer, *args, **kwargs):
        'Derived classes drawing method'
        if not self.get_visible():
            return
        self.stale = False

    def set_alpha(self, alpha):
        """
        Set the alpha value used for blending - not supported on
        all backends.

        Parameters
        ----------
        alpha : float
            .. ACCEPTS: float (0.0 transparent through 1.0 opaque)
        """
        self._alpha = alpha
        self.pchanged()
        self.stale = True

    def set_visible(self, b):
        """
        Set the artist's visibility.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        self._visible = b
        self.pchanged()
        self.stale = True

    def set_animated(self, b):
        """
        Set the artist's animation state.

        Parameters
        ----------
        b : bool
            .. ACCEPTS: bool
        """
        if self._animated != b:
            self._animated = b
            self.pchanged()

    def update(self, props):
        """
        Update this artist's properties from the dictionary *prop*.
        """
        def _update_property(self, k, v):
            """Sorting out how to update property (setter or setattr).

            Parameters
            ----------
            k : str
                The name of property to update
            v : obj
                The value to assign to the property

            Returns
            -------
            ret : obj or None
                If using a `set_*` method return it's return, else None.
            """
            k = k.lower()
            # white list attributes we want to be able to update through
            # art.update, art.set, setp
            if k in {'axes'}:
                return setattr(self, k, v)
            else:
                func = getattr(self, 'set_' + k, None)
                if not callable(func):
                    raise AttributeError('Unknown property %s' % k)
                return func(v)

        store = self.eventson
        self.eventson = False
        try:
            ret = [_update_property(self, k, v)
                   for k, v in props.items()]
        finally:
            self.eventson = store

        if len(ret):
            self.pchanged()
            self.stale = True
        return ret

    def get_label(self):
        """Get the label used for this artist in the legend."""
        return self._label

    def set_label(self, s):
        """
        Set the label to *s* for auto legend.

        Parameters
        ----------
        s : object
            *s* will be converted to a string by calling `str` (`unicode` on
            Py2).

            .. ACCEPTS: object
        """
        if s is not None:
            self._label = six.text_type(s)
        else:
            self._label = None
        self.pchanged()
        self.stale = True

    def get_zorder(self):
        """Return the artist's zorder."""
        return self.zorder

    def set_zorder(self, level):
        """
        Set the zorder for the artist.  Artists with lower zorder
        values are drawn first.

        Parameters
        ----------
        level : float
            .. ACCEPTS: float
        """
        if level is None:
            level = self.__class__.zorder
        self.zorder = level
        self.pchanged()
        self.stale = True

    @property
    def sticky_edges(self):
        """
        `x` and `y` sticky edge lists.

        When performing autoscaling, if a data limit coincides with a value in
        the corresponding sticky_edges list, then no margin will be added--the
        view limit "sticks" to the edge. A typical usecase is histograms,
        where one usually expects no margin on the bottom edge (0) of the
        histogram.

        This attribute cannot be assigned to; however, the `x` and `y` lists
        can be modified in place as needed.

        Examples
        --------

        >>> artist.sticky_edges.x[:] = (xmin, xmax)
        >>> artist.sticky_edges.y[:] = (ymin, ymax)

        """
        return self._sticky_edges

    def update_from(self, other):
        'Copy properties from *other* to *self*.'
        self._transform = other._transform
        self._transformSet = other._transformSet
        self._visible = other._visible
        self._alpha = other._alpha
        self.clipbox = other.clipbox
        self._clipon = other._clipon
        self._clippath = other._clippath
        self._label = other._label
        self._sketch = other._sketch
        self._path_effects = other._path_effects
        self.sticky_edges.x[:] = other.sticky_edges.x[:]
        self.sticky_edges.y[:] = other.sticky_edges.y[:]
        self.pchanged()
        self.stale = True

    def properties(self):
        """
        return a dictionary mapping property name -> value for all Artist props
        """
        return ArtistInspector(self).properties()

    def set(self, **kwargs):
        """A property batch setter. Pass *kwargs* to set properties.
        """
        props = OrderedDict(
            sorted(kwargs.items(), reverse=True,
                   key=lambda x: (self._prop_order.get(x[0], 0), x[0])))

        return self.update(props)

    def findobj(self, match=None, include_self=True):
        """
        Find artist objects.

        Recursively find all :class:`~matplotlib.artist.Artist` instances
        contained in self.

        *match* can be

          - None: return all objects contained in artist.

          - function with signature ``boolean = match(artist)``
            used to filter matches

          - class instance: e.g., Line2D.  Only return artists of class type.

        If *include_self* is True (default), include self in the list to be
        checked for a match.

        """
        if match is None:  # always return True
            def matchfunc(x):
                return True
        elif isinstance(match, type) and issubclass(match, Artist):
            def matchfunc(x):
                return isinstance(x, match)
        elif callable(match):
            matchfunc = match
        else:
            raise ValueError('match must be None, a matplotlib.artist.Artist '
                             'subclass, or a callable')

        artists = sum([c.findobj(matchfunc) for c in self.get_children()], [])
        if include_self and matchfunc(self):
            artists.append(self)
        return artists

    def get_cursor_data(self, event):
        """
        Get the cursor data for a given event.
        """
        return None

    def format_cursor_data(self, data):
        """
        Return *cursor data* string formatted.
        """
        try:
            data[0]
        except (TypeError, IndexError):
            data = [data]
        return ', '.join('{:0.3g}'.format(item) for item in data if
                isinstance(item, (np.floating, np.integer, int, float)))

    @property
    def mouseover(self):
        return self._mouseover

    @mouseover.setter
    def mouseover(self, val):
        val = bool(val)
        self._mouseover = val
        ax = self.axes
        if ax:
            if val:
                ax.mouseover_set.add(self)
            else:
                ax.mouseover_set.discard(self)


class ArtistInspector(object):
    """
    A helper class to inspect an :class:`~matplotlib.artist.Artist`
    and return information about it's settable properties and their
    current values.
    """
    def __init__(self, o):
        """
        Initialize the artist inspector with an
        :class:`~matplotlib.artist.Artist` or iterable of :class:`Artists`.
        If an iterable is used, we assume it is a homogeneous sequence (all
        :class:`Artists` are of the same type) and it is your responsibility
        to make sure this is so.
        """
        if not isinstance(o, Artist):
            if cbook.iterable(o):
                o = list(o)
                if len(o):
                    o = o[0]

        self.oorig = o
        if not inspect.isclass(o):
            o = type(o)
        self.o = o

        self.aliasd = self.get_aliases()

    def get_aliases(self):
        """
        Get a dict mapping *fullname* -> *alias* for each *alias* in
        the :class:`~matplotlib.artist.ArtistInspector`.

        e.g., for lines::

          {'markerfacecolor': 'mfc',
           'linewidth'      : 'lw',
          }

        """
        names = [name for name in dir(self.o)
                 if name.startswith(('set_', 'get_'))
                    and callable(getattr(self.o, name))]
        aliases = {}
        for name in names:
            func = getattr(self.o, name)
            if not self.is_alias(func):
                continue
            docstring = func.__doc__
            fullname = docstring[10:]
            aliases.setdefault(fullname[4:], {})[name[4:]] = None
        return aliases

    _get_valid_values_regex = re.compile(
        r"\n\s*(?:\.\.\s+)?ACCEPTS:\s*((?:.|\n)*?)(?:$|(?:\n\n))"
    )

    def get_valid_values(self, attr):
        """
        Get the legal arguments for the setter associated with *attr*.

        This is done by querying the docstring of the function *set_attr*
        for a line that begins with "ACCEPTS" or ".. ACCEPTS":

        e.g., for a line linestyle, return
        "[ ``'-'`` | ``'--'`` | ``'-.'`` | ``':'`` | ``'steps'`` | ``'None'``
        ]"
        """

        name = 'set_%s' % attr
        if not hasattr(self.o, name):
            raise AttributeError('%s has no function %s' % (self.o, name))
        func = getattr(self.o, name)

        docstring = func.__doc__
        if docstring is None:
            return 'unknown'

        if docstring.startswith('alias for '):
            return None

        match = self._get_valid_values_regex.search(docstring)
        if match is not None:
            return re.sub("\n *", " ", match.group(1))
        return 'unknown'

    def _get_setters_and_targets(self):
        """
        Get the attribute strings and a full path to where the setter
        is defined for all setters in an object.
        """

        setters = []
        for name in dir(self.o):
            if not name.startswith('set_'):
                continue
            func = getattr(self.o, name)
            if not callable(func):
                continue
            if six.PY2:
                nargs = len(inspect.getargspec(func)[0])
            else:
                nargs = len(inspect.getfullargspec(func)[0])
            if nargs < 2 or self.is_alias(func):
                continue
            source_class = self.o.__module__ + "." + self.o.__name__
            for cls in self.o.mro():
                if name in cls.__dict__:
                    source_class = cls.__module__ + "." + cls.__name__
                    break
            setters.append((name[4:], source_class + "." + name))
        return setters

    def get_setters(self):
        """
        Get the attribute strings with setters for object.  e.g., for a line,
        return ``['markerfacecolor', 'linewidth', ....]``.
        """

        return [prop for prop, target in self._get_setters_and_targets()]

    def is_alias(self, o):
        """
        Return *True* if method object *o* is an alias for another
        function.
        """
        ds = o.__doc__
        if ds is None:
            return False
        return ds.startswith('alias for ')

    def aliased_name(self, s):
        """
        return 'PROPNAME or alias' if *s* has an alias, else return
        PROPNAME.

        e.g., for the line markerfacecolor property, which has an
        alias, return 'markerfacecolor or mfc' and for the transform
        property, which does not, return 'transform'
        """

        if s in self.aliasd:
            return s + ''.join([' or %s' % x
                                for x in sorted(self.aliasd[s])])
        else:
            return s

    def aliased_name_rest(self, s, target):
        """
        return 'PROPNAME or alias' if *s* has an alias, else return
        PROPNAME formatted for ReST

        e.g., for the line markerfacecolor property, which has an
        alias, return 'markerfacecolor or mfc' and for the transform
        property, which does not, return 'transform'
        """

        if s in self.aliasd:
            aliases = ''.join([' or %s' % x
                               for x in sorted(self.aliasd[s])])
        else:
            aliases = ''
        return ':meth:`%s <%s>`%s' % (s, target, aliases)

    def pprint_setters(self, prop=None, leadingspace=2):
        """
        If *prop* is *None*, return a list of strings of all settable
        properties and their valid values.

        If *prop* is not *None*, it is a valid property name and that
        property will be returned as a string of property : valid
        values.
        """
        if leadingspace:
            pad = ' ' * leadingspace
        else:
            pad = ''
        if prop is not None:
            accepts = self.get_valid_values(prop)
            return '%s%s: %s' % (pad, prop, accepts)

        attrs = self._get_setters_and_targets()
        attrs.sort()
        lines = []

        for prop, path in attrs:
            accepts = self.get_valid_values(prop)
            name = self.aliased_name(prop)

            lines.append('%s%s: %s' % (pad, name, accepts))
        return lines

    def pprint_setters_rest(self, prop=None, leadingspace=4):
        """
        If *prop* is *None*, return a list of strings of all settable
        properties and their valid values.  Format the output for ReST

        If *prop* is not *None*, it is a valid property name and that
        property will be returned as a string of property : valid
        values.
        """
        if leadingspace:
            pad = ' ' * leadingspace
        else:
            pad = ''
        if prop is not None:
            accepts = self.get_valid_values(prop)
            return '%s%s: %s' % (pad, prop, accepts)

        attrs = self._get_setters_and_targets()
        attrs.sort()
        lines = []

        ########
        names = [self.aliased_name_rest(prop, target)
                 for prop, target in attrs]
        accepts = [self.get_valid_values(prop) for prop, target in attrs]

        col0_len = max(len(n) for n in names)
        col1_len = max(len(a) for a in accepts)

        lines.append('')
        lines.append(pad + '.. table::')
        lines.append(pad + '   :class: property-table')
        pad += '   '

        table_formatstr = pad + '=' * col0_len + '   ' + '=' * col1_len

        lines.append('')
        lines.append(table_formatstr)
        lines.append(pad + 'Property'.ljust(col0_len + 3) +
                     'Description'.ljust(col1_len))
        lines.append(table_formatstr)

        lines.extend([pad + n.ljust(col0_len + 3) + a.ljust(col1_len)
                      for n, a in zip(names, accepts)])

        lines.append(table_formatstr)
        lines.append('')
        return lines

    def properties(self):
        """
        return a dictionary mapping property name -> value
        """
        o = self.oorig
        getters = [name for name in dir(o)
                   if name.startswith('get_') and callable(getattr(o, name))]
        getters.sort()
        d = dict()
        for name in getters:
            func = getattr(o, name)
            if self.is_alias(func):
                continue

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    val = func()
            except:
                continue
            else:
                d[name[4:]] = val

        return d

    def pprint_getters(self):
        """
        Return the getters and actual values as list of strings.
        """

        lines = []
        for name, val in sorted(six.iteritems(self.properties())):
            if getattr(val, 'shape', ()) != () and len(val) > 6:
                s = str(val[:6]) + '...'
            else:
                s = str(val)
            s = s.replace('\n', ' ')
            if len(s) > 50:
                s = s[:50] + '...'
            name = self.aliased_name(name)
            lines.append('    %s = %s' % (name, s))
        return lines


def getp(obj, property=None):
    """
    Return the value of object's property.  *property* is an optional string
    for the property you want to return

    Example usage::

        getp(obj)  # get all the object properties
        getp(obj, 'linestyle')  # get the linestyle property

    *obj* is a :class:`Artist` instance, e.g.,
    :class:`~matplotllib.lines.Line2D` or an instance of a
    :class:`~matplotlib.axes.Axes` or :class:`matplotlib.text.Text`.
    If the *property* is 'somename', this function returns

      obj.get_somename()

    :func:`getp` can be used to query all the gettable properties with
    ``getp(obj)``. Many properties have aliases for shorter typing, e.g.
    'lw' is an alias for 'linewidth'.  In the output, aliases and full
    property names will be listed as:

      property or alias = value

    e.g.:

      linewidth or lw = 2
    """
    if property is None:
        insp = ArtistInspector(obj)
        ret = insp.pprint_getters()
        print('\n'.join(ret))
        return

    func = getattr(obj, 'get_' + property)
    return func()

# alias
get = getp


def setp(obj, *args, **kwargs):
    """
    Set a property on an artist object.

    matplotlib supports the use of :func:`setp` ("set property") and
    :func:`getp` to set and get object properties, as well as to do
    introspection on the object.  For example, to set the linestyle of a
    line to be dashed, you can do::

      >>> line, = plot([1,2,3])
      >>> setp(line, linestyle='--')

    If you want to know the valid types of arguments, you can provide
    the name of the property you want to set without a value::

      >>> setp(line, 'linestyle')
          linestyle: [ '-' | '--' | '-.' | ':' | 'steps' | 'None' ]

    If you want to see all the properties that can be set, and their
    possible values, you can do::

      >>> setp(line)
          ... long output listing omitted

    You may specify another output file to `setp` if `sys.stdout` is not
    acceptable for some reason using the `file` keyword-only argument::

      >>> with fopen('output.log') as f:
      >>>     setp(line, file=f)

    :func:`setp` operates on a single instance or a iterable of
    instances. If you are in query mode introspecting the possible
    values, only the first instance in the sequence is used. When
    actually setting values, all the instances will be set.  e.g.,
    suppose you have a list of two lines, the following will make both
    lines thicker and red::

      >>> x = arange(0,1.0,0.01)
      >>> y1 = sin(2*pi*x)
      >>> y2 = sin(4*pi*x)
      >>> lines = plot(x, y1, x, y2)
      >>> setp(lines, linewidth=2, color='r')

    :func:`setp` works with the MATLAB style string/value pairs or
    with python kwargs.  For example, the following are equivalent::

      >>> setp(lines, 'linewidth', 2, 'color', 'r')  # MATLAB style
      >>> setp(lines, linewidth=2, color='r')        # python style
    """

    if isinstance(obj, Artist):
        objs = [obj]
    else:
        objs = list(cbook.flatten(obj))

    if not objs:
        return

    insp = ArtistInspector(objs[0])

    # file has to be popped before checking if kwargs is empty
    printArgs = {}
    if 'file' in kwargs:
        printArgs['file'] = kwargs.pop('file')

    if not kwargs and len(args) < 2:
        if args:
            print(insp.pprint_setters(prop=args[0]), **printArgs)
        else:
            print('\n'.join(insp.pprint_setters()), **printArgs)
        return

    if len(args) % 2:
        raise ValueError('The set args must be string, value pairs')

    # put args into ordereddict to maintain order
    funcvals = OrderedDict()
    for i in range(0, len(args) - 1, 2):
        funcvals[args[i]] = args[i + 1]

    ret = [o.update(funcvals) for o in objs]
    ret.extend([o.set(**kwargs) for o in objs])
    return [x for x in cbook.flatten(ret)]


def kwdoc(a):
    hardcopy = matplotlib.rcParams['docstring.hardcopy']
    if hardcopy:
        return '\n'.join(ArtistInspector(a).pprint_setters_rest(
                         leadingspace=4))
    else:
        return '\n'.join(ArtistInspector(a).pprint_setters(leadingspace=2))

docstring.interpd.update(Artist=kwdoc(Artist))
