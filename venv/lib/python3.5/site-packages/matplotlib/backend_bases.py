"""
Abstract base classes define the primitives that renderers and
graphics contexts must implement to serve as a matplotlib backend

:class:`RendererBase`
    An abstract base class to handle drawing/rendering operations.

:class:`FigureCanvasBase`
    The abstraction layer that separates the
    :class:`matplotlib.figure.Figure` from the backend specific
    details like a user interface drawing area

:class:`GraphicsContextBase`
    An abstract base class that provides color, line styles, etc...

:class:`Event`
    The base class for all of the matplotlib event
    handling.  Derived classes such as :class:`KeyEvent` and
    :class:`MouseEvent` store the meta data like keys and buttons
    pressed, x and y locations in pixel and
    :class:`~matplotlib.axes.Axes` coordinates.

:class:`ShowBase`
    The base class for the Show class of each interactive backend;
    the 'show' callable is then set to Show.__call__, inherited from
    ShowBase.

:class:`ToolContainerBase`
     The base class for the Toolbar class of each interactive backend.

:class:`StatusbarBase`
    The base class for the messaging area.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange

from contextlib import contextmanager
from functools import partial
import importlib
import io
import os
import sys
import time
import warnings
from weakref import WeakKeyDictionary

import numpy as np

from matplotlib import (
    backend_tools as tools, cbook, colors, textpath, tight_bbox, transforms,
    widgets, get_backend, is_interactive, rcParams)
from matplotlib._pylab_helpers import Gcf
from matplotlib.transforms import Bbox, TransformedBbox, Affine2D
from matplotlib.path import Path

try:
    from PIL import Image
    _has_pil = True
    del Image
except ImportError:
    _has_pil = False


_default_filetypes = {
    'ps': 'Postscript',
    'eps': 'Encapsulated Postscript',
    'pdf': 'Portable Document Format',
    'pgf': 'PGF code for LaTeX',
    'png': 'Portable Network Graphics',
    'raw': 'Raw RGBA bitmap',
    'rgba': 'Raw RGBA bitmap',
    'svg': 'Scalable Vector Graphics',
    'svgz': 'Scalable Vector Graphics'
}


_default_backends = {
    'ps': 'matplotlib.backends.backend_ps',
    'eps': 'matplotlib.backends.backend_ps',
    'pdf': 'matplotlib.backends.backend_pdf',
    'pgf': 'matplotlib.backends.backend_pgf',
    'png': 'matplotlib.backends.backend_agg',
    'raw': 'matplotlib.backends.backend_agg',
    'rgba': 'matplotlib.backends.backend_agg',
    'svg': 'matplotlib.backends.backend_svg',
    'svgz': 'matplotlib.backends.backend_svg',
}


def register_backend(format, backend, description=None):
    """
    Register a backend for saving to a given file format.

    Parameters
    ----------
    format : str
        File extension

    backend : module string or canvas class
        Backend for handling file output

    description : str, optional
        Description of the file type.  Defaults to an empty string
    """
    if description is None:
        description = ''
    _default_backends[format] = backend
    _default_filetypes[format] = description


def get_registered_canvas_class(format):
    """
    Return the registered default canvas for given file format.
    Handles deferred import of required backend.
    """
    if format not in _default_backends:
        return None
    backend_class = _default_backends[format]
    if isinstance(backend_class, six.string_types):
        backend_class = importlib.import_module(backend_class).FigureCanvas
        _default_backends[format] = backend_class
    return backend_class


class _Backend(object):
    # A backend can be defined by using the following pattern:
    #
    # @_Backend.export
    # class FooBackend(_Backend):
    #     # override the attributes and methods documented below.

    # The following attributes and methods must be overridden by subclasses.

    # The `FigureCanvas` and `FigureManager` classes must be defined.
    FigureCanvas = None
    FigureManager = None

    # The following methods must be left as None for non-interactive backends.
    # For interactive backends, `trigger_manager_draw` should be a function
    # taking a manager as argument and triggering a canvas draw, and `mainloop`
    # should be a function taking no argument and starting the backend main
    # loop.
    trigger_manager_draw = None
    mainloop = None

    # The following methods will be automatically defined and exported, but
    # can be overridden.

    @classmethod
    def new_figure_manager(cls, num, *args, **kwargs):
        """Create a new figure manager instance.
        """
        # This import needs to happen here due to circular imports.
        from matplotlib.figure import Figure
        fig_cls = kwargs.pop('FigureClass', Figure)
        fig = fig_cls(*args, **kwargs)
        return cls.new_figure_manager_given_figure(num, fig)

    @classmethod
    def new_figure_manager_given_figure(cls, num, figure):
        """Create a new figure manager instance for the given figure.
        """
        canvas = cls.FigureCanvas(figure)
        manager = cls.FigureManager(canvas, num)
        return manager

    @classmethod
    def draw_if_interactive(cls):
        if cls.trigger_manager_draw is not None and is_interactive():
            manager = Gcf.get_active()
            if manager:
                cls.trigger_manager_draw(manager)

    @classmethod
    def show(cls, block=None):
        """Show all figures.

        `show` blocks by calling `mainloop` if *block* is ``True``, or if it
        is ``None`` and we are neither in IPython's ``%pylab`` mode, nor in
        `interactive` mode.
        """
        if cls.mainloop is None:
            return
        managers = Gcf.get_all_fig_managers()
        if not managers:
            return
        for manager in managers:
            manager.show()
        if block is None:
            # Hack: Are we in IPython's pylab mode?
            from matplotlib import pyplot
            try:
                # IPython versions >= 0.10 tack the _needmain attribute onto
                # pyplot.show, and always set it to False, when in %pylab mode.
                ipython_pylab = not pyplot.show._needmain
            except AttributeError:
                ipython_pylab = False
            block = not ipython_pylab and not is_interactive()
            # TODO: The above is a hack to get the WebAgg backend working with
            # ipython's `%pylab` mode until proper integration is implemented.
            if get_backend() == "WebAgg":
                block = True
        if block:
            cls.mainloop()

    # This method is the one actually exporting the required methods.

    @staticmethod
    def export(cls):
        for name in ["FigureCanvas",
                     "FigureManager",
                     "new_figure_manager",
                     "new_figure_manager_given_figure",
                     "draw_if_interactive",
                     "show"]:
            setattr(sys.modules[cls.__module__], name, getattr(cls, name))

        # For back-compatibility, generate a shim `Show` class.

        class Show(ShowBase):
            def mainloop(self):
                return cls.mainloop()

        setattr(sys.modules[cls.__module__], "Show", Show)
        return cls


class ShowBase(_Backend):
    """
    Simple base class to generate a show() callable in backends.

    Subclass must override mainloop() method.
    """

    def __call__(self, block=None):
        return self.show(block=block)


class RendererBase(object):
    """An abstract base class to handle drawing/rendering operations.

    The following methods must be implemented in the backend for full
    functionality (though just implementing :meth:`draw_path` alone would
    give a highly capable backend):

    * :meth:`draw_path`
    * :meth:`draw_image`
    * :meth:`draw_gouraud_triangle`

    The following methods *should* be implemented in the backend for
    optimization reasons:

    * :meth:`draw_text`
    * :meth:`draw_markers`
    * :meth:`draw_path_collection`
    * :meth:`draw_quad_mesh`

    """
    def __init__(self):
        self._texmanager = None
        self._text2path = textpath.TextToPath()

    def open_group(self, s, gid=None):
        """
        Open a grouping element with label *s*. If *gid* is given, use
        *gid* as the id of the group. Is only currently used by
        :mod:`~matplotlib.backends.backend_svg`.
        """

    def close_group(self, s):
        """
        Close a grouping element with label *s*
        Is only currently used by :mod:`~matplotlib.backends.backend_svg`
        """

    def draw_path(self, gc, path, transform, rgbFace=None):
        """
        Draws a :class:`~matplotlib.path.Path` instance using the
        given affine transform.
        """
        raise NotImplementedError

    def draw_markers(self, gc, marker_path, marker_trans, path,
                     trans, rgbFace=None):
        """
        Draws a marker at each of the vertices in path.  This includes
        all vertices, including control points on curves.  To avoid
        that behavior, those vertices should be removed before calling
        this function.

        This provides a fallback implementation of draw_markers that
        makes multiple calls to :meth:`draw_path`.  Some backends may
        want to override this method in order to draw the marker only
        once and reuse it multiple times.

        Parameters
        ----------
        gc : `GraphicsContextBase`
            The graphics context

        marker_trans : `matplotlib.transforms.Transform`
            An affine transform applied to the marker.

        trans : `matplotlib.transforms.Transform`
            An affine transform applied to the path.

        """
        for vertices, codes in path.iter_segments(trans, simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                self.draw_path(gc, marker_path,
                               marker_trans +
                               transforms.Affine2D().translate(x, y),
                               rgbFace)

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offsetTrans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position):
        """
        Draws a collection of paths selecting drawing properties from
        the lists *facecolors*, *edgecolors*, *linewidths*,
        *linestyles* and *antialiaseds*. *offsets* is a list of
        offsets to apply to each of the paths.  The offsets in
        *offsets* are first transformed by *offsetTrans* before being
        applied.  *offset_position* may be either "screen" or "data"
        depending on the space that the offsets are in.

        This provides a fallback implementation of
        :meth:`draw_path_collection` that makes multiple calls to
        :meth:`draw_path`.  Some backends may want to override this in
        order to render each set of path data only once, and then
        reference that path multiple times with the different offsets,
        colors, styles etc.  The generator methods
        :meth:`_iter_collection_raw_paths` and
        :meth:`_iter_collection` are provided to help with (and
        standardize) the implementation across backends.  It is highly
        recommended to use those generators, so that changes to the
        behavior of :meth:`draw_path_collection` can be made globally.
        """
        path_ids = []
        for path, transform in self._iter_collection_raw_paths(
                master_transform, paths, all_transforms):
            path_ids.append((path, transforms.Affine2D(transform)))

        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
                gc, master_transform, all_transforms, path_ids, offsets,
                offsetTrans, facecolors, edgecolors, linewidths, linestyles,
                antialiaseds, urls, offset_position):
            path, transform = path_id
            transform = transforms.Affine2D(
                            transform.get_matrix()).translate(xo, yo)
            self.draw_path(gc0, path, transform, rgbFace)

    def draw_quad_mesh(self, gc, master_transform, meshWidth, meshHeight,
                       coordinates, offsets, offsetTrans, facecolors,
                       antialiased, edgecolors):
        """
        This provides a fallback implementation of
        :meth:`draw_quad_mesh` that generates paths and then calls
        :meth:`draw_path_collection`.
        """

        from matplotlib.collections import QuadMesh
        paths = QuadMesh.convert_mesh_to_paths(
            meshWidth, meshHeight, coordinates)

        if edgecolors is None:
            edgecolors = facecolors
        linewidths = np.array([gc.get_linewidth()], float)

        return self.draw_path_collection(
            gc, master_transform, paths, [], offsets, offsetTrans, facecolors,
            edgecolors, linewidths, [], [antialiased], [None], 'screen')

    def draw_gouraud_triangle(self, gc, points, colors, transform):
        """
        Draw a Gouraud-shaded triangle.

        Parameters
        ----------
        points : array_like, shape=(3, 2)
            Array of (x, y) points for the triangle.

        colors : array_like, shape=(3, 4)
            RGBA colors for each point of the triangle.

        transform : `matplotlib.transforms.Transform`
            An affine transform to apply to the points.

        """
        raise NotImplementedError

    def draw_gouraud_triangles(self, gc, triangles_array, colors_array,
                               transform):
        """
        Draws a series of Gouraud triangles.

        Parameters
        ----------
        points : array_like, shape=(N, 3, 2)
            Array of *N* (x, y) points for the triangles.

        colors : array_like, shape=(N, 3, 4)
            Array of *N* RGBA colors for each point of the triangles.

        transform : `matplotlib.transforms.Transform`
            An affine transform to apply to the points.
        """
        transform = transform.frozen()
        for tri, col in zip(triangles_array, colors_array):
            self.draw_gouraud_triangle(gc, tri, col, transform)

    def _iter_collection_raw_paths(self, master_transform, paths,
                                   all_transforms):
        """
        This is a helper method (along with :meth:`_iter_collection`) to make
        it easier to write a space-efficient :meth:`draw_path_collection`
        implementation in a backend.

        This method yields all of the base path/transform
        combinations, given a master transform, a list of paths and
        list of transforms.

        The arguments should be exactly what is passed in to
        :meth:`draw_path_collection`.

        The backend should take each yielded path and transform and
        create an object that can be referenced (reused) later.
        """
        Npaths = len(paths)
        Ntransforms = len(all_transforms)
        N = max(Npaths, Ntransforms)

        if Npaths == 0:
            return

        transform = transforms.IdentityTransform()
        for i in xrange(N):
            path = paths[i % Npaths]
            if Ntransforms:
                transform = Affine2D(all_transforms[i % Ntransforms])
            yield path, transform + master_transform

    def _iter_collection_uses_per_path(self, paths, all_transforms,
                                       offsets, facecolors, edgecolors):
        """
        Compute how many times each raw path object returned by
        _iter_collection_raw_paths would be used when calling
        _iter_collection. This is intended for the backend to decide
        on the tradeoff between using the paths in-line and storing
        them once and reusing. Rounds up in case the number of uses
        is not the same for every path.
        """
        Npaths = len(paths)
        if Npaths == 0 or (len(facecolors) == 0 and len(edgecolors) == 0):
            return 0
        Npath_ids = max(Npaths, len(all_transforms))
        N = max(Npath_ids, len(offsets))
        return (N + Npath_ids - 1) // Npath_ids

    def _iter_collection(self, gc, master_transform, all_transforms,
                         path_ids, offsets, offsetTrans, facecolors,
                         edgecolors, linewidths, linestyles,
                         antialiaseds, urls, offset_position):
        """
        This is a helper method (along with
        :meth:`_iter_collection_raw_paths`) to make it easier to write
        a space-efficient :meth:`draw_path_collection` implementation in a
        backend.

        This method yields all of the path, offset and graphics
        context combinations to draw the path collection.  The caller
        should already have looped over the results of
        :meth:`_iter_collection_raw_paths` to draw this collection.

        The arguments should be the same as that passed into
        :meth:`draw_path_collection`, with the exception of
        *path_ids*, which is a list of arbitrary objects that the
        backend will use to reference one of the paths created in the
        :meth:`_iter_collection_raw_paths` stage.

        Each yielded result is of the form::

           xo, yo, path_id, gc, rgbFace

        where *xo*, *yo* is an offset; *path_id* is one of the elements of
        *path_ids*; *gc* is a graphics context and *rgbFace* is a color to
        use for filling the path.
        """
        Ntransforms = len(all_transforms)
        Npaths = len(path_ids)
        Noffsets = len(offsets)
        N = max(Npaths, Noffsets)
        Nfacecolors = len(facecolors)
        Nedgecolors = len(edgecolors)
        Nlinewidths = len(linewidths)
        Nlinestyles = len(linestyles)
        Naa = len(antialiaseds)
        Nurls = len(urls)

        if (Nfacecolors == 0 and Nedgecolors == 0) or Npaths == 0:
            return
        if Noffsets:
            toffsets = offsetTrans.transform(offsets)

        gc0 = self.new_gc()
        gc0.copy_properties(gc)

        if Nfacecolors == 0:
            rgbFace = None

        if Nedgecolors == 0:
            gc0.set_linewidth(0.0)

        xo, yo = 0, 0
        for i in xrange(N):
            path_id = path_ids[i % Npaths]
            if Noffsets:
                xo, yo = toffsets[i % Noffsets]
                if offset_position == 'data':
                    if Ntransforms:
                        transform = (
                            Affine2D(all_transforms[i % Ntransforms]) +
                            master_transform)
                    else:
                        transform = master_transform
                    xo, yo = transform.transform_point((xo, yo))
                    xp, yp = transform.transform_point((0, 0))
                    xo = -(xp - xo)
                    yo = -(yp - yo)
            if not (np.isfinite(xo) and np.isfinite(yo)):
                continue
            if Nfacecolors:
                rgbFace = facecolors[i % Nfacecolors]
            if Nedgecolors:
                if Nlinewidths:
                    gc0.set_linewidth(linewidths[i % Nlinewidths])
                if Nlinestyles:
                    gc0.set_dashes(*linestyles[i % Nlinestyles])
                fg = edgecolors[i % Nedgecolors]
                if len(fg) == 4:
                    if fg[3] == 0.0:
                        gc0.set_linewidth(0)
                    else:
                        gc0.set_foreground(fg)
                else:
                    gc0.set_foreground(fg)
            if rgbFace is not None and len(rgbFace) == 4:
                if rgbFace[3] == 0:
                    rgbFace = None
            gc0.set_antialiased(antialiaseds[i % Naa])
            if Nurls:
                gc0.set_url(urls[i % Nurls])

            yield xo, yo, path_id, gc0, rgbFace
        gc0.restore()

    def get_image_magnification(self):
        """
        Get the factor by which to magnify images passed to :meth:`draw_image`.
        Allows a backend to have images at a different resolution to other
        artists.
        """
        return 1.0

    def draw_image(self, gc, x, y, im, transform=None):
        """
        Draw an RGBA image.

        Parameters
        ----------
        gc : `GraphicsContextBase`
            a graphics context with clipping information.

        x : scalar
            the distance in physical units (i.e., dots or pixels) from the left
            hand side of the canvas.

        y : scalar
            the distance in physical units (i.e., dots or pixels) from the
            bottom side of the canvas.

        im : array_like, shape=(N, M, 4), dtype=np.uint8
            An array of RGBA pixels.

        transform : `matplotlib.transforms.Affine2DBase`
            If and only if the concrete backend is written such that
            :meth:`option_scale_image` returns ``True``, an affine
            transformation *may* be passed to :meth:`draw_image`. It takes the
            form of a :class:`~matplotlib.transforms.Affine2DBase` instance.
            The translation vector of the transformation is given in physical
            units (i.e., dots or pixels). Note that the transformation does not
            override `x` and `y`, and has to be applied *before* translating
            the result by `x` and `y` (this can be accomplished by adding `x`
            and `y` to the translation vector defined by `transform`).
        """
        raise NotImplementedError

    def option_image_nocomposite(self):
        """
        override this method for renderers that do not necessarily always
        want to rescale and composite raster images. (like SVG, PDF, or PS)
        """
        return False

    def option_scale_image(self):
        """
        override this method for renderers that support arbitrary affine
        transformations in :meth:`draw_image` (most vector backends).
        """
        return False

    def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
        """
        """
        self._draw_text_as_path(gc, x, y, s, prop, angle, ismath="TeX")

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        """
        Draw the text instance

        Parameters
        ----------
        gc : `GraphicsContextBase`
            the graphics context

        x : scalar
            the x location of the text in display coords

        y : scalar
            the y location of the text baseline in display coords

        s : str
            the text string

        prop : `matplotlib.font_manager.FontProperties`
            font properties

        angle : scalar
            the rotation angle in degrees

        mtext : `matplotlib.text.Text`
            the original text object to be rendered

        Notes
        -----
        **backend implementers note**

        When you are trying to determine if you have gotten your bounding box
        right (which is what enables the text layout/alignment to work
        properly), it helps to change the line in text.py::

            if 0: bbox_artist(self, renderer)

        to if 1, and then the actual bounding box will be plotted along with
        your text.
        """

        self._draw_text_as_path(gc, x, y, s, prop, angle, ismath)

    def _get_text_path_transform(self, x, y, s, prop, angle, ismath):
        """
        return the text path and transform

        Parameters
        ----------
        prop : `matplotlib.font_manager.FontProperties`
          font property

        s : str
          text to be converted

        usetex : bool
          If True, use matplotlib usetex mode.

        ismath : bool
          If True, use mathtext parser. If "TeX", use *usetex* mode.
        """

        text2path = self._text2path
        fontsize = self.points_to_pixels(prop.get_size_in_points())

        if ismath == "TeX":
            verts, codes = text2path.get_text_path(prop, s, ismath=False,
                                                   usetex=True)
        else:
            verts, codes = text2path.get_text_path(prop, s, ismath=ismath,
                                                   usetex=False)

        path = Path(verts, codes)
        angle = np.deg2rad(angle)
        if self.flipy():
            transform = Affine2D().scale(fontsize / text2path.FONT_SCALE,
                                         fontsize / text2path.FONT_SCALE)
            transform = transform.rotate(angle).translate(x, self.height - y)
        else:
            transform = Affine2D().scale(fontsize / text2path.FONT_SCALE,
                                         fontsize / text2path.FONT_SCALE)
            transform = transform.rotate(angle).translate(x, y)

        return path, transform

    def _draw_text_as_path(self, gc, x, y, s, prop, angle, ismath):
        """
        draw the text by converting them to paths using textpath module.

        Parameters
        ----------
        prop : `matplotlib.font_manager.FontProperties`
          font property

        s : str
          text to be converted

        usetex : bool
          If True, use matplotlib usetex mode.

        ismath : bool
          If True, use mathtext parser. If "TeX", use *usetex* mode.
        """
        path, transform = self._get_text_path_transform(
            x, y, s, prop, angle, ismath)
        color = gc.get_rgb()

        gc.set_linewidth(0.0)
        self.draw_path(gc, path, transform, rgbFace=color)

    def get_text_width_height_descent(self, s, prop, ismath):
        """
        Get the width, height, and descent (offset from the bottom
        to the baseline), in display coords, of the string *s* with
        :class:`~matplotlib.font_manager.FontProperties` *prop*
        """
        if ismath == 'TeX':
            # todo: handle props
            size = prop.get_size_in_points()
            texmanager = self._text2path.get_texmanager()
            fontsize = prop.get_size_in_points()
            w, h, d = texmanager.get_text_width_height_descent(
                s, fontsize, renderer=self)
            return w, h, d

        dpi = self.points_to_pixels(72)
        if ismath:
            dims = self._text2path.mathtext_parser.parse(s, dpi, prop)
            return dims[0:3]  # return width, height, descent

        flags = self._text2path._get_hinting_flag()
        font = self._text2path._get_font(prop)
        size = prop.get_size_in_points()
        font.set_size(size, dpi)
        # the width and height of unrotated string
        font.set_text(s, 0.0, flags=flags)
        w, h = font.get_width_height()
        d = font.get_descent()
        w /= 64.0  # convert from subpixels
        h /= 64.0
        d /= 64.0
        return w, h, d

    def flipy(self):
        """
        Return true if y small numbers are top for renderer Is used
        for drawing text (:mod:`matplotlib.text`) and images
        (:mod:`matplotlib.image`) only
        """
        return True

    def get_canvas_width_height(self):
        'return the canvas width and height in display coords'
        return 1, 1

    def get_texmanager(self):
        """
        return the :class:`matplotlib.texmanager.TexManager` instance
        """
        if self._texmanager is None:
            from matplotlib.texmanager import TexManager
            self._texmanager = TexManager()
        return self._texmanager

    def new_gc(self):
        """
        Return an instance of a :class:`GraphicsContextBase`
        """
        return GraphicsContextBase()

    def points_to_pixels(self, points):
        """
        Convert points to display units

        You need to override this function (unless your backend
        doesn't have a dpi, e.g., postscript or svg).  Some imaging
        systems assume some value for pixels per inch::

            points to pixels = points * pixels_per_inch/72.0 * dpi/72.0

        Parameters
        ----------
        points : scalar or array_like
            a float or a numpy array of float

        Returns
        -------
        Points converted to pixels
        """
        return points

    def strip_math(self, s):
        return cbook.strip_math(s)

    def start_rasterizing(self):
        """
        Used in MixedModeRenderer. Switch to the raster renderer.
        """

    def stop_rasterizing(self):
        """
        Used in MixedModeRenderer. Switch back to the vector renderer
        and draw the contents of the raster renderer as an image on
        the vector renderer.
        """

    def start_filter(self):
        """
        Used in AggRenderer. Switch to a temporary renderer for image
        filtering effects.
        """

    def stop_filter(self, filter_func):
        """
        Used in AggRenderer. Switch back to the original renderer.
        The contents of the temporary renderer is processed with the
        *filter_func* and is drawn on the original renderer as an
        image.
        """


class GraphicsContextBase(object):
    """
    An abstract base class that provides color, line styles, etc...
    """

    def __init__(self):
        self._alpha = 1.0
        self._forced_alpha = False  # if True, _alpha overrides A from RGBA
        self._antialiased = 1  # use 0,1 not True, False for extension code
        self._capstyle = 'butt'
        self._cliprect = None
        self._clippath = None
        self._dashes = None, None
        self._joinstyle = 'round'
        self._linestyle = 'solid'
        self._linewidth = 1
        self._rgb = (0.0, 0.0, 0.0, 1.0)
        self._hatch = None
        self._hatch_color = colors.to_rgba(rcParams['hatch.color'])
        self._hatch_linewidth = rcParams['hatch.linewidth']
        self._url = None
        self._gid = None
        self._snap = None
        self._sketch = None

    def copy_properties(self, gc):
        'Copy properties from gc to self'
        self._alpha = gc._alpha
        self._forced_alpha = gc._forced_alpha
        self._antialiased = gc._antialiased
        self._capstyle = gc._capstyle
        self._cliprect = gc._cliprect
        self._clippath = gc._clippath
        self._dashes = gc._dashes
        self._joinstyle = gc._joinstyle
        self._linestyle = gc._linestyle
        self._linewidth = gc._linewidth
        self._rgb = gc._rgb
        self._hatch = gc._hatch
        self._hatch_color = gc._hatch_color
        self._hatch_linewidth = gc._hatch_linewidth
        self._url = gc._url
        self._gid = gc._gid
        self._snap = gc._snap
        self._sketch = gc._sketch

    def restore(self):
        """
        Restore the graphics context from the stack - needed only
        for backends that save graphics contexts on a stack
        """

    def get_alpha(self):
        """
        Return the alpha value used for blending - not supported on
        all backends
        """
        return self._alpha

    def get_antialiased(self):
        "Return true if the object should try to do antialiased rendering"
        return self._antialiased

    def get_capstyle(self):
        """
        Return the capstyle as a string in ('butt', 'round', 'projecting')
        """
        return self._capstyle

    def get_clip_rectangle(self):
        """
        Return the clip rectangle as a :class:`~matplotlib.transforms.Bbox`
        instance
        """
        return self._cliprect

    def get_clip_path(self):
        """
        Return the clip path in the form (path, transform), where path
        is a :class:`~matplotlib.path.Path` instance, and transform is
        an affine transform to apply to the path before clipping.
        """
        if self._clippath is not None:
            return self._clippath.get_transformed_path_and_affine()
        return None, None

    def get_dashes(self):
        """
        Return the dash information as an offset dashlist tuple.

        The dash list is a even size list that gives the ink on, ink
        off in pixels.

        See p107 of to PostScript `BLUEBOOK
        <https://www-cdf.fnal.gov/offline/PostScript/BLUEBOOK.PDF>`_
        for more info.

        Default value is None
        """
        return self._dashes

    def get_forced_alpha(self):
        """
        Return whether the value given by get_alpha() should be used to
        override any other alpha-channel values.
        """
        return self._forced_alpha

    def get_joinstyle(self):
        """
        Return the line join style as one of ('miter', 'round', 'bevel')
        """
        return self._joinstyle

    @cbook.deprecated("2.1")
    def get_linestyle(self):
        """
        Return the linestyle: one of ('solid', 'dashed', 'dashdot',
        'dotted').
        """
        return self._linestyle

    def get_linewidth(self):
        """
        Return the line width in points as a scalar
        """
        return self._linewidth

    def get_rgb(self):
        """
        returns a tuple of three or four floats from 0-1.
        """
        return self._rgb

    def get_url(self):
        """
        returns a url if one is set, None otherwise
        """
        return self._url

    def get_gid(self):
        """
        Return the object identifier if one is set, None otherwise.
        """
        return self._gid

    def get_snap(self):
        """
        returns the snap setting which may be:

          * True: snap vertices to the nearest pixel center

          * False: leave vertices as-is

          * None: (auto) If the path contains only rectilinear line
            segments, round to the nearest pixel center
        """
        return self._snap

    def set_alpha(self, alpha):
        """
        Set the alpha value used for blending - not supported on all backends.
        If ``alpha=None`` (the default), the alpha components of the
        foreground and fill colors will be used to set their respective
        transparencies (where applicable); otherwise, ``alpha`` will override
        them.
        """
        if alpha is not None:
            self._alpha = alpha
            self._forced_alpha = True
        else:
            self._alpha = 1.0
            self._forced_alpha = False
        self.set_foreground(self._rgb, isRGBA=True)

    def set_antialiased(self, b):
        """
        True if object should be drawn with antialiased rendering
        """

        # use 0, 1 to make life easier on extension code trying to read the gc
        if b:
            self._antialiased = 1
        else:
            self._antialiased = 0

    def set_capstyle(self, cs):
        """
        Set the capstyle as a string in ('butt', 'round', 'projecting')
        """
        if cs in ('butt', 'round', 'projecting'):
            self._capstyle = cs
        else:
            raise ValueError('Unrecognized cap style.  Found %s' % cs)

    def set_clip_rectangle(self, rectangle):
        """
        Set the clip rectangle with sequence (left, bottom, width, height)
        """
        self._cliprect = rectangle

    def set_clip_path(self, path):
        """
        Set the clip path and transformation.  Path should be a
        :class:`~matplotlib.transforms.TransformedPath` instance.
        """
        if (path is not None
                and not isinstance(path, transforms.TransformedPath)):
            raise ValueError("Path should be a "
                             "matplotlib.transforms.TransformedPath instance")
        self._clippath = path

    def set_dashes(self, dash_offset, dash_list):
        """
        Set the dash style for the gc.

        Parameters
        ----------
        dash_offset : float
            is the offset (usually 0).

        dash_list : array_like
            specifies the on-off sequence as points.
            ``(None, None)`` specifies a solid line

        """
        if dash_list is not None:
            dl = np.asarray(dash_list)
            if np.any(dl < 0.0):
                raise ValueError(
                    "All values in the dash list must be positive")
        self._dashes = dash_offset, dash_list

    def set_foreground(self, fg, isRGBA=False):
        """
        Set the foreground color.  fg can be a MATLAB format string, a
        html hex color string, an rgb or rgba unit tuple, or a float between 0
        and 1.  In the latter case, grayscale is used.

        If you know fg is rgba, set ``isRGBA=True`` for efficiency.
        """
        if self._forced_alpha and isRGBA:
            self._rgb = fg[:3] + (self._alpha,)
        elif self._forced_alpha:
            self._rgb = colors.to_rgba(fg, self._alpha)
        elif isRGBA:
            self._rgb = fg
        else:
            self._rgb = colors.to_rgba(fg)

    def set_joinstyle(self, js):
        """
        Set the join style to be one of ('miter', 'round', 'bevel')
        """
        if js in ('miter', 'round', 'bevel'):
            self._joinstyle = js
        else:
            raise ValueError('Unrecognized join style.  Found %s' % js)

    def set_linewidth(self, w):
        """
        Set the linewidth in points
        """
        self._linewidth = float(w)

    @cbook.deprecated("2.1")
    def set_linestyle(self, style):
        """
        Set the linestyle to be one of ('solid', 'dashed', 'dashdot',
        'dotted'). These are defined in the rcParams
        `lines.dashed_pattern`, `lines.dashdot_pattern` and
        `lines.dotted_pattern`.  One may also specify customized dash
        styles by providing a tuple of (offset, dash pairs).
        """
        self._linestyle = style

    def set_url(self, url):
        """
        Sets the url for links in compatible backends
        """
        self._url = url

    def set_gid(self, id):
        """
        Sets the id.
        """
        self._gid = id

    def set_snap(self, snap):
        """
        Sets the snap setting which may be:

          * True: snap vertices to the nearest pixel center

          * False: leave vertices as-is

          * None: (auto) If the path contains only rectilinear line
            segments, round to the nearest pixel center
        """
        self._snap = snap

    def set_hatch(self, hatch):
        """
        Sets the hatch style for filling
        """
        self._hatch = hatch

    def get_hatch(self):
        """
        Gets the current hatch style
        """
        return self._hatch

    def get_hatch_path(self, density=6.0):
        """
        Returns a Path for the current hatch.
        """
        hatch = self.get_hatch()
        if hatch is None:
            return None
        return Path.hatch(hatch, density)

    def get_hatch_color(self):
        """
        Gets the color to use for hatching.
        """
        return self._hatch_color

    def set_hatch_color(self, hatch_color):
        """
        sets the color to use for hatching.
        """
        self._hatch_color = hatch_color

    def get_hatch_linewidth(self):
        """
        Gets the linewidth to use for hatching.
        """
        return self._hatch_linewidth

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
             (default 128)

        randomness : float, optional
            The scale factor by which the length is shrunken or
            expanded (default 16)
        """
        self._sketch = (
            None if scale is None
            else (scale, length or 128., randomness or 16.))


class TimerBase(object):
    '''
    A base class for providing timer events, useful for things animations.
    Backends need to implement a few specific methods in order to use their
    own timing mechanisms so that the timer events are integrated into their
    event loops.

    Mandatory functions that must be implemented:

        * `_timer_start`: Contains backend-specific code for starting
          the timer

        * `_timer_stop`: Contains backend-specific code for stopping
          the timer

    Optional overrides:

        * `_timer_set_single_shot`: Code for setting the timer to
          single shot operating mode, if supported by the timer
          object. If not, the `Timer` class itself will store the flag
          and the `_on_timer` method should be overridden to support
          such behavior.

        * `_timer_set_interval`: Code for setting the interval on the
          timer, if there is a method for doing so on the timer
          object.

        * `_on_timer`: This is the internal function that any timer
          object should call, which will handle the task of running
          all callbacks that have been set.

    Attributes
    ----------
    interval : scalar
        The time between timer events in milliseconds. Default is 1000 ms.

    single_shot : bool
        Boolean flag indicating whether this timer should operate as single
        shot (run once and then stop). Defaults to `False`.

    callbacks : List[Tuple[callable, Tuple, Dict]]
        Stores list of (func, args, kwargs) tuples that will be called upon
        timer events. This list can be manipulated directly, or the
        functions `add_callback` and `remove_callback` can be used.

    '''
    def __init__(self, interval=None, callbacks=None):
        #Initialize empty callbacks list and setup default settings if necssary
        if callbacks is None:
            self.callbacks = []
        else:
            self.callbacks = callbacks[:]  # Create a copy

        if interval is None:
            self._interval = 1000
        else:
            self._interval = interval

        self._single = False

        # Default attribute for holding the GUI-specific timer object
        self._timer = None

    def __del__(self):
        'Need to stop timer and possibly disconnect timer.'
        self._timer_stop()

    def start(self, interval=None):
        '''
        Start the timer object. `interval` is optional and will be used
        to reset the timer interval first if provided.
        '''
        if interval is not None:
            self._set_interval(interval)
        self._timer_start()

    def stop(self):
        '''
        Stop the timer.
        '''
        self._timer_stop()

    def _timer_start(self):
        pass

    def _timer_stop(self):
        pass

    def _get_interval(self):
        return self._interval

    def _set_interval(self, interval):
        # Force to int since none of the backends actually support fractional
        # milliseconds, and some error or give warnings.
        interval = int(interval)
        self._interval = interval
        self._timer_set_interval()

    interval = property(_get_interval, _set_interval)

    def _get_single_shot(self):
        return self._single

    def _set_single_shot(self, ss=True):
        self._single = ss
        self._timer_set_single_shot()

    single_shot = property(_get_single_shot, _set_single_shot)

    def add_callback(self, func, *args, **kwargs):
        '''
        Register `func` to be called by timer when the event fires. Any
        additional arguments provided will be passed to `func`.
        '''
        self.callbacks.append((func, args, kwargs))

    def remove_callback(self, func, *args, **kwargs):
        '''
        Remove `func` from list of callbacks. `args` and `kwargs` are optional
        and used to distinguish between copies of the same function registered
        to be called with different arguments.
        '''
        if args or kwargs:
            self.callbacks.remove((func, args, kwargs))
        else:
            funcs = [c[0] for c in self.callbacks]
            if func in funcs:
                self.callbacks.pop(funcs.index(func))

    def _timer_set_interval(self):
        """Used to set interval on underlying timer object."""

    def _timer_set_single_shot(self):
        """Used to set single shot on underlying timer object."""

    def _on_timer(self):
        '''
        Runs all function that have been registered as callbacks. Functions
        can return False (or 0) if they should not be called any more. If there
        are no callbacks, the timer is automatically stopped.
        '''
        for func, args, kwargs in self.callbacks:
            ret = func(*args, **kwargs)
            # docstring above explains why we use `if ret == 0` here,
            # instead of `if not ret`.
            # This will also catch `ret == False` as `False == 0`
            # but does not annoy the linters
            # https://docs.python.org/3/library/stdtypes.html#boolean-values
            if ret == 0:
                self.callbacks.remove((func, args, kwargs))

        if len(self.callbacks) == 0:
            self.stop()


class Event(object):
    """
    A matplotlib event.  Attach additional attributes as defined in
    :meth:`FigureCanvasBase.mpl_connect`.  The following attributes
    are defined and shown with their default values

    Attributes
    ----------
    name : str
        the event name

    canvas : `FigureCanvasBase`
        the backend-specific canvas instance generating the event

    guiEvent
        the GUI event that triggered the matplotlib event

    """
    def __init__(self, name, canvas, guiEvent=None):
        self.name = name
        self.canvas = canvas
        self.guiEvent = guiEvent


@cbook.deprecated("2.1")
class IdleEvent(Event):
    """
    An event triggered by the GUI backend when it is idle -- useful
    for passive animation
    """


class DrawEvent(Event):
    """
    An event triggered by a draw operation on the canvas

    In most backends callbacks subscribed to this callback will be
    fired after the rendering is complete but before the screen is
    updated.  Any extra artists drawn to the canvas's renderer will
    be reflected without an explicit call to ``blit``.

    .. warning ::

       Calling ``canvas.draw`` and ``canvas.blit`` in these callbacks may
       not be safe with all backends and may cause infinite recursion.

    In addition to the :class:`Event` attributes, the following event
    attributes are defined:

    Attributes
    ----------
    renderer : `RendererBase`
        the renderer for the draw event

    """
    def __init__(self, name, canvas, renderer):
        Event.__init__(self, name, canvas)
        self.renderer = renderer


class ResizeEvent(Event):
    """
    An event triggered by a canvas resize

    In addition to the :class:`Event` attributes, the following event
    attributes are defined:

    Attributes
    ----------
    width : scalar
        width of the canvas in pixels

    height : scalar
        height of the canvas in pixels

    """
    def __init__(self, name, canvas):
        Event.__init__(self, name, canvas)
        self.width, self.height = canvas.get_width_height()


class CloseEvent(Event):
    """
    An event triggered by a figure being closed

    """
    def __init__(self, name, canvas, guiEvent=None):
        Event.__init__(self, name, canvas, guiEvent)


class LocationEvent(Event):
    """
    An event that has a screen location

    The following additional attributes are defined and shown with
    their default values.

    In addition to the :class:`Event` attributes, the following
    event attributes are defined:

    Attributes
    ----------
    x : scalar
        x position - pixels from left of canvas

    y : scalar
        y position - pixels from bottom of canvas

    inaxes : bool
        the :class:`~matplotlib.axes.Axes` instance if mouse is over axes

    xdata : scalar
        x coord of mouse in data coords

    ydata : scalar
        y coord of mouse in data coords

    """
    x = None       # x position - pixels from left of canvas
    y = None       # y position - pixels from right of canvas
    inaxes = None  # the Axes instance if mouse us over axes
    xdata = None   # x coord of mouse in data coords
    ydata = None   # y coord of mouse in data coords

    # the last event that was triggered before this one
    lastevent = None

    def __init__(self, name, canvas, x, y, guiEvent=None):
        """
        *x*, *y* in figure coords, 0,0 = bottom, left
        """
        Event.__init__(self, name, canvas, guiEvent=guiEvent)
        self.x = x
        self.y = y

        if x is None or y is None:
            # cannot check if event was in axes if no x,y info
            self.inaxes = None
            self._update_enter_leave()
            return

        # Find all axes containing the mouse
        if self.canvas.mouse_grabber is None:
            axes_list = [a for a in self.canvas.figure.get_axes()
                         if a.in_axes(self)]
        else:
            axes_list = [self.canvas.mouse_grabber]

        if axes_list:
            self.inaxes = cbook._topmost_artist(axes_list)
            try:
                trans = self.inaxes.transData.inverted()
                xdata, ydata = trans.transform_point((x, y))
            except ValueError:
                self.xdata = None
                self.ydata = None
            else:
                self.xdata = xdata
                self.ydata = ydata
        else:
            self.inaxes = None

        self._update_enter_leave()

    def _update_enter_leave(self):
        'process the figure/axes enter leave events'
        if LocationEvent.lastevent is not None:
            last = LocationEvent.lastevent
            if last.inaxes != self.inaxes:
                # process axes enter/leave events
                try:
                    if last.inaxes is not None:
                        last.canvas.callbacks.process('axes_leave_event', last)
                except:
                    pass
                    # See ticket 2901582.
                    # I think this is a valid exception to the rule
                    # against catching all exceptions; if anything goes
                    # wrong, we simply want to move on and process the
                    # current event.
                if self.inaxes is not None:
                    self.canvas.callbacks.process('axes_enter_event', self)

        else:
            # process a figure enter event
            if self.inaxes is not None:
                self.canvas.callbacks.process('axes_enter_event', self)

        LocationEvent.lastevent = self


class MouseEvent(LocationEvent):
    """
    A mouse event ('button_press_event',
                   'button_release_event',
                   'scroll_event',
                   'motion_notify_event').

    In addition to the :class:`Event` and :class:`LocationEvent`
    attributes, the following attributes are defined:

    Attributes
    ----------
    button : None, scalar, or str
        button pressed None, 1, 2, 3, 'up', 'down' (up and down are used
        for scroll events).  Note that in the nbagg backend, both the
        middle and right clicks return 3 since right clicking will bring
        up the context menu in some browsers.

    key : None, or str
        the key depressed when the mouse event triggered (see
        :class:`KeyEvent`)

    step : scalar
        number of scroll steps (positive for 'up', negative for 'down')

    Examples
    --------
    Usage::

        def on_press(event):
            print('you pressed', event.button, event.xdata, event.ydata)

        cid = fig.canvas.mpl_connect('button_press_event', on_press)

    """
    x = None         # x position - pixels from left of canvas
    y = None         # y position - pixels from right of canvas
    button = None    # button pressed None, 1, 2, 3
    dblclick = None  # whether or not the event is the result of a double click
    inaxes = None    # the Axes instance if mouse us over axes
    xdata = None     # x coord of mouse in data coords
    ydata = None     # y coord of mouse in data coords
    step = None      # scroll steps for scroll events

    def __init__(self, name, canvas, x, y, button=None, key=None,
                 step=0, dblclick=False, guiEvent=None):
        """
        x, y in figure coords, 0,0 = bottom, left
        button pressed None, 1, 2, 3, 'up', 'down'
        """
        LocationEvent.__init__(self, name, canvas, x, y, guiEvent=guiEvent)
        self.button = button
        self.key = key
        self.step = step
        self.dblclick = dblclick

    def __str__(self):
        return ("MPL MouseEvent: xy=(%d,%d) xydata=(%s,%s) button=%s " +
                "dblclick=%s inaxes=%s") % (self.x, self.y, self.xdata,
                                            self.ydata, self.button,
                                            self.dblclick, self.inaxes)


class PickEvent(Event):
    """
    a pick event, fired when the user picks a location on the canvas
    sufficiently close to an artist.

    Attrs: all the :class:`Event` attributes plus

    Attributes
    ----------
    mouseevent : `MouseEvent`
        the mouse event that generated the pick

    artist : `matplotlib.artist.Artist`
        the picked artist

    other
        extra class dependent attrs -- e.g., a
        :class:`~matplotlib.lines.Line2D` pick may define different
        extra attributes than a
        :class:`~matplotlib.collections.PatchCollection` pick event

    Examples
    --------
    Usage::

        ax.plot(np.rand(100), 'o', picker=5)  # 5 points tolerance

        def on_pick(event):
            line = event.artist
            xdata, ydata = line.get_data()
            ind = event.ind
            print('on pick line:', np.array([xdata[ind], ydata[ind]]).T)

        cid = fig.canvas.mpl_connect('pick_event', on_pick)

    """
    def __init__(self, name, canvas, mouseevent, artist,
                 guiEvent=None, **kwargs):
        Event.__init__(self, name, canvas, guiEvent)
        self.mouseevent = mouseevent
        self.artist = artist
        self.__dict__.update(kwargs)


class KeyEvent(LocationEvent):
    """
    A key event (key press, key release).

    Attach additional attributes as defined in
    :meth:`FigureCanvasBase.mpl_connect`.

    In addition to the :class:`Event` and :class:`LocationEvent`
    attributes, the following attributes are defined:

    Attributes
    ----------
    key : None or str
        the key(s) pressed. Could be **None**, a single case sensitive ascii
        character ("g", "G", "#", etc.), a special key
        ("control", "shift", "f1", "up", etc.) or a
        combination of the above (e.g., "ctrl+alt+g", "ctrl+alt+G").

    Notes
    -----
    Modifier keys will be prefixed to the pressed key and will be in the order
    "ctrl", "alt", "super". The exception to this rule is when the pressed key
    is itself a modifier key, therefore "ctrl+alt" and "alt+control" can both
    be valid key values.

    Examples
    --------
    Usage::

        def on_key(event):
            print('you pressed', event.key, event.xdata, event.ydata)

        cid = fig.canvas.mpl_connect('key_press_event', on_key)

    """
    def __init__(self, name, canvas, key, x=0, y=0, guiEvent=None):
        LocationEvent.__init__(self, name, canvas, x, y, guiEvent=guiEvent)
        self.key = key


class FigureCanvasBase(object):
    """
    The canvas the figure renders into.

    Public attributes

    Attributes
    ----------
    figure : `matplotlib.figure.Figure`
        A high-level figure instance

    """
    events = [
        'resize_event',
        'draw_event',
        'key_press_event',
        'key_release_event',
        'button_press_event',
        'button_release_event',
        'scroll_event',
        'motion_notify_event',
        'pick_event',
        'idle_event',
        'figure_enter_event',
        'figure_leave_event',
        'axes_enter_event',
        'axes_leave_event',
        'close_event'
    ]

    supports_blit = True
    fixed_dpi = None

    filetypes = _default_filetypes
    if _has_pil:
        # JPEG support
        register_backend('jpg', 'matplotlib.backends.backend_agg',
                         'Joint Photographic Experts Group')
        register_backend('jpeg', 'matplotlib.backends.backend_agg',
                         'Joint Photographic Experts Group')
        # TIFF support
        register_backend('tif', 'matplotlib.backends.backend_agg',
                         'Tagged Image File Format')
        register_backend('tiff', 'matplotlib.backends.backend_agg',
                         'Tagged Image File Format')

    def __init__(self, figure):
        self._is_idle_drawing = True
        self._is_saving = False
        figure.set_canvas(self)
        self.figure = figure
        # a dictionary from event name to a dictionary that maps cid->func
        self.callbacks = cbook.CallbackRegistry()
        self.widgetlock = widgets.LockDraw()
        self._button = None  # the button pressed
        self._key = None  # the key pressed
        self._lastx, self._lasty = None, None
        self.button_pick_id = self.mpl_connect('button_press_event', self.pick)
        self.scroll_pick_id = self.mpl_connect('scroll_event', self.pick)
        self.mouse_grabber = None  # the axes currently grabbing mouse
        self.toolbar = None  # NavigationToolbar2 will set me
        self._is_idle_drawing = False

    @contextmanager
    def _idle_draw_cntx(self):
        self._is_idle_drawing = True
        yield
        self._is_idle_drawing = False

    def is_saving(self):
        """
        Returns whether the renderer is in the process of saving
        to a file, rather than rendering for an on-screen buffer.
        """
        return self._is_saving

    @cbook.deprecated("2.2")
    def onRemove(self, ev):
        """
        Mouse event processor which removes the top artist
        under the cursor.  Connect this to the 'mouse_press_event'
        using::

            canvas.mpl_connect('mouse_press_event',canvas.onRemove)
        """
        # Find the top artist under the cursor
        under = cbook._topmost_artist(self.figure.hitlist(ev))
        h = None
        if under:
            h = under[-1]

        # Try deleting that artist, or its parent if you
        # can't delete the artist
        while h:
            if h.remove():
                self.draw_idle()
                break
            parent = None
            for p in under:
                if h in p.get_children():
                    parent = p
                    break
            h = parent

    def pick(self, mouseevent):
        if not self.widgetlock.locked():
            self.figure.pick(mouseevent)

    def blit(self, bbox=None):
        """Blit the canvas in bbox (default entire canvas)."""

    def resize(self, w, h):
        """Set the canvas size in pixels."""

    def draw_event(self, renderer):
        """Pass a `DrawEvent` to all functions connected to ``draw_event``."""
        s = 'draw_event'
        event = DrawEvent(s, self, renderer)
        self.callbacks.process(s, event)

    def resize_event(self):
        """Pass a `ResizeEvent` to all functions connected to ``resize_event``.
        """
        s = 'resize_event'
        event = ResizeEvent(s, self)
        self.callbacks.process(s, event)
        self.draw_idle()

    def close_event(self, guiEvent=None):
        """Pass a `CloseEvent` to all functions connected to ``close_event``.
        """
        s = 'close_event'
        try:
            event = CloseEvent(s, self, guiEvent=guiEvent)
            self.callbacks.process(s, event)
        except (TypeError, AttributeError):
            pass
            # Suppress the TypeError when the python session is being killed.
            # It may be that a better solution would be a mechanism to
            # disconnect all callbacks upon shutdown.
            # AttributeError occurs on OSX with qt4agg upon exiting
            # with an open window; 'callbacks' attribute no longer exists.

    def key_press_event(self, key, guiEvent=None):
        """Pass a `KeyEvent` to all functions connected to ``key_press_event``.
        """
        self._key = key
        s = 'key_press_event'
        event = KeyEvent(
            s, self, key, self._lastx, self._lasty, guiEvent=guiEvent)
        self.callbacks.process(s, event)

    def key_release_event(self, key, guiEvent=None):
        """
        Pass a `KeyEvent` to all functions connected to ``key_release_event``.
        """
        s = 'key_release_event'
        event = KeyEvent(
            s, self, key, self._lastx, self._lasty, guiEvent=guiEvent)
        self.callbacks.process(s, event)
        self._key = None

    def pick_event(self, mouseevent, artist, **kwargs):
        """
        This method will be called by artists who are picked and will
        fire off :class:`PickEvent` callbacks registered listeners
        """
        s = 'pick_event'
        event = PickEvent(s, self, mouseevent, artist,
                          guiEvent=mouseevent.guiEvent,
                          **kwargs)
        self.callbacks.process(s, event)

    def scroll_event(self, x, y, step, guiEvent=None):
        """
        Backend derived classes should call this function on any
        scroll wheel event.  x,y are the canvas coords: 0,0 is lower,
        left.  button and key are as defined in MouseEvent.

        This method will be call all functions connected to the
        'scroll_event' with a :class:`MouseEvent` instance.
        """
        if step >= 0:
            self._button = 'up'
        else:
            self._button = 'down'
        s = 'scroll_event'
        mouseevent = MouseEvent(s, self, x, y, self._button, self._key,
                                step=step, guiEvent=guiEvent)
        self.callbacks.process(s, mouseevent)

    def button_press_event(self, x, y, button, dblclick=False, guiEvent=None):
        """
        Backend derived classes should call this function on any mouse
        button press.  x,y are the canvas coords: 0,0 is lower, left.
        button and key are as defined in :class:`MouseEvent`.

        This method will be call all functions connected to the
        'button_press_event' with a :class:`MouseEvent` instance.
        """
        self._button = button
        s = 'button_press_event'
        mouseevent = MouseEvent(s, self, x, y, button, self._key,
                                dblclick=dblclick, guiEvent=guiEvent)
        self.callbacks.process(s, mouseevent)

    def button_release_event(self, x, y, button, guiEvent=None):
        """
        Backend derived classes should call this function on any mouse
        button release.

        This method will call all functions connected to the
        'button_release_event' with a :class:`MouseEvent` instance.

        Parameters
        ----------
        x : scalar
            the canvas coordinates where 0=left

        y : scalar
            the canvas coordinates where 0=bottom

        guiEvent
            the native UI event that generated the mpl event

        """
        s = 'button_release_event'
        event = MouseEvent(s, self, x, y, button, self._key, guiEvent=guiEvent)
        self.callbacks.process(s, event)
        self._button = None

    def motion_notify_event(self, x, y, guiEvent=None):
        """
        Backend derived classes should call this function on any
        motion-notify-event.

        This method will call all functions connected to the
        'motion_notify_event' with a :class:`MouseEvent` instance.

        Parameters
        ----------
        x : scalar
            the canvas coordinates where 0=left

        y : scalar
            the canvas coordinates where 0=bottom

        guiEvent
            the native UI event that generated the mpl event

        """
        self._lastx, self._lasty = x, y
        s = 'motion_notify_event'
        event = MouseEvent(s, self, x, y, self._button, self._key,
                           guiEvent=guiEvent)
        self.callbacks.process(s, event)

    def leave_notify_event(self, guiEvent=None):
        """
        Backend derived classes should call this function when leaving
        canvas

        Parameters
        ----------
        guiEvent
            the native UI event that generated the mpl event

        """

        self.callbacks.process('figure_leave_event', LocationEvent.lastevent)
        LocationEvent.lastevent = None
        self._lastx, self._lasty = None, None

    def enter_notify_event(self, guiEvent=None, xy=None):
        """
        Backend derived classes should call this function when entering
        canvas

        Parameters
        ----------
        guiEvent
            the native UI event that generated the mpl event
        xy : tuple of 2 scalars
            the coordinate location of the pointer when the canvas is
            entered

        """
        if xy is not None:
            x, y = xy
            self._lastx, self._lasty = x, y

        event = Event('figure_enter_event', self, guiEvent)
        self.callbacks.process('figure_enter_event', event)

    @cbook.deprecated("2.1")
    def idle_event(self, guiEvent=None):
        """Called when GUI is idle."""
        s = 'idle_event'
        event = IdleEvent(s, self, guiEvent=guiEvent)
        self.callbacks.process(s, event)

    def grab_mouse(self, ax):
        """
        Set the child axes which are currently grabbing the mouse events.
        Usually called by the widgets themselves.
        It is an error to call this if the mouse is already grabbed by
        another axes.
        """
        if self.mouse_grabber not in (None, ax):
            raise RuntimeError("Another Axes already grabs mouse input")
        self.mouse_grabber = ax

    def release_mouse(self, ax):
        """
        Release the mouse grab held by the axes, ax.
        Usually called by the widgets.
        It is ok to call this even if you ax doesn't have the mouse
        grab currently.
        """
        if self.mouse_grabber is ax:
            self.mouse_grabber = None

    def draw(self, *args, **kwargs):
        """Render the :class:`~matplotlib.figure.Figure`."""

    def draw_idle(self, *args, **kwargs):
        """
        :meth:`draw` only if idle; defaults to draw but backends can override
        """
        if not self._is_idle_drawing:
            with self._idle_draw_cntx():
                self.draw(*args, **kwargs)

    def draw_cursor(self, event):
        """
        Draw a cursor in the event.axes if inaxes is not None.  Use
        native GUI drawing for efficiency if possible
        """

    def get_width_height(self):
        """
        Return the figure width and height in points or pixels
        (depending on the backend), truncated to integers
        """
        return int(self.figure.bbox.width), int(self.figure.bbox.height)

    @classmethod
    def get_supported_filetypes(cls):
        """Return dict of savefig file formats supported by this backend"""
        return cls.filetypes

    @classmethod
    def get_supported_filetypes_grouped(cls):
        """Return a dict of savefig file formats supported by this backend,
        where the keys are a file type name, such as 'Joint Photographic
        Experts Group', and the values are a list of filename extensions used
        for that filetype, such as ['jpg', 'jpeg']."""
        groupings = {}
        for ext, name in six.iteritems(cls.filetypes):
            groupings.setdefault(name, []).append(ext)
            groupings[name].sort()
        return groupings

    def _get_output_canvas(self, fmt):
        """
        Return a canvas suitable for saving figures to a specified file format.

        If necessary, this function will switch to a registered backend that
        supports the format.
        """
        method_name = 'print_%s' % fmt
        # Return the current canvas if it supports the requested format.
        if hasattr(self, method_name):
            return self
        # Return a default canvas for the requested format, if it exists.
        canvas_class = get_registered_canvas_class(fmt)
        if canvas_class:
            return self.switch_backends(canvas_class)
        # Else report error for unsupported format.
        raise ValueError(
            "Format {!r} is not supported (supported formats: {})"
            .format(fmt, ", ".join(sorted(self.get_supported_filetypes()))))

    def print_figure(self, filename, dpi=None, facecolor=None, edgecolor=None,
                     orientation='portrait', format=None, **kwargs):
        """
        Render the figure to hardcopy. Set the figure patch face and edge
        colors.  This is useful because some of the GUIs have a gray figure
        face color background and you'll probably want to override this on
        hardcopy.

        Parameters
        ----------
        filename
            can also be a file object on image backends

        orientation : {'landscape', 'portrait'}, optional
            only currently applies to PostScript printing.

        dpi : scalar, optional
            the dots per inch to save the figure in; if None, use savefig.dpi

        facecolor : color spec or None, optional
            the facecolor of the figure; if None, defaults to savefig.facecolor

        edgecolor : color spec or None, optional
            the edgecolor of the figure; if None, defaults to savefig.edgecolor

        format : str, optional
            when set, forcibly set the file format to save to

        bbox_inches : str or `~matplotlib.transforms.Bbox`, optional
            Bbox in inches. Only the given portion of the figure is
            saved. If 'tight', try to figure out the tight bbox of
            the figure. If None, use savefig.bbox

        pad_inches : scalar, optional
            Amount of padding around the figure when bbox_inches is
            'tight'. If None, use savefig.pad_inches

        bbox_extra_artists : list of `~matplotlib.artist.Artist`, optional
            A list of extra artists that will be considered when the
            tight bbox is calculated.

        """
        self._is_saving = True
        # Remove the figure manager, if any, to avoid resizing the GUI widget.
        # Having *no* manager and a *None* manager are currently different (see
        # Figure.show); should probably be normalized to None at some point.
        _no_manager = object()
        if hasattr(self, 'manager'):
            manager = self.manager
            del self.manager
        else:
            manager = _no_manager

        if format is None:
            # get format from filename, or from backend's default filetype
            if isinstance(filename, getattr(os, "PathLike", ())):
                filename = os.fspath(filename)
            if isinstance(filename, six.string_types):
                format = os.path.splitext(filename)[1][1:]
            if format is None or format == '':
                format = self.get_default_filetype()
                if isinstance(filename, six.string_types):
                    filename = filename.rstrip('.') + '.' + format
        format = format.lower()

        # get canvas object and print method for format
        canvas = self._get_output_canvas(format)
        print_method = getattr(canvas, 'print_%s' % format)

        if dpi is None:
            dpi = rcParams['savefig.dpi']

        if dpi == 'figure':
            dpi = getattr(self.figure, '_original_dpi', self.figure.dpi)

        if facecolor is None:
            facecolor = rcParams['savefig.facecolor']
        if edgecolor is None:
            edgecolor = rcParams['savefig.edgecolor']

        origDPI = self.figure.dpi
        origfacecolor = self.figure.get_facecolor()
        origedgecolor = self.figure.get_edgecolor()

        self.figure.dpi = dpi
        self.figure.set_facecolor(facecolor)
        self.figure.set_edgecolor(edgecolor)

        bbox_inches = kwargs.pop("bbox_inches", None)
        if bbox_inches is None:
            bbox_inches = rcParams['savefig.bbox']

        if bbox_inches:
            # call adjust_bbox to save only the given area
            if bbox_inches == "tight":
                # When bbox_inches == "tight", it saves the figure twice.  The
                # first save command (to a BytesIO) is just to estimate the
                # bounding box of the figure.
                result = print_method(
                    io.BytesIO(),
                    dpi=dpi,
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    orientation=orientation,
                    dryrun=True,
                    **kwargs)
                renderer = self.figure._cachedRenderer
                bbox_inches = self.figure.get_tightbbox(renderer)

                bbox_artists = kwargs.pop("bbox_extra_artists", None)
                if bbox_artists is None:
                    bbox_artists = self.figure.get_default_bbox_extra_artists()

                bbox_filtered = []
                for a in bbox_artists:
                    bbox = a.get_window_extent(renderer)
                    if a.get_clip_on():
                        clip_box = a.get_clip_box()
                        if clip_box is not None:
                            bbox = Bbox.intersection(bbox, clip_box)
                        clip_path = a.get_clip_path()
                        if clip_path is not None and bbox is not None:
                            clip_path = clip_path.get_fully_transformed_path()
                            bbox = Bbox.intersection(bbox,
                                                     clip_path.get_extents())
                    if bbox is not None and (bbox.width != 0 or
                                             bbox.height != 0):
                        bbox_filtered.append(bbox)

                if bbox_filtered:
                    _bbox = Bbox.union(bbox_filtered)
                    trans = Affine2D().scale(1.0 / self.figure.dpi)
                    bbox_extra = TransformedBbox(_bbox, trans)
                    bbox_inches = Bbox.union([bbox_inches, bbox_extra])

                pad = kwargs.pop("pad_inches", None)
                if pad is None:
                    pad = rcParams['savefig.pad_inches']

                bbox_inches = bbox_inches.padded(pad)

            restore_bbox = tight_bbox.adjust_bbox(self.figure, bbox_inches,
                                                  canvas.fixed_dpi)

            _bbox_inches_restore = (bbox_inches, restore_bbox)
        else:
            _bbox_inches_restore = None

        try:
            result = print_method(
                filename,
                dpi=dpi,
                facecolor=facecolor,
                edgecolor=edgecolor,
                orientation=orientation,
                bbox_inches_restore=_bbox_inches_restore,
                **kwargs)
        finally:
            if bbox_inches and restore_bbox:
                restore_bbox()

            self.figure.dpi = origDPI
            self.figure.set_facecolor(origfacecolor)
            self.figure.set_edgecolor(origedgecolor)
            self.figure.set_canvas(self)
            if manager is not _no_manager:
                self.manager = manager
            self._is_saving = False
        return result

    @classmethod
    def get_default_filetype(cls):
        """
        Get the default savefig file format as specified in rcParam
        ``savefig.format``. Returned string excludes period. Overridden
        in backends that only support a single file type.
        """
        return rcParams['savefig.format']

    def get_window_title(self):
        """
        Get the title text of the window containing the figure.
        Return None if there is no window (e.g., a PS backend).
        """
        if hasattr(self, "manager"):
            return self.manager.get_window_title()

    def set_window_title(self, title):
        """
        Set the title text of the window containing the figure.  Note that
        this has no effect if there is no window (e.g., a PS backend).
        """
        if hasattr(self, "manager"):
            self.manager.set_window_title(title)

    def get_default_filename(self):
        """
        Return a string, which includes extension, suitable for use as
        a default filename.
        """
        default_basename = self.get_window_title() or 'image'
        default_basename = default_basename.replace(' ', '_')
        default_filetype = self.get_default_filetype()
        default_filename = default_basename + '.' + default_filetype

        save_dir = os.path.expanduser(rcParams['savefig.directory'])

        # ensure non-existing filename in save dir
        i = 1
        while os.path.isfile(os.path.join(save_dir, default_filename)):
            # attach numerical count to basename
            default_filename = (
                '{}-{}.{}'.format(default_basename, i, default_filetype))
            i += 1

        return default_filename

    def switch_backends(self, FigureCanvasClass):
        """
        Instantiate an instance of FigureCanvasClass

        This is used for backend switching, e.g., to instantiate a
        FigureCanvasPS from a FigureCanvasGTK.  Note, deep copying is
        not done, so any changes to one of the instances (e.g., setting
        figure size or line props), will be reflected in the other
        """
        newCanvas = FigureCanvasClass(self.figure)
        newCanvas._is_saving = self._is_saving
        return newCanvas

    def mpl_connect(self, s, func):
        """
        Connect event with string *s* to *func*.  The signature of *func* is::

          def func(event)

        where event is a :class:`matplotlib.backend_bases.Event`.  The
        following events are recognized

        - 'button_press_event'
        - 'button_release_event'
        - 'draw_event'
        - 'key_press_event'
        - 'key_release_event'
        - 'motion_notify_event'
        - 'pick_event'
        - 'resize_event'
        - 'scroll_event'
        - 'figure_enter_event',
        - 'figure_leave_event',
        - 'axes_enter_event',
        - 'axes_leave_event'
        - 'close_event'

        For the location events (button and key press/release), if the
        mouse is over the axes, the variable ``event.inaxes`` will be
        set to the :class:`~matplotlib.axes.Axes` the event occurs is
        over, and additionally, the variables ``event.xdata`` and
        ``event.ydata`` will be defined.  This is the mouse location
        in data coords.  See
        :class:`~matplotlib.backend_bases.KeyEvent` and
        :class:`~matplotlib.backend_bases.MouseEvent` for more info.

        Return value is a connection id that can be used with
        :meth:`~matplotlib.backend_bases.Event.mpl_disconnect`.

        Examples
        --------
        Usage::

            def on_press(event):
                print('you pressed', event.button, event.xdata, event.ydata)

            cid = canvas.mpl_connect('button_press_event', on_press)

        """
        if s == 'idle_event':
            cbook.warn_deprecated(1.5,
                "idle_event is only implemented for the wx backend, and will "
                "be removed in matplotlib 2.1. Use the animations module "
                "instead.")

        return self.callbacks.connect(s, func)

    def mpl_disconnect(self, cid):
        """
        Disconnect callback id cid

        Examples
        --------
        Usage::

            cid = canvas.mpl_connect('button_press_event', on_press)
            #...later
            canvas.mpl_disconnect(cid)
        """
        return self.callbacks.disconnect(cid)

    def new_timer(self, *args, **kwargs):
        """
        Creates a new backend-specific subclass of
        :class:`backend_bases.Timer`. This is useful for getting periodic
        events through the backend's native event loop. Implemented only for
        backends with GUIs.

        Other Parameters
        ----------------
        interval : scalar
            Timer interval in milliseconds

        callbacks : List[Tuple[callable, Tuple, Dict]]
            Sequence of (func, args, kwargs) where ``func(*args, **kwargs)``
            will be executed by the timer every *interval*.

            callbacks which return ``False`` or ``0`` will be removed from the
            timer.

        Examples
        --------

        >>> timer = fig.canvas.new_timer(callbacks=[(f1, (1, ), {'a': 3}),])

        """
        return TimerBase(*args, **kwargs)

    def flush_events(self):
        """Flush the GUI events for the figure.

        Interactive backends need to reimplement this method.
        """

    def start_event_loop(self, timeout=0):
        """Start a blocking event loop.

        Such an event loop is used by interactive functions, such as `ginput`
        and `waitforbuttonpress`, to wait for events.

        The event loop blocks until a callback function triggers
        `stop_event_loop`, or *timeout* is reached.

        If *timeout* is negative, never timeout.

        Only interactive backends need to reimplement this method and it relies
        on `flush_events` being properly implemented.

        Interactive backends should implement this in a more native way.
        """
        if timeout <= 0:
            timeout = np.inf
        timestep = 0.01
        counter = 0
        self._looping = True
        while self._looping and counter * timestep < timeout:
            self.flush_events()
            time.sleep(timestep)
            counter += 1

    def stop_event_loop(self):
        """Stop the current blocking event loop.

        Interactive backends need to reimplement this to match
        `start_event_loop`
        """
        self._looping = False

    start_event_loop_default = cbook.deprecated(
        "2.1", name="start_event_loop_default")(start_event_loop)
    stop_event_loop_default = cbook.deprecated(
        "2.1", name="stop_event_loop_default")(stop_event_loop)


def key_press_handler(event, canvas, toolbar=None):
    """
    Implement the default mpl key bindings for the canvas and toolbar
    described at :ref:`key-event-handling`

    Parameters
    ----------
    event : :class:`KeyEvent`
        a key press/release event
    canvas : :class:`FigureCanvasBase`
        the backend-specific canvas instance
    toolbar : :class:`NavigationToolbar2`
        the navigation cursor toolbar

    """
    # these bindings happen whether you are over an axes or not

    if event.key is None:
        return

    # Load key-mappings from your matplotlibrc file.
    fullscreen_keys = rcParams['keymap.fullscreen']
    home_keys = rcParams['keymap.home']
    back_keys = rcParams['keymap.back']
    forward_keys = rcParams['keymap.forward']
    pan_keys = rcParams['keymap.pan']
    zoom_keys = rcParams['keymap.zoom']
    save_keys = rcParams['keymap.save']
    quit_keys = rcParams['keymap.quit']
    grid_keys = rcParams['keymap.grid']
    grid_minor_keys = rcParams['keymap.grid_minor']
    toggle_yscale_keys = rcParams['keymap.yscale']
    toggle_xscale_keys = rcParams['keymap.xscale']
    all_keys = rcParams['keymap.all_axes']

    # toggle fullscreen mode ('f', 'ctrl + f')
    if event.key in fullscreen_keys:
        try:
            canvas.manager.full_screen_toggle()
        except AttributeError:
            pass

    # quit the figure (default key 'ctrl+w')
    if event.key in quit_keys:
        Gcf.destroy_fig(canvas.figure)

    if toolbar is not None:
        # home or reset mnemonic  (default key 'h', 'home' and 'r')
        if event.key in home_keys:
            toolbar.home()
        # forward / backward keys to enable left handed quick navigation
        # (default key for backward: 'left', 'backspace' and 'c')
        elif event.key in back_keys:
            toolbar.back()
        # (default key for forward: 'right' and 'v')
        elif event.key in forward_keys:
            toolbar.forward()
        # pan mnemonic (default key 'p')
        elif event.key in pan_keys:
            toolbar.pan()
            toolbar._set_cursor(event)
        # zoom mnemonic (default key 'o')
        elif event.key in zoom_keys:
            toolbar.zoom()
            toolbar._set_cursor(event)
        # saving current figure (default key 's')
        elif event.key in save_keys:
            toolbar.save_figure()

    if event.inaxes is None:
        return

    # these bindings require the mouse to be over an axes to trigger
    def _get_uniform_gridstate(ticks):
        # Return True/False if all grid lines are on or off, None if they are
        # not all in the same state.
        if all(tick.gridOn for tick in ticks):
            return True
        elif not any(tick.gridOn for tick in ticks):
            return False
        else:
            return None

    ax = event.inaxes
    # toggle major grids in current axes (default key 'g')
    # Both here and below (for 'G'), we do nothing if *any* grid (major or
    # minor, x or y) is not in a uniform state, to avoid messing up user
    # customization.
    if (event.key in grid_keys
            # Exclude minor grids not in a uniform state.
            and None not in [_get_uniform_gridstate(ax.xaxis.minorTicks),
                             _get_uniform_gridstate(ax.yaxis.minorTicks)]):
        x_state = _get_uniform_gridstate(ax.xaxis.majorTicks)
        y_state = _get_uniform_gridstate(ax.yaxis.majorTicks)
        cycle = [(False, False), (True, False), (True, True), (False, True)]
        try:
            x_state, y_state = (
                cycle[(cycle.index((x_state, y_state)) + 1) % len(cycle)])
        except ValueError:
            # Exclude major grids not in a uniform state.
            pass
        else:
            # If turning major grids off, also turn minor grids off.
            ax.grid(x_state, which="major" if x_state else "both", axis="x")
            ax.grid(y_state, which="major" if y_state else "both", axis="y")
            canvas.draw_idle()
    # toggle major and minor grids in current axes (default key 'G')
    if (event.key in grid_minor_keys
            # Exclude major grids not in a uniform state.
            and None not in [_get_uniform_gridstate(ax.xaxis.majorTicks),
                             _get_uniform_gridstate(ax.yaxis.majorTicks)]):
        x_state = _get_uniform_gridstate(ax.xaxis.minorTicks)
        y_state = _get_uniform_gridstate(ax.yaxis.minorTicks)
        cycle = [(False, False), (True, False), (True, True), (False, True)]
        try:
            x_state, y_state = (
                cycle[(cycle.index((x_state, y_state)) + 1) % len(cycle)])
        except ValueError:
            # Exclude minor grids not in a uniform state.
            pass
        else:
            ax.grid(x_state, which="both", axis="x")
            ax.grid(y_state, which="both", axis="y")
            canvas.draw_idle()
    # toggle scaling of y-axes between 'log and 'linear' (default key 'l')
    elif event.key in toggle_yscale_keys:
        scale = ax.get_yscale()
        if scale == 'log':
            ax.set_yscale('linear')
            ax.figure.canvas.draw_idle()
        elif scale == 'linear':
            try:
                ax.set_yscale('log')
            except ValueError as exc:
                warnings.warn(str(exc))
                ax.set_yscale('linear')
            ax.figure.canvas.draw_idle()
    # toggle scaling of x-axes between 'log and 'linear' (default key 'k')
    elif event.key in toggle_xscale_keys:
        scalex = ax.get_xscale()
        if scalex == 'log':
            ax.set_xscale('linear')
            ax.figure.canvas.draw_idle()
        elif scalex == 'linear':
            try:
                ax.set_xscale('log')
            except ValueError as exc:
                warnings.warn(str(exc))
                ax.set_xscale('linear')
            ax.figure.canvas.draw_idle()

    elif (event.key.isdigit() and event.key != '0') or event.key in all_keys:
        # keys in list 'all' enables all axes (default key 'a'),
        # otherwise if key is a number only enable this particular axes
        # if it was the axes, where the event was raised
        if not (event.key in all_keys):
            n = int(event.key) - 1
        for i, a in enumerate(canvas.figure.get_axes()):
            # consider axes, in which the event was raised
            # FIXME: Why only this axes?
            if event.x is not None and event.y is not None \
                    and a.in_axes(event):
                if event.key in all_keys:
                    a.set_navigate(True)
                else:
                    a.set_navigate(i == n)


class NonGuiException(Exception):
    pass


class FigureManagerBase(object):
    """
    Helper class for pyplot mode, wraps everything up into a neat bundle

    Attributes
    ----------
    canvas : :class:`FigureCanvasBase`
        The backend-specific canvas instance

    num : int or str
        The figure number

    key_press_handler_id : int
        The default key handler cid, when using the toolmanager.  Can be used
        to disable default key press handling ::

            figure.canvas.mpl_disconnect(
                figure.canvas.manager.key_press_handler_id)
    """
    def __init__(self, canvas, num):
        self.canvas = canvas
        canvas.manager = self  # store a pointer to parent
        self.num = num

        self.key_press_handler_id = None
        if rcParams['toolbar'] != 'toolmanager':
            self.key_press_handler_id = self.canvas.mpl_connect(
                'key_press_event',
                self.key_press)

    def show(self):
        """
        For GUI backends, show the figure window and redraw.
        For non-GUI backends, raise an exception to be caught
        by :meth:`~matplotlib.figure.Figure.show`, for an
        optional warning.
        """
        raise NonGuiException()

    def destroy(self):
        pass

    def full_screen_toggle(self):
        pass

    def resize(self, w, h):
        """"For GUI backends, resize the window (in pixels)."""

    def key_press(self, event):
        """
        Implement the default mpl key bindings defined at
        :ref:`key-event-handling`
        """
        if rcParams['toolbar'] != 'toolmanager':
            key_press_handler(event, self.canvas, self.canvas.toolbar)

    @cbook.deprecated("2.2")
    def show_popup(self, msg):
        """Display message in a popup -- GUI only."""

    def get_window_title(self):
        """Get the title text of the window containing the figure.

        Return None for non-GUI (e.g., PS) backends.
        """
        return 'image'

    def set_window_title(self, title):
        """Set the title text of the window containing the figure.

        This has no effect for non-GUI (e.g., PS) backends.
        """


cursors = tools.cursors


class NavigationToolbar2(object):
    """
    Base class for the navigation cursor, version 2

    backends must implement a canvas that handles connections for
    'button_press_event' and 'button_release_event'.  See
    :meth:`FigureCanvasBase.mpl_connect` for more information


    They must also define

      :meth:`save_figure`
         save the current figure

      :meth:`set_cursor`
         if you want the pointer icon to change

      :meth:`_init_toolbar`
         create your toolbar widget

      :meth:`draw_rubberband` (optional)
         draw the zoom to rect "rubberband" rectangle

      :meth:`press`  (optional)
         whenever a mouse button is pressed, you'll be notified with
         the event

      :meth:`release` (optional)
         whenever a mouse button is released, you'll be notified with
         the event

      :meth:`set_message` (optional)
         display message

      :meth:`set_history_buttons` (optional)
         you can change the history back / forward buttons to
         indicate disabled / enabled state.

    That's it, we'll do the rest!
    """

    # list of toolitems to add to the toolbar, format is:
    # (
    #   text, # the text of the button (often not visible to users)
    #   tooltip_text, # the tooltip shown on hover (where possible)
    #   image_file, # name of the image for the button (without the extension)
    #   name_of_method, # name of the method in NavigationToolbar2 to call
    # )
    toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to  previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
      )

    def __init__(self, canvas):
        self.canvas = canvas
        canvas.toolbar = self
        self._nav_stack = cbook.Stack()
        self._xypress = None  # the location and axis info at the time
                              # of the press
        self._idPress = None
        self._idRelease = None
        self._active = None
        # This cursor will be set after the initial draw.
        self._lastCursor = cursors.POINTER
        self._init_toolbar()
        self._idDrag = self.canvas.mpl_connect(
            'motion_notify_event', self.mouse_move)

        self._ids_zoom = []
        self._zoom_mode = None

        self._button_pressed = None  # determined by the button pressed
                                     # at start

        self.mode = ''  # a mode string for the status bar
        self.set_history_buttons()

    def set_message(self, s):
        """Display a message on toolbar or in status bar."""

    def back(self, *args):
        """move back up the view lim stack"""
        self._nav_stack.back()
        self.set_history_buttons()
        self._update_view()

    @cbook.deprecated("2.1", alternative="canvas.draw_idle")
    def dynamic_update(self):
        self.canvas.draw_idle()

    def draw_rubberband(self, event, x0, y0, x1, y1):
        """Draw a rectangle rubberband to indicate zoom limits.

        Note that it is not guaranteed that ``x0 <= x1`` and ``y0 <= y1``.
        """

    def remove_rubberband(self):
        """Remove the rubberband."""

    def forward(self, *args):
        """Move forward in the view lim stack."""
        self._nav_stack.forward()
        self.set_history_buttons()
        self._update_view()

    def home(self, *args):
        """Restore the original view."""
        self._nav_stack.home()
        self.set_history_buttons()
        self._update_view()

    def _init_toolbar(self):
        """
        This is where you actually build the GUI widgets (called by
        __init__).  The icons ``home.xpm``, ``back.xpm``, ``forward.xpm``,
        ``hand.xpm``, ``zoom_to_rect.xpm`` and ``filesave.xpm`` are standard
        across backends (there are ppm versions in CVS also).

        You just need to set the callbacks

        home         : self.home
        back         : self.back
        forward      : self.forward
        hand         : self.pan
        zoom_to_rect : self.zoom
        filesave     : self.save_figure

        You only need to define the last one - the others are in the base
        class implementation.

        """
        raise NotImplementedError

    def _set_cursor(self, event):
        if not event.inaxes or not self._active:
            if self._lastCursor != cursors.POINTER:
                self.set_cursor(cursors.POINTER)
                self._lastCursor = cursors.POINTER
        else:
            if (self._active == 'ZOOM'
                    and self._lastCursor != cursors.SELECT_REGION):
                self.set_cursor(cursors.SELECT_REGION)
                self._lastCursor = cursors.SELECT_REGION
            elif (self._active == 'PAN' and
                  self._lastCursor != cursors.MOVE):
                self.set_cursor(cursors.MOVE)
                self._lastCursor = cursors.MOVE

    def mouse_move(self, event):
        self._set_cursor(event)

        if event.inaxes and event.inaxes.get_navigate():

            try:
                s = event.inaxes.format_coord(event.xdata, event.ydata)
            except (ValueError, OverflowError):
                pass
            else:
                artists = [a for a in event.inaxes.mouseover_set
                           if a.contains(event) and a.get_visible()]

                if artists:
                    a = cbook._topmost_artist(artists)
                    if a is not event.inaxes.patch:
                        data = a.get_cursor_data(event)
                        if data is not None:
                            s += ' [%s]' % a.format_cursor_data(data)

                if len(self.mode):
                    self.set_message('%s, %s' % (self.mode, s))
                else:
                    self.set_message(s)
        else:
            self.set_message(self.mode)

    def pan(self, *args):
        """Activate the pan/zoom tool. pan with left button, zoom with right"""
        # set the pointer icon and button press funcs to the
        # appropriate callbacks

        if self._active == 'PAN':
            self._active = None
        else:
            self._active = 'PAN'
        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''

        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''

        if self._active:
            self._idPress = self.canvas.mpl_connect(
                'button_press_event', self.press_pan)
            self._idRelease = self.canvas.mpl_connect(
                'button_release_event', self.release_pan)
            self.mode = 'pan/zoom'
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)

        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)

        self.set_message(self.mode)

    def press(self, event):
        """Called whenever a mouse button is pressed."""

    def press_pan(self, event):
        """Callback for mouse button press in pan/zoom mode."""

        if event.button == 1:
            self._button_pressed = 1
        elif event.button == 3:
            self._button_pressed = 3
        else:
            self._button_pressed = None
            return

        if self._nav_stack() is None:
            # set the home button to this view
            self.push_current()

        x, y = event.x, event.y
        self._xypress = []
        for i, a in enumerate(self.canvas.figure.get_axes()):
            if (x is not None and y is not None and a.in_axes(event) and
                    a.get_navigate() and a.can_pan()):
                a.start_pan(x, y, event.button)
                self._xypress.append((a, i))
                self.canvas.mpl_disconnect(self._idDrag)
                self._idDrag = self.canvas.mpl_connect('motion_notify_event',
                                                       self.drag_pan)

        self.press(event)

    def press_zoom(self, event):
        """Callback for mouse button press in zoom to rect mode."""
        # If we're already in the middle of a zoom, pressing another
        # button works to "cancel"
        if self._ids_zoom != []:
            for zoom_id in self._ids_zoom:
                self.canvas.mpl_disconnect(zoom_id)
            self.release(event)
            self.draw()
            self._xypress = None
            self._button_pressed = None
            self._ids_zoom = []
            return

        if event.button == 1:
            self._button_pressed = 1
        elif event.button == 3:
            self._button_pressed = 3
        else:
            self._button_pressed = None
            return

        if self._nav_stack() is None:
            # set the home button to this view
            self.push_current()

        x, y = event.x, event.y
        self._xypress = []
        for i, a in enumerate(self.canvas.figure.get_axes()):
            if (x is not None and y is not None and a.in_axes(event) and
                    a.get_navigate() and a.can_zoom()):
                self._xypress.append((x, y, a, i, a._get_view()))

        id1 = self.canvas.mpl_connect('motion_notify_event', self.drag_zoom)
        id2 = self.canvas.mpl_connect('key_press_event',
                                      self._switch_on_zoom_mode)
        id3 = self.canvas.mpl_connect('key_release_event',
                                      self._switch_off_zoom_mode)

        self._ids_zoom = id1, id2, id3
        self._zoom_mode = event.key

        self.press(event)

    def _switch_on_zoom_mode(self, event):
        self._zoom_mode = event.key
        self.mouse_move(event)

    def _switch_off_zoom_mode(self, event):
        self._zoom_mode = None
        self.mouse_move(event)

    def push_current(self):
        """Push the current view limits and position onto the stack."""
        self._nav_stack.push(
            WeakKeyDictionary(
                {ax: (ax._get_view(),
                      # Store both the original and modified positions.
                      (ax.get_position(True).frozen(),
                       ax.get_position().frozen()))
                 for ax in self.canvas.figure.axes}))
        self.set_history_buttons()

    def release(self, event):
        """Callback for mouse button release."""

    def release_pan(self, event):
        """Callback for mouse button release in pan/zoom mode."""

        if self._button_pressed is None:
            return
        self.canvas.mpl_disconnect(self._idDrag)
        self._idDrag = self.canvas.mpl_connect(
            'motion_notify_event', self.mouse_move)
        for a, ind in self._xypress:
            a.end_pan()
        if not self._xypress:
            return
        self._xypress = []
        self._button_pressed = None
        self.push_current()
        self.release(event)
        self.draw()

    def drag_pan(self, event):
        """Callback for dragging in pan/zoom mode."""
        for a, ind in self._xypress:
            #safer to use the recorded button at the press than current button:
            #multiple button can get pressed during motion...
            a.drag_pan(self._button_pressed, event.key, event.x, event.y)
        self.canvas.draw_idle()

    def drag_zoom(self, event):
        """Callback for dragging in zoom mode."""
        if self._xypress:
            x, y = event.x, event.y
            lastx, lasty, a, ind, view = self._xypress[0]
            (x1, y1), (x2, y2) = np.clip(
                [[lastx, lasty], [x, y]], a.bbox.min, a.bbox.max)
            if self._zoom_mode == "x":
                y1, y2 = a.bbox.intervaly
            elif self._zoom_mode == "y":
                x1, x2 = a.bbox.intervalx
            self.draw_rubberband(event, x1, y1, x2, y2)

    def release_zoom(self, event):
        """Callback for mouse button release in zoom to rect mode."""
        for zoom_id in self._ids_zoom:
            self.canvas.mpl_disconnect(zoom_id)
        self._ids_zoom = []

        self.remove_rubberband()

        if not self._xypress:
            return

        last_a = []

        for cur_xypress in self._xypress:
            x, y = event.x, event.y
            lastx, lasty, a, ind, view = cur_xypress
            # ignore singular clicks - 5 pixels is a threshold
            # allows the user to "cancel" a zoom action
            # by zooming by less than 5 pixels
            if ((abs(x - lastx) < 5 and self._zoom_mode!="y") or
                    (abs(y - lasty) < 5 and self._zoom_mode!="x")):
                self._xypress = None
                self.release(event)
                self.draw()
                return

            # detect twinx,y axes and avoid double zooming
            twinx, twiny = False, False
            if last_a:
                for la in last_a:
                    if a.get_shared_x_axes().joined(a, la):
                        twinx = True
                    if a.get_shared_y_axes().joined(a, la):
                        twiny = True
            last_a.append(a)

            if self._button_pressed == 1:
                direction = 'in'
            elif self._button_pressed == 3:
                direction = 'out'
            else:
                continue

            a._set_view_from_bbox((lastx, lasty, x, y), direction,
                                  self._zoom_mode, twinx, twiny)

        self.draw()
        self._xypress = None
        self._button_pressed = None

        self._zoom_mode = None

        self.push_current()
        self.release(event)

    def draw(self):
        """Redraw the canvases, update the locators."""
        for a in self.canvas.figure.get_axes():
            xaxis = getattr(a, 'xaxis', None)
            yaxis = getattr(a, 'yaxis', None)
            locators = []
            if xaxis is not None:
                locators.append(xaxis.get_major_locator())
                locators.append(xaxis.get_minor_locator())
            if yaxis is not None:
                locators.append(yaxis.get_major_locator())
                locators.append(yaxis.get_minor_locator())

            for loc in locators:
                loc.refresh()
        self.canvas.draw_idle()

    def _update_view(self):
        """Update the viewlim and position from the view and
        position stack for each axes.
        """
        nav_info = self._nav_stack()
        if nav_info is None:
            return
        # Retrieve all items at once to avoid any risk of GC deleting an Axes
        # while in the middle of the loop below.
        items = list(nav_info.items())
        for ax, (view, (pos_orig, pos_active)) in items:
            ax._set_view(view)
            # Restore both the original and modified positions
            ax._set_position(pos_orig, 'original')
            ax._set_position(pos_active, 'active')
        self.canvas.draw_idle()

    def save_figure(self, *args):
        """Save the current figure."""
        raise NotImplementedError

    def set_cursor(self, cursor):
        """Set the current cursor to one of the :class:`Cursors` enums values.

        If required by the backend, this method should trigger an update in
        the backend event loop after the cursor is set, as this method may be
        called e.g. before a long-running task during which the GUI is not
        updated.
        """

    def update(self):
        """Reset the axes stack."""
        self._nav_stack.clear()
        self.set_history_buttons()

    def zoom(self, *args):
        """Activate zoom to rect mode."""
        if self._active == 'ZOOM':
            self._active = None
        else:
            self._active = 'ZOOM'

        if self._idPress is not None:
            self._idPress = self.canvas.mpl_disconnect(self._idPress)
            self.mode = ''

        if self._idRelease is not None:
            self._idRelease = self.canvas.mpl_disconnect(self._idRelease)
            self.mode = ''

        if self._active:
            self._idPress = self.canvas.mpl_connect('button_press_event',
                                                    self.press_zoom)
            self._idRelease = self.canvas.mpl_connect('button_release_event',
                                                      self.release_zoom)
            self.mode = 'zoom rect'
            self.canvas.widgetlock(self)
        else:
            self.canvas.widgetlock.release(self)

        for a in self.canvas.figure.get_axes():
            a.set_navigate_mode(self._active)

        self.set_message(self.mode)

    def set_history_buttons(self):
        """Enable or disable the back/forward button."""


class ToolContainerBase(object):
    """
    Base class for all tool containers, e.g. toolbars.

    Attributes
    ----------
    toolmanager : `ToolManager`
        The tools with which this `ToolContainer` wants to communicate.
    """

    _icon_extension = '.png'
    """
    Toolcontainer button icon image format extension

    **String**: Image extension
    """

    def __init__(self, toolmanager):
        self.toolmanager = toolmanager
        self.toolmanager.toolmanager_connect('tool_removed_event',
                                             self._remove_tool_cbk)

    def _tool_toggled_cbk(self, event):
        """
        Captures the 'tool_trigger_[name]'

        This only gets used for toggled tools
        """
        self.toggle_toolitem(event.tool.name, event.tool.toggled)

    def add_tool(self, tool, group, position=-1):
        """
        Adds a tool to this container

        Parameters
        ----------
        tool : tool_like
            The tool to add, see `ToolManager.get_tool`.
        group : str
            The name of the group to add this tool to.
        position : int (optional)
            The position within the group to place this tool.  Defaults to end.
        """
        tool = self.toolmanager.get_tool(tool)
        image = self._get_image_filename(tool.image)
        toggle = getattr(tool, 'toggled', None) is not None
        self.add_toolitem(tool.name, group, position,
                          image, tool.description, toggle)
        if toggle:
            self.toolmanager.toolmanager_connect('tool_trigger_%s' % tool.name,
                                                 self._tool_toggled_cbk)
            # If initially toggled
            if tool.toggled:
                self.toggle_toolitem(tool.name, True)

    def _remove_tool_cbk(self, event):
        """Captures the 'tool_removed_event' signal and removes the tool."""
        self.remove_toolitem(event.tool.name)

    def _get_image_filename(self, image):
        """Find the image based on its name."""
        if not image:
            return None

        basedir = os.path.join(rcParams['datapath'], 'images')
        possible_images = (
            image,
            image + self._icon_extension,
            os.path.join(basedir, image),
            os.path.join(basedir, image) + self._icon_extension)

        for fname in possible_images:
            if os.path.isfile(fname):
                return fname

    def trigger_tool(self, name):
        """
        Trigger the tool

        Parameters
        ----------
        name : String
            Name (id) of the tool triggered from within the container
        """
        self.toolmanager.trigger_tool(name, sender=self)

    def add_toolitem(self, name, group, position, image, description, toggle):
        """
        Add a toolitem to the container

        This method must get implemented per backend

        The callback associated with the button click event,
        must be **EXACTLY** `self.trigger_tool(name)`

        Parameters
        ----------
        name : string
            Name of the tool to add, this gets used as the tool's ID and as the
            default label of the buttons
        group : String
            Name of the group that this tool belongs to
        position : Int
            Position of the tool within its group, if -1 it goes at the End
        image_file : String
            Filename of the image for the button or `None`
        description : String
            Description of the tool, used for the tooltips
        toggle : Bool
            * `True` : The button is a toggle (change the pressed/unpressed
              state between consecutive clicks)
            * `False` : The button is a normal button (returns to unpressed
              state after release)
        """
        raise NotImplementedError

    def toggle_toolitem(self, name, toggled):
        """
        Toggle the toolitem without firing event

        Parameters
        ----------
        name : String
            Id of the tool to toggle
        toggled : bool
            Whether to set this tool as toggled or not.
        """
        raise NotImplementedError

    def remove_toolitem(self, name):
        """
        Remove a toolitem from the `ToolContainer`

        This method must get implemented per backend

        Called when `ToolManager` emits a `tool_removed_event`

        Parameters
        ----------
        name : string
            Name of the tool to remove
        """
        raise NotImplementedError


class StatusbarBase(object):
    """Base class for the statusbar"""
    def __init__(self, toolmanager):
        self.toolmanager = toolmanager
        self.toolmanager.toolmanager_connect('tool_message_event',
                                             self._message_cbk)

    def _message_cbk(self, event):
        """Captures the 'tool_message_event' and set the message"""
        self.set_message(event.message)

    def set_message(self, s):
        """
        Display a message on toolbar or in status bar

        Parameters
        ----------
        s : str
            Message text
        """
        pass
