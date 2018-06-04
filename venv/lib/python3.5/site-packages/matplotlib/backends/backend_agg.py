"""
An agg http://antigrain.com/ backend

Features that are implemented

 * capstyles and join styles
 * dashes
 * linewidth
 * lines, rectangles, ellipses
 * clipping to a rectangle
 * output to RGBA and PNG, optionally JPEG and TIFF
 * alpha blending
 * DPI scaling properly - everything scales properly (dashes, linewidths, etc)
 * draw polygon
 * freetype2 w/ ft2font

TODO:

  * integrate screen dpi w/ ppi and text

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import threading
import numpy as np
from collections import OrderedDict
from math import radians, cos, sin
from matplotlib import cbook, rcParams, __version__
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, RendererBase, cursors)
from matplotlib.cbook import maxdict
from matplotlib.figure import Figure
from matplotlib.font_manager import findfont, get_font
from matplotlib.ft2font import (LOAD_FORCE_AUTOHINT, LOAD_NO_HINTING,
                                LOAD_DEFAULT, LOAD_NO_AUTOHINT)
from matplotlib.mathtext import MathTextParser
from matplotlib.path import Path
from matplotlib.transforms import Bbox, BboxBase
from matplotlib import colors as mcolors

from matplotlib.backends._backend_agg import RendererAgg as _RendererAgg
from matplotlib import _png

try:
    from PIL import Image
    _has_pil = True
except ImportError:
    _has_pil = False

backend_version = 'v2.2'

def get_hinting_flag():
    mapping = {
        True: LOAD_FORCE_AUTOHINT,
        False: LOAD_NO_HINTING,
        'either': LOAD_DEFAULT,
        'native': LOAD_NO_AUTOHINT,
        'auto': LOAD_FORCE_AUTOHINT,
        'none': LOAD_NO_HINTING
        }
    return mapping[rcParams['text.hinting']]


class RendererAgg(RendererBase):
    """
    The renderer handles all the drawing primitives using a graphics
    context instance that controls the colors/styles
    """

    @property
    @cbook.deprecated("2.2")
    def debug(self):
        return 1

    # we want to cache the fonts at the class level so that when
    # multiple figures are created we can reuse them.  This helps with
    # a bug on windows where the creation of too many figures leads to
    # too many open file handles.  However, storing them at the class
    # level is not thread safe.  The solution here is to let the
    # FigureCanvas acquire a lock on the fontd at the start of the
    # draw, and release it when it is done.  This allows multiple
    # renderers to share the cached fonts, but only one figure can
    # draw at time and so the font cache is used by only one
    # renderer at a time.

    lock = threading.RLock()

    def __init__(self, width, height, dpi):
        RendererBase.__init__(self)

        self.dpi = dpi
        self.width = width
        self.height = height
        self._renderer = _RendererAgg(int(width), int(height), dpi)
        self._filter_renderers = []

        self._update_methods()
        self.mathtext_parser = MathTextParser('Agg')

        self.bbox = Bbox.from_bounds(0, 0, self.width, self.height)

    def __getstate__(self):
        # We only want to preserve the init keywords of the Renderer.
        # Anything else can be re-created.
        return {'width': self.width, 'height': self.height, 'dpi': self.dpi}

    def __setstate__(self, state):
        self.__init__(state['width'], state['height'], state['dpi'])

    def _get_hinting_flag(self):
        if rcParams['text.hinting']:
            return LOAD_FORCE_AUTOHINT
        else:
            return LOAD_NO_HINTING

    # for filtering to work with rasterization, methods needs to be wrapped.
    # maybe there is better way to do it.
    def draw_markers(self, *kl, **kw):
        return self._renderer.draw_markers(*kl, **kw)

    def draw_path_collection(self, *kl, **kw):
        return self._renderer.draw_path_collection(*kl, **kw)

    def _update_methods(self):
        self.draw_quad_mesh = self._renderer.draw_quad_mesh
        self.draw_gouraud_triangle = self._renderer.draw_gouraud_triangle
        self.draw_gouraud_triangles = self._renderer.draw_gouraud_triangles
        self.draw_image = self._renderer.draw_image
        self.copy_from_bbox = self._renderer.copy_from_bbox
        self.get_content_extents = self._renderer.get_content_extents

    def tostring_rgba_minimized(self):
        extents = self.get_content_extents()
        bbox = [[extents[0], self.height - (extents[1] + extents[3])],
                [extents[0] + extents[2], self.height - extents[1]]]
        region = self.copy_from_bbox(bbox)
        return np.array(region), extents

    def draw_path(self, gc, path, transform, rgbFace=None):
        """
        Draw the path
        """
        nmax = rcParams['agg.path.chunksize'] # here at least for testing
        npts = path.vertices.shape[0]

        if (nmax > 100 and npts > nmax and path.should_simplify and
                rgbFace is None and gc.get_hatch() is None):
            nch = np.ceil(npts / nmax)
            chsize = int(np.ceil(npts / nch))
            i0 = np.arange(0, npts, chsize)
            i1 = np.zeros_like(i0)
            i1[:-1] = i0[1:] - 1
            i1[-1] = npts
            for ii0, ii1 in zip(i0, i1):
                v = path.vertices[ii0:ii1, :]
                c = path.codes
                if c is not None:
                    c = c[ii0:ii1]
                    c[0] = Path.MOVETO  # move to end of last chunk
                p = Path(v, c)
                try:
                    self._renderer.draw_path(gc, p, transform, rgbFace)
                except OverflowError:
                    raise OverflowError("Exceeded cell block limit (set "
                                        "'agg.path.chunksize' rcparam)")
        else:
            try:
                self._renderer.draw_path(gc, path, transform, rgbFace)
            except OverflowError:
                raise OverflowError("Exceeded cell block limit (set "
                                    "'agg.path.chunksize' rcparam)")


    def draw_mathtext(self, gc, x, y, s, prop, angle):
        """
        Draw the math text using matplotlib.mathtext
        """
        ox, oy, width, height, descent, font_image, used_characters = \
            self.mathtext_parser.parse(s, self.dpi, prop)

        xd = descent * sin(radians(angle))
        yd = descent * cos(radians(angle))
        x = np.round(x + ox + xd)
        y = np.round(y - oy + yd)
        self._renderer.draw_text_image(font_image, x, y + 1, angle, gc)

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        """
        Render the text
        """
        if ismath:
            return self.draw_mathtext(gc, x, y, s, prop, angle)

        flags = get_hinting_flag()
        font = self._get_agg_font(prop)

        if font is None:
            return None
        if len(s) == 1 and ord(s) > 127:
            font.load_char(ord(s), flags=flags)
        else:
            # We pass '0' for angle here, since it will be rotated (in raster
            # space) in the following call to draw_text_image).
            font.set_text(s, 0, flags=flags)
        font.draw_glyphs_to_bitmap(antialiased=rcParams['text.antialiased'])
        d = font.get_descent() / 64.0
        # The descent needs to be adjusted for the angle.
        xo, yo = font.get_bitmap_offset()
        xo /= 64.0
        yo /= 64.0
        xd = -d * sin(radians(angle))
        yd = d * cos(radians(angle))

        self._renderer.draw_text_image(
            font, np.round(x - xd + xo), np.round(y + yd + yo) + 1, angle, gc)

    def get_text_width_height_descent(self, s, prop, ismath):
        """
        Get the width, height, and descent (offset from the bottom
        to the baseline), in display coords, of the string *s* with
        :class:`~matplotlib.font_manager.FontProperties` *prop*
        """
        if ismath in ["TeX", "TeX!"]:
            # todo: handle props
            size = prop.get_size_in_points()
            texmanager = self.get_texmanager()
            fontsize = prop.get_size_in_points()
            w, h, d = texmanager.get_text_width_height_descent(
                s, fontsize, renderer=self)
            return w, h, d

        if ismath:
            ox, oy, width, height, descent, fonts, used_characters = \
                self.mathtext_parser.parse(s, self.dpi, prop)
            return width, height, descent

        flags = get_hinting_flag()
        font = self._get_agg_font(prop)
        font.set_text(s, 0.0, flags=flags)
        w, h = font.get_width_height()  # width and height of unrotated string
        d = font.get_descent()
        w /= 64.0  # convert from subpixels
        h /= 64.0
        d /= 64.0
        return w, h, d

    def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
        # todo, handle props, angle, origins
        size = prop.get_size_in_points()

        texmanager = self.get_texmanager()

        Z = texmanager.get_grey(s, size, self.dpi)
        Z = np.array(Z * 255.0, np.uint8)

        w, h, d = self.get_text_width_height_descent(s, prop, ismath)
        xd = d * sin(radians(angle))
        yd = d * cos(radians(angle))
        x = np.round(x + xd)
        y = np.round(y + yd)

        self._renderer.draw_text_image(Z, x, y, angle, gc)

    def get_canvas_width_height(self):
        'return the canvas width and height in display coords'
        return self.width, self.height

    def _get_agg_font(self, prop):
        """
        Get the font for text instance t, cacheing for efficiency
        """
        fname = findfont(prop)
        font = get_font(fname)

        font.clear()
        size = prop.get_size_in_points()
        font.set_size(size, self.dpi)

        return font

    def points_to_pixels(self, points):
        """
        convert point measures to pixes using dpi and the pixels per
        inch of the display
        """
        return points*self.dpi/72.0

    def tostring_rgb(self):
        return self._renderer.tostring_rgb()

    def tostring_argb(self):
        return self._renderer.tostring_argb()

    def buffer_rgba(self):
        return self._renderer.buffer_rgba()

    def clear(self):
        self._renderer.clear()

    def option_image_nocomposite(self):
        # It is generally faster to composite each image directly to
        # the Figure, and there's no file size benefit to compositing
        # with the Agg backend
        return True

    def option_scale_image(self):
        """
        agg backend doesn't support arbitrary scaling of image.
        """
        return False

    def restore_region(self, region, bbox=None, xy=None):
        """
        Restore the saved region. If bbox (instance of BboxBase, or
        its extents) is given, only the region specified by the bbox
        will be restored. *xy* (a tuple of two floasts) optionally
        specifies the new position (the LLC of the original region,
        not the LLC of the bbox) where the region will be restored.

        >>> region = renderer.copy_from_bbox()
        >>> x1, y1, x2, y2 = region.get_extents()
        >>> renderer.restore_region(region, bbox=(x1+dx, y1, x2, y2),
        ...                         xy=(x1-dx, y1))

        """
        if bbox is not None or xy is not None:
            if bbox is None:
                x1, y1, x2, y2 = region.get_extents()
            elif isinstance(bbox, BboxBase):
                x1, y1, x2, y2 = bbox.extents
            else:
                x1, y1, x2, y2 = bbox

            if xy is None:
                ox, oy = x1, y1
            else:
                ox, oy = xy

            # The incoming data is float, but the _renderer type-checking wants
            # to see integers.
            self._renderer.restore_region(region, int(x1), int(y1),
                                          int(x2), int(y2), int(ox), int(oy))

        else:
            self._renderer.restore_region(region)

    def start_filter(self):
        """
        Start filtering. It simply create a new canvas (the old one is saved).
        """
        self._filter_renderers.append(self._renderer)
        self._renderer = _RendererAgg(int(self.width), int(self.height),
                                      self.dpi)
        self._update_methods()

    def stop_filter(self, post_processing):
        """
        Save the plot in the current canvas as a image and apply
        the *post_processing* function.

           def post_processing(image, dpi):
             # ny, nx, depth = image.shape
             # image (numpy array) has RGBA channels and has a depth of 4.
             ...
             # create a new_image (numpy array of 4 channels, size can be
             # different). The resulting image may have offsets from
             # lower-left corner of the original image
             return new_image, offset_x, offset_y

        The saved renderer is restored and the returned image from
        post_processing is plotted (using draw_image) on it.
        """

        # WARNING:  For agg_filter to work, the renderer's method need to
        # overridden in the class. See draw_markers and draw_path_collections.

        width, height = int(self.width), int(self.height)

        buffer, bounds = self.tostring_rgba_minimized()

        l, b, w, h = bounds

        self._renderer = self._filter_renderers.pop()
        self._update_methods()

        if w > 0 and h > 0:
            img = np.fromstring(buffer, np.uint8)
            img, ox, oy = post_processing(img.reshape((h, w, 4)) / 255.,
                                          self.dpi)
            gc = self.new_gc()
            if img.dtype.kind == 'f':
                img = np.asarray(img * 255., np.uint8)
            img = img[::-1]
            self._renderer.draw_image(
                gc, l + ox, height - b - h + oy, img)


class FigureCanvasAgg(FigureCanvasBase):
    """
    The canvas the figure renders into.  Calls the draw and print fig
    methods, creates the renderers, etc...

    Attributes
    ----------
    figure : `matplotlib.figure.Figure`
        A high-level Figure instance

    """

    def copy_from_bbox(self, bbox):
        renderer = self.get_renderer()
        return renderer.copy_from_bbox(bbox)

    def restore_region(self, region, bbox=None, xy=None):
        renderer = self.get_renderer()
        return renderer.restore_region(region, bbox, xy)

    def draw(self):
        """
        Draw the figure using the renderer
        """
        self.renderer = self.get_renderer(cleared=True)
        # acquire a lock on the shared font cache
        RendererAgg.lock.acquire()

        toolbar = self.toolbar
        try:
            # if toolbar:
            #     toolbar.set_cursor(cursors.WAIT)
            self.figure.draw(self.renderer)
            # A GUI class may be need to update a window using this draw, so
            # don't forget to call the superclass.
            super(FigureCanvasAgg, self).draw()
        finally:
            # if toolbar:
            #     toolbar.set_cursor(toolbar._lastCursor)
            RendererAgg.lock.release()

    def get_renderer(self, cleared=False):
        l, b, w, h = self.figure.bbox.bounds
        key = w, h, self.figure.dpi
        try: self._lastKey, self.renderer
        except AttributeError: need_new_renderer = True
        else:  need_new_renderer = (self._lastKey != key)

        if need_new_renderer:
            self.renderer = RendererAgg(w, h, self.figure.dpi)
            self._lastKey = key
        elif cleared:
            self.renderer.clear()
        return self.renderer

    def tostring_rgb(self):
        '''Get the image as an RGB byte string

        `draw` must be called at least once before this function will work and
        to update the renderer for any subsequent changes to the Figure.

        Returns
        -------
        bytes
        '''
        return self.renderer.tostring_rgb()

    def tostring_argb(self):
        '''Get the image as an ARGB byte string

        `draw` must be called at least once before this function will work and
        to update the renderer for any subsequent changes to the Figure.

        Returns
        -------
        bytes

        '''
        return self.renderer.tostring_argb()

    def buffer_rgba(self):
        '''Get the image as an RGBA byte string

        `draw` must be called at least once before this function will work and
        to update the renderer for any subsequent changes to the Figure.

        Returns
        -------
        bytes
        '''
        return self.renderer.buffer_rgba()

    def print_raw(self, filename_or_obj, *args, **kwargs):
        FigureCanvasAgg.draw(self)
        renderer = self.get_renderer()
        original_dpi = renderer.dpi
        renderer.dpi = self.figure.dpi
        if isinstance(filename_or_obj, six.string_types):
            fileobj = open(filename_or_obj, 'wb')
            close = True
        else:
            fileobj = filename_or_obj
            close = False
        try:
            fileobj.write(renderer._renderer.buffer_rgba())
        finally:
            if close:
                fileobj.close()
            renderer.dpi = original_dpi
    print_rgba = print_raw

    def print_png(self, filename_or_obj, *args, **kwargs):
        FigureCanvasAgg.draw(self)
        renderer = self.get_renderer()
        original_dpi = renderer.dpi
        renderer.dpi = self.figure.dpi

        version_str = 'matplotlib version ' + __version__ + \
            ', http://matplotlib.org/'
        metadata = OrderedDict({'Software': version_str})
        user_metadata = kwargs.pop("metadata", None)
        if user_metadata is not None:
            metadata.update(user_metadata)

        try:
            with cbook.open_file_cm(filename_or_obj, "wb") as fh:
                _png.write_png(renderer._renderer, fh,
                               self.figure.dpi, metadata=metadata)
        finally:
            renderer.dpi = original_dpi

    def print_to_buffer(self):
        FigureCanvasAgg.draw(self)
        renderer = self.get_renderer()
        original_dpi = renderer.dpi
        renderer.dpi = self.figure.dpi
        try:
            result = (renderer._renderer.buffer_rgba(),
                      (int(renderer.width), int(renderer.height)))
        finally:
            renderer.dpi = original_dpi
        return result

    if _has_pil:
        # add JPEG support
        def print_jpg(self, filename_or_obj, *args, **kwargs):
            """
            Other Parameters
            ----------------
            quality : int
                The image quality, on a scale from 1 (worst) to
                95 (best). The default is 95, if not given in the
                matplotlibrc file in the savefig.jpeg_quality parameter.
                Values above 95 should be avoided; 100 completely
                disables the JPEG quantization stage.

            optimize : bool
                If present, indicates that the encoder should
                make an extra pass over the image in order to select
                optimal encoder settings.

            progressive : bool
                If present, indicates that this image
                should be stored as a progressive JPEG file.
            """
            buf, size = self.print_to_buffer()
            if kwargs.pop("dryrun", False):
                return
            # The image is "pasted" onto a white background image to safely
            # handle any transparency
            image = Image.frombuffer('RGBA', size, buf, 'raw', 'RGBA', 0, 1)
            rgba = mcolors.to_rgba(rcParams['savefig.facecolor'])
            color = tuple([int(x * 255.0) for x in rgba[:3]])
            background = Image.new('RGB', size, color)
            background.paste(image, image)
            options = {k: kwargs[k]
                       for k in ['quality', 'optimize', 'progressive', 'dpi']
                       if k in kwargs}
            options.setdefault('quality', rcParams['savefig.jpeg_quality'])
            if 'dpi' in options:
                # Set the same dpi in both x and y directions
                options['dpi'] = (options['dpi'], options['dpi'])

            return background.save(filename_or_obj, format='jpeg', **options)
        print_jpeg = print_jpg

        # add TIFF support
        def print_tif(self, filename_or_obj, *args, **kwargs):
            buf, size = self.print_to_buffer()
            if kwargs.pop("dryrun", False):
                return
            image = Image.frombuffer('RGBA', size, buf, 'raw', 'RGBA', 0, 1)
            dpi = (self.figure.dpi, self.figure.dpi)
            return image.save(filename_or_obj, format='tiff',
                              dpi=dpi)
        print_tiff = print_tif


@_Backend.export
class _BackendAgg(_Backend):
    FigureCanvas = FigureCanvasAgg
    FigureManager = FigureManagerBase
