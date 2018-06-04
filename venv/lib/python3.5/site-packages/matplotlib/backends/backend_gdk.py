from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import warnings

import gobject
import gtk; gdk = gtk.gdk
import pango
pygtk_version_required = (2,2,0)
if gtk.pygtk_version < pygtk_version_required:
    raise ImportError ("PyGTK %d.%d.%d is installed\n"
                      "PyGTK %d.%d.%d or later is required"
                      % (gtk.pygtk_version + pygtk_version_required))
del pygtk_version_required

import numpy as np

import matplotlib
from matplotlib import rcParams
from matplotlib._pylab_helpers import Gcf
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, GraphicsContextBase,
    RendererBase)
from matplotlib.cbook import warn_deprecated
from matplotlib.mathtext import MathTextParser
from matplotlib.transforms import Affine2D
from matplotlib.backends._backend_gdk import pixbuf_get_pixels_array

backend_version = "%d.%d.%d" % gtk.pygtk_version

# Image formats that this backend supports - for FileChooser and print_figure()
IMAGE_FORMAT = sorted(['bmp', 'eps', 'jpg', 'png', 'ps', 'svg']) # 'raw', 'rgb'
IMAGE_FORMAT_DEFAULT = 'png'


class RendererGDK(RendererBase):
    fontweights = {
        100          : pango.WEIGHT_ULTRALIGHT,
        200          : pango.WEIGHT_LIGHT,
        300          : pango.WEIGHT_LIGHT,
        400          : pango.WEIGHT_NORMAL,
        500          : pango.WEIGHT_NORMAL,
        600          : pango.WEIGHT_BOLD,
        700          : pango.WEIGHT_BOLD,
        800          : pango.WEIGHT_HEAVY,
        900          : pango.WEIGHT_ULTRABOLD,
        'ultralight' : pango.WEIGHT_ULTRALIGHT,
        'light'      : pango.WEIGHT_LIGHT,
        'normal'     : pango.WEIGHT_NORMAL,
        'medium'     : pango.WEIGHT_NORMAL,
        'semibold'   : pango.WEIGHT_BOLD,
        'bold'       : pango.WEIGHT_BOLD,
        'heavy'      : pango.WEIGHT_HEAVY,
        'ultrabold'  : pango.WEIGHT_ULTRABOLD,
        'black'      : pango.WEIGHT_ULTRABOLD,
                   }

    # cache for efficiency, these must be at class, not instance level
    layoutd = {}  # a map from text prop tups to pango layouts
    rotated = {}  # a map from text prop tups to rotated text pixbufs

    def __init__(self, gtkDA, dpi):
        # widget gtkDA is used for:
        #  '<widget>.create_pango_layout(s)'
        #  cmap line below)
        self.gtkDA = gtkDA
        self.dpi   = dpi
        self._cmap = gtkDA.get_colormap()
        self.mathtext_parser = MathTextParser("Agg")

    def set_pixmap (self, pixmap):
        self.gdkDrawable = pixmap

    def set_width_height (self, width, height):
        """w,h is the figure w,h not the pixmap w,h
        """
        self.width, self.height = width, height

    def draw_path(self, gc, path, transform, rgbFace=None):
        transform = transform + Affine2D(). \
            scale(1.0, -1.0).translate(0, self.height)
        polygons = path.to_polygons(transform, self.width, self.height)
        for polygon in polygons:
            # draw_polygon won't take an arbitrary sequence -- it must be a list
            # of tuples
            polygon = [(int(np.round(x)), int(np.round(y))) for x, y in polygon]
            if rgbFace is not None:
                saveColor = gc.gdkGC.foreground
                gc.gdkGC.foreground = gc.rgb_to_gdk_color(rgbFace)
                self.gdkDrawable.draw_polygon(gc.gdkGC, True, polygon)
                gc.gdkGC.foreground = saveColor
            if gc.gdkGC.line_width > 0:
                self.gdkDrawable.draw_lines(gc.gdkGC, polygon)

    def draw_image(self, gc, x, y, im):
        bbox = gc.get_clip_rectangle()

        if bbox != None:
            l,b,w,h = bbox.bounds
            #rectangle = (int(l), self.height-int(b+h),
            #             int(w), int(h))
            # set clip rect?

        rows, cols = im.shape[:2]

        pixbuf = gtk.gdk.Pixbuf(gtk.gdk.COLORSPACE_RGB,
                                has_alpha=True, bits_per_sample=8,
                                width=cols, height=rows)

        array = pixbuf_get_pixels_array(pixbuf)
        array[:, :, :] = im[::-1]

        gc = self.new_gc()


        y = self.height-y-rows

        try: # new in 2.2
            # can use None instead of gc.gdkGC, if don't need clipping
            self.gdkDrawable.draw_pixbuf (gc.gdkGC, pixbuf, 0, 0,
                                          int(x), int(y), cols, rows,
                                          gdk.RGB_DITHER_NONE, 0, 0)
        except AttributeError:
            # deprecated in 2.2
            pixbuf.render_to_drawable(self.gdkDrawable, gc.gdkGC, 0, 0,
                                  int(x), int(y), cols, rows,
                                  gdk.RGB_DITHER_NONE, 0, 0)

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        x, y = int(x), int(y)

        if x < 0 or y < 0: # window has shrunk and text is off the edge
            return

        if angle not in (0,90):
            warnings.warn('backend_gdk: unable to draw text at angles ' +
                          'other than 0 or 90')
        elif ismath:
            self._draw_mathtext(gc, x, y, s, prop, angle)

        elif angle==90:
            self._draw_rotated_text(gc, x, y, s, prop, angle)

        else:
            layout, inkRect, logicalRect = self._get_pango_layout(s, prop)
            l, b, w, h = inkRect
            if (x + w > self.width or y + h > self.height):
                return

            self.gdkDrawable.draw_layout(gc.gdkGC, x, y-h-b, layout)

    def _draw_mathtext(self, gc, x, y, s, prop, angle):
        ox, oy, width, height, descent, font_image, used_characters = \
            self.mathtext_parser.parse(s, self.dpi, prop)

        if angle == 90:
            width, height = height, width
            x -= width
        y -= height

        imw = font_image.get_width()
        imh = font_image.get_height()

        pixbuf = gtk.gdk.Pixbuf(gtk.gdk.COLORSPACE_RGB, has_alpha=True,
                                bits_per_sample=8, width=imw, height=imh)

        array = pixbuf_get_pixels_array(pixbuf)

        rgb = gc.get_rgb()
        array[:,:,0] = int(rgb[0]*255)
        array[:,:,1] = int(rgb[1]*255)
        array[:,:,2] = int(rgb[2]*255)
        array[:,:,3] = (
            np.fromstring(font_image.as_str(), np.uint8).reshape((imh, imw)))

        # can use None instead of gc.gdkGC, if don't need clipping
        self.gdkDrawable.draw_pixbuf(gc.gdkGC, pixbuf, 0, 0,
                                     int(x), int(y), imw, imh,
                                     gdk.RGB_DITHER_NONE, 0, 0)

    def _draw_rotated_text(self, gc, x, y, s, prop, angle):
        """
        Draw the text rotated 90 degrees, other angles are not supported
        """
        # this function (and its called functions) is a bottleneck
        # Pango 1.6 supports rotated text, but pygtk 2.4.0 does not yet have
        # wrapper functions
        # GTK+ 2.6 pixbufs support rotation

        gdrawable = self.gdkDrawable
        ggc = gc.gdkGC

        layout, inkRect, logicalRect = self._get_pango_layout(s, prop)
        l, b, w, h = inkRect
        x = int(x-h)
        y = int(y-w)

        if (x < 0 or y < 0 or # window has shrunk and text is off the edge
            x + w > self.width or y + h > self.height):
            return

        key = (x,y,s,angle,hash(prop))
        imageVert = self.rotated.get(key)
        if imageVert != None:
            gdrawable.draw_image(ggc, imageVert, 0, 0, x, y, h, w)
            return

        imageBack = gdrawable.get_image(x, y, w, h)
        imageVert = gdrawable.get_image(x, y, h, w)
        imageFlip = gtk.gdk.Image(type=gdk.IMAGE_FASTEST,
                                  visual=gdrawable.get_visual(),
                                  width=w, height=h)
        if imageFlip == None or imageBack == None or imageVert == None:
            warnings.warn("Could not renderer vertical text")
            return
        imageFlip.set_colormap(self._cmap)
        for i in range(w):
            for j in range(h):
                imageFlip.put_pixel(i, j, imageVert.get_pixel(j,w-i-1) )

        gdrawable.draw_image(ggc, imageFlip, 0, 0, x, y, w, h)
        gdrawable.draw_layout(ggc, x, y-b, layout)

        imageIn  = gdrawable.get_image(x, y, w, h)
        for i in range(w):
            for j in range(h):
                imageVert.put_pixel(j, i, imageIn.get_pixel(w-i-1,j) )

        gdrawable.draw_image(ggc, imageBack, 0, 0, x, y, w, h)
        gdrawable.draw_image(ggc, imageVert, 0, 0, x, y, h, w)
        self.rotated[key] = imageVert

    def _get_pango_layout(self, s, prop):
        """
        Create a pango layout instance for Text 's' with properties 'prop'.
        Return - pango layout (from cache if already exists)

        Note that pango assumes a logical DPI of 96
        Ref: pango/fonts.c/pango_font_description_set_size() manual page
        """
        # problem? - cache gets bigger and bigger, is never cleared out
        # two (not one) layouts are created for every text item s (then they
        # are cached) - why?

        key = self.dpi, s, hash(prop)
        value = self.layoutd.get(key)
        if value != None:
            return value

        size = prop.get_size_in_points() * self.dpi / 96.0
        size = np.round(size)

        font_str = '%s, %s %i' % (prop.get_name(), prop.get_style(), size,)
        font = pango.FontDescription(font_str)

        # later - add fontweight to font_str
        font.set_weight(self.fontweights[prop.get_weight()])

        layout = self.gtkDA.create_pango_layout(s)
        layout.set_font_description(font)
        inkRect, logicalRect = layout.get_pixel_extents()

        self.layoutd[key] = layout, inkRect, logicalRect
        return layout, inkRect, logicalRect

    def flipy(self):
        return True

    def get_canvas_width_height(self):
        return self.width, self.height

    def get_text_width_height_descent(self, s, prop, ismath):
        if ismath:
            ox, oy, width, height, descent, font_image, used_characters = \
                self.mathtext_parser.parse(s, self.dpi, prop)
            return width, height, descent

        layout, inkRect, logicalRect = self._get_pango_layout(s, prop)
        l, b, w, h = inkRect
        ll, lb, lw, lh = logicalRect

        return w, h + 1, h - lh

    def new_gc(self):
        return GraphicsContextGDK(renderer=self)

    def points_to_pixels(self, points):
        return points/72.0 * self.dpi


class GraphicsContextGDK(GraphicsContextBase):
    # a cache shared by all class instances
    _cached = {}  # map: rgb color -> gdk.Color

    _joind = {
        'bevel' : gdk.JOIN_BEVEL,
        'miter' : gdk.JOIN_MITER,
        'round' : gdk.JOIN_ROUND,
        }

    _capd = {
        'butt'       : gdk.CAP_BUTT,
        'projecting' : gdk.CAP_PROJECTING,
        'round'      : gdk.CAP_ROUND,
        }


    def __init__(self, renderer):
        GraphicsContextBase.__init__(self)
        self.renderer = renderer
        self.gdkGC    = gtk.gdk.GC(renderer.gdkDrawable)
        self._cmap    = renderer._cmap


    def rgb_to_gdk_color(self, rgb):
        """
        rgb - an RGB tuple (three 0.0-1.0 values)
        return an allocated gtk.gdk.Color
        """
        try:
            return self._cached[tuple(rgb)]
        except KeyError:
            color = self._cached[tuple(rgb)] = \
                    self._cmap.alloc_color(
                        int(rgb[0]*65535),int(rgb[1]*65535),int(rgb[2]*65535))
            return color


    #def set_antialiased(self, b):
        # anti-aliasing is not supported by GDK

    def set_capstyle(self, cs):
        GraphicsContextBase.set_capstyle(self, cs)
        self.gdkGC.cap_style = self._capd[self._capstyle]


    def set_clip_rectangle(self, rectangle):
        GraphicsContextBase.set_clip_rectangle(self, rectangle)
        if rectangle is None:
            return
        l,b,w,h = rectangle.bounds
        rectangle = (int(l), self.renderer.height-int(b+h)+1,
                     int(w), int(h))
        #rectangle = (int(l), self.renderer.height-int(b+h),
        #             int(w+1), int(h+2))
        self.gdkGC.set_clip_rectangle(rectangle)

    def set_dashes(self, dash_offset, dash_list):
        GraphicsContextBase.set_dashes(self, dash_offset, dash_list)

        if dash_list == None:
            self.gdkGC.line_style = gdk.LINE_SOLID
        else:
            pixels = self.renderer.points_to_pixels(np.asarray(dash_list))
            dl = [max(1, int(np.round(val))) for val in pixels]
            self.gdkGC.set_dashes(dash_offset, dl)
            self.gdkGC.line_style = gdk.LINE_ON_OFF_DASH


    def set_foreground(self, fg, isRGBA=False):
        GraphicsContextBase.set_foreground(self, fg, isRGBA)
        self.gdkGC.foreground = self.rgb_to_gdk_color(self.get_rgb())


    def set_joinstyle(self, js):
        GraphicsContextBase.set_joinstyle(self, js)
        self.gdkGC.join_style = self._joind[self._joinstyle]


    def set_linewidth(self, w):
        GraphicsContextBase.set_linewidth(self, w)
        if w == 0:
            self.gdkGC.line_width = 0
        else:
            pixels = self.renderer.points_to_pixels(w)
            self.gdkGC.line_width = max(1, int(np.round(pixels)))


class FigureCanvasGDK (FigureCanvasBase):
    def __init__(self, figure):
        FigureCanvasBase.__init__(self, figure)
        if self.__class__ == matplotlib.backends.backend_gdk.FigureCanvasGDK:
            warn_deprecated('2.0', message="The GDK backend is "
                            "deprecated. It is untested, known to be "
                            "broken and will be removed in Matplotlib 3.0. "
                            "Use the Agg backend instead. "
                            "See Matplotlib usage FAQ for"
                            " more info on backends.",
                            alternative="Agg")
        self._renderer_init()

    def _renderer_init(self):
        self._renderer = RendererGDK (gtk.DrawingArea(), self.figure.dpi)

    def _render_figure(self, pixmap, width, height):
        self._renderer.set_pixmap (pixmap)
        self._renderer.set_width_height (width, height)
        self.figure.draw (self._renderer)

    filetypes = FigureCanvasBase.filetypes.copy()
    filetypes['jpg'] = 'JPEG'
    filetypes['jpeg'] = 'JPEG'

    def print_jpeg(self, filename, *args, **kwargs):
        return self._print_image(filename, 'jpeg')
    print_jpg = print_jpeg

    def print_png(self, filename, *args, **kwargs):
        return self._print_image(filename, 'png')

    def _print_image(self, filename, format, *args, **kwargs):
        width, height = self.get_width_height()
        pixmap = gtk.gdk.Pixmap (None, width, height, depth=24)
        self._render_figure(pixmap, width, height)

        # jpg colors don't match the display very well, png colors match
        # better
        pixbuf = gtk.gdk.Pixbuf(gtk.gdk.COLORSPACE_RGB, 0, 8,
                                width, height)
        pixbuf.get_from_drawable(pixmap, pixmap.get_colormap(),
                                 0, 0, 0, 0, width, height)

        # set the default quality, if we are writing a JPEG.
        # http://www.pygtk.org/docs/pygtk/class-gdkpixbuf.html#method-gdkpixbuf--save
        options = {k: kwargs[k] for k in ['quality'] if k in kwargs}
        if format in ['jpg', 'jpeg']:
            options.setdefault('quality', rcParams['savefig.jpeg_quality'])
            options['quality'] = str(options['quality'])

        pixbuf.save(filename, format, options=options)


@_Backend.export
class _BackendGDK(_Backend):
    FigureCanvas = FigureCanvasGDK
    FigureManager = FigureManagerBase
