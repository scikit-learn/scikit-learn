from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict

import six
from six import unichr
from six.moves import xrange

import base64
import codecs
import gzip
import hashlib
import io
import logging
import re
import uuid

import numpy as np

from matplotlib import cbook, __version__, rcParams
from matplotlib.backend_bases import (
     _Backend, FigureCanvasBase, FigureManagerBase, RendererBase)
from matplotlib.backends.backend_mixed import MixedModeRenderer
from matplotlib.cbook import is_writable_file_like, maxdict
from matplotlib.colors import rgb2hex
from matplotlib.font_manager import findfont, get_font
from matplotlib.ft2font import LOAD_NO_HINTING
from matplotlib.mathtext import MathTextParser
from matplotlib.path import Path
from matplotlib import _path
from matplotlib.transforms import Affine2D, Affine2DBase
from matplotlib import _png

_log = logging.getLogger(__name__)

backend_version = __version__

# ----------------------------------------------------------------------
# SimpleXMLWriter class
#
# Based on an original by Fredrik Lundh, but modified here to:
#   1. Support modern Python idioms
#   2. Remove encoding support (it's handled by the file writer instead)
#   3. Support proper indentation
#   4. Minify things a little bit

# --------------------------------------------------------------------
# The SimpleXMLWriter module is
#
# Copyright (c) 2001-2004 by Fredrik Lundh
#
# By obtaining, using, and/or copying this software and/or its
# associated documentation, you agree that you have read, understood,
# and will comply with the following terms and conditions:
#
# Permission to use, copy, modify, and distribute this software and
# its associated documentation for any purpose and without fee is
# hereby granted, provided that the above copyright notice appears in
# all copies, and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of
# Secret Labs AB or the author not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD
# TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANT-
# ABILITY AND FITNESS.  IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR
# BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY
# DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
# WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
# ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
# OF THIS SOFTWARE.
# --------------------------------------------------------------------

def escape_cdata(s):
    s = s.replace("&", "&amp;")
    s = s.replace("<", "&lt;")
    s = s.replace(">", "&gt;")
    return s

_escape_xml_comment = re.compile(r'-(?=-)')
def escape_comment(s):
    s = escape_cdata(s)
    return _escape_xml_comment.sub('- ', s)

def escape_attrib(s):
    s = s.replace("&", "&amp;")
    s = s.replace("'", "&apos;")
    s = s.replace("\"", "&quot;")
    s = s.replace("<", "&lt;")
    s = s.replace(">", "&gt;")
    return s

def short_float_fmt(x):
    """
    Create a short string representation of a float, which is %f
    formatting with trailing zeros and the decimal point removed.
    """
    return '{0:f}'.format(x).rstrip('0').rstrip('.')

##
# XML writer class.
#
# @param file A file or file-like object.  This object must implement
#    a <b>write</b> method that takes an 8-bit string.

class XMLWriter(object):
    def __init__(self, file):
        self.__write = file.write
        if hasattr(file, "flush"):
            self.flush = file.flush
        self.__open = 0 # true if start tag is open
        self.__tags = []
        self.__data = []
        self.__indentation = " " * 64

    def __flush(self, indent=True):
        # flush internal buffers
        if self.__open:
            if indent:
                self.__write(">\n")
            else:
                self.__write(">")
            self.__open = 0
        if self.__data:
            data = ''.join(self.__data)
            self.__write(escape_cdata(data))
            self.__data = []

    ## Opens a new element.  Attributes can be given as keyword
    # arguments, or as a string/string dictionary. The method returns
    # an opaque identifier that can be passed to the <b>close</b>
    # method, to close all open elements up to and including this one.
    #
    # @param tag Element tag.
    # @param attrib Attribute dictionary.  Alternatively, attributes
    #    can be given as keyword arguments.
    # @return An element identifier.

    def start(self, tag, attrib={}, **extra):
        self.__flush()
        tag = escape_cdata(tag)
        self.__data = []
        self.__tags.append(tag)
        self.__write(self.__indentation[:len(self.__tags) - 1])
        self.__write("<%s" % tag)
        if attrib or extra:
            attrib = attrib.copy()
            attrib.update(extra)
            attrib = sorted(six.iteritems(attrib))
            for k, v in attrib:
                if not v == '':
                    k = escape_cdata(k)
                    v = escape_attrib(v)
                    self.__write(" %s=\"%s\"" % (k, v))
        self.__open = 1
        return len(self.__tags)-1

    ##
    # Adds a comment to the output stream.
    #
    # @param comment Comment text, as a Unicode string.

    def comment(self, comment):
        self.__flush()
        self.__write(self.__indentation[:len(self.__tags)])
        self.__write("<!-- %s -->\n" % escape_comment(comment))

    ##
    # Adds character data to the output stream.
    #
    # @param text Character data, as a Unicode string.

    def data(self, text):
        self.__data.append(text)

    ##
    # Closes the current element (opened by the most recent call to
    # <b>start</b>).
    #
    # @param tag Element tag.  If given, the tag must match the start
    #    tag.  If omitted, the current element is closed.

    def end(self, tag=None, indent=True):
        if tag:
            assert self.__tags, "unbalanced end(%s)" % tag
            assert escape_cdata(tag) == self.__tags[-1],\
                   "expected end(%s), got %s" % (self.__tags[-1], tag)
        else:
            assert self.__tags, "unbalanced end()"
        tag = self.__tags.pop()
        if self.__data:
            self.__flush(indent)
        elif self.__open:
            self.__open = 0
            self.__write("/>\n")
            return
        if indent:
            self.__write(self.__indentation[:len(self.__tags)])
        self.__write("</%s>\n" % tag)

    ##
    # Closes open elements, up to (and including) the element identified
    # by the given identifier.
    #
    # @param id Element identifier, as returned by the <b>start</b> method.

    def close(self, id):
        while len(self.__tags) > id:
            self.end()

    ##
    # Adds an entire element.  This is the same as calling <b>start</b>,
    # <b>data</b>, and <b>end</b> in sequence. The <b>text</b> argument
    # can be omitted.

    def element(self, tag, text=None, attrib={}, **extra):
        self.start(*(tag, attrib), **extra)
        if text:
            self.data(text)
        self.end(indent=False)

    ##
    # Flushes the output stream.

    def flush(self):
        pass # replaced by the constructor

# ----------------------------------------------------------------------

def generate_transform(transform_list=[]):
    if len(transform_list):
        output = io.StringIO()
        for type, value in transform_list:
            if type == 'scale' and (value == (1.0,) or value == (1.0, 1.0)):
                continue
            if type == 'translate' and value == (0.0, 0.0):
                continue
            if type == 'rotate' and value == (0.0,):
                continue
            if type == 'matrix' and isinstance(value, Affine2DBase):
                value = value.to_values()

            output.write('%s(%s)' % (
                type, ' '.join(short_float_fmt(x) for x in value)))
        return output.getvalue()
    return ''

def generate_css(attrib={}):
    if attrib:
        output = io.StringIO()
        attrib = sorted(six.iteritems(attrib))
        for k, v in attrib:
            k = escape_attrib(k)
            v = escape_attrib(v)
            output.write("%s:%s;" % (k, v))
        return output.getvalue()
    return ''

_capstyle_d = {'projecting' : 'square', 'butt' : 'butt', 'round': 'round',}
class RendererSVG(RendererBase):
    FONT_SCALE = 100.0
    fontd = maxdict(50)

    def __init__(self, width, height, svgwriter, basename=None, image_dpi=72):
        self.width = width
        self.height = height
        self.writer = XMLWriter(svgwriter)
        self.image_dpi = image_dpi  # the actual dpi we want to rasterize stuff with

        self._groupd = {}
        if not rcParams['svg.image_inline']:
            assert basename is not None
            self.basename = basename
            self._imaged = {}
        self._clipd = OrderedDict()
        self._char_defs = {}
        self._markers = {}
        self._path_collection_id = 0
        self._imaged = {}
        self._hatchd = OrderedDict()
        self._has_gouraud = False
        self._n_gradients = 0
        self._fonts = OrderedDict()
        self.mathtext_parser = MathTextParser('SVG')

        RendererBase.__init__(self)
        self._glyph_map = dict()
        str_height = short_float_fmt(height)
        str_width = short_float_fmt(width)
        svgwriter.write(svgProlog)
        self._start_id = self.writer.start(
            'svg',
            width='%spt' % str_width,
            height='%spt' % str_height,
            viewBox='0 0 %s %s' % (str_width, str_height),
            xmlns="http://www.w3.org/2000/svg",
            version="1.1",
            attrib={'xmlns:xlink': "http://www.w3.org/1999/xlink"})
        self._write_default_style()

    def finalize(self):
        self._write_clips()
        self._write_hatches()
        self._write_svgfonts()
        self.writer.close(self._start_id)
        self.writer.flush()

    def _write_default_style(self):
        writer = self.writer
        default_style = generate_css({
            'stroke-linejoin': 'round',
            'stroke-linecap': 'butt'})
        writer.start('defs')
        writer.start('style', type='text/css')
        writer.data('*{%s}\n' % default_style)
        writer.end('style')
        writer.end('defs')

    def _make_id(self, type, content):
        content = str(content)
        if rcParams['svg.hashsalt'] is None:
            salt = str(uuid.uuid4())
        else:
            salt = rcParams['svg.hashsalt']
        if six.PY3:
            content = content.encode('utf8')
            salt = salt.encode('utf8')
        m = hashlib.md5()
        m.update(salt)
        m.update(content)
        return '%s%s' % (type, m.hexdigest()[:10])

    def _make_flip_transform(self, transform):
        return (transform +
                Affine2D()
                .scale(1.0, -1.0)
                .translate(0.0, self.height))

    def _get_font(self, prop):
        fname = findfont(prop)
        font = get_font(fname)
        font.clear()
        size = prop.get_size_in_points()
        font.set_size(size, 72.0)
        return font

    def _get_hatch(self, gc, rgbFace):
        """
        Create a new hatch pattern
        """
        if rgbFace is not None:
            rgbFace = tuple(rgbFace)
        edge = gc.get_hatch_color()
        if edge is not None:
            edge = tuple(edge)
        dictkey = (gc.get_hatch(), rgbFace, edge)
        oid = self._hatchd.get(dictkey)
        if oid is None:
            oid = self._make_id('h', dictkey)
            self._hatchd[dictkey] = ((gc.get_hatch_path(), rgbFace, edge), oid)
        else:
            _, oid = oid
        return oid

    def _write_hatches(self):
        if not len(self._hatchd):
            return
        HATCH_SIZE = 72
        writer = self.writer
        writer.start('defs')
        for ((path, face, stroke), oid) in six.itervalues(self._hatchd):
            writer.start(
                'pattern',
                id=oid,
                patternUnits="userSpaceOnUse",
                x="0", y="0", width=six.text_type(HATCH_SIZE),
                height=six.text_type(HATCH_SIZE))
            path_data = self._convert_path(
                path,
                Affine2D().scale(HATCH_SIZE).scale(1.0, -1.0).translate(0, HATCH_SIZE),
                simplify=False)
            if face is None:
                fill = 'none'
            else:
                fill = rgb2hex(face)
            writer.element(
                'rect',
                x="0", y="0", width=six.text_type(HATCH_SIZE+1),
                height=six.text_type(HATCH_SIZE+1),
                fill=fill)
            writer.element(
                'path',
                d=path_data,
                style=generate_css({
                    'fill': rgb2hex(stroke),
                    'stroke': rgb2hex(stroke),
                    'stroke-width': six.text_type(rcParams['hatch.linewidth']),
                    'stroke-linecap': 'butt',
                    'stroke-linejoin': 'miter'
                    })
                )
            writer.end('pattern')
        writer.end('defs')

    def _get_style_dict(self, gc, rgbFace):
        """
        return the style string.  style is generated from the
        GraphicsContext and rgbFace
        """
        attrib = {}

        forced_alpha = gc.get_forced_alpha()

        if gc.get_hatch() is not None:
            attrib['fill'] = "url(#%s)" % self._get_hatch(gc, rgbFace)
            if rgbFace is not None and len(rgbFace) == 4 and rgbFace[3] != 1.0 and not forced_alpha:
                attrib['fill-opacity'] = short_float_fmt(rgbFace[3])
        else:
            if rgbFace is None:
                attrib['fill'] = 'none'
            else:
                if tuple(rgbFace[:3]) != (0, 0, 0):
                    attrib['fill'] = rgb2hex(rgbFace)
                if len(rgbFace) == 4 and rgbFace[3] != 1.0 and not forced_alpha:
                    attrib['fill-opacity'] = short_float_fmt(rgbFace[3])

        if forced_alpha and gc.get_alpha() != 1.0:
            attrib['opacity'] = short_float_fmt(gc.get_alpha())

        offset, seq = gc.get_dashes()
        if seq is not None:
            attrib['stroke-dasharray'] = ','.join([short_float_fmt(val) for val in seq])
            attrib['stroke-dashoffset'] = short_float_fmt(float(offset))

        linewidth = gc.get_linewidth()
        if linewidth:
            rgb = gc.get_rgb()
            attrib['stroke'] = rgb2hex(rgb)
            if not forced_alpha and rgb[3] != 1.0:
                attrib['stroke-opacity'] = short_float_fmt(rgb[3])
            if linewidth != 1.0:
                attrib['stroke-width'] = short_float_fmt(linewidth)
            if gc.get_joinstyle() != 'round':
                attrib['stroke-linejoin'] = gc.get_joinstyle()
            if gc.get_capstyle() != 'butt':
                attrib['stroke-linecap'] = _capstyle_d[gc.get_capstyle()]

        return attrib

    def _get_style(self, gc, rgbFace):
        return generate_css(self._get_style_dict(gc, rgbFace))

    def _get_clip(self, gc):
        cliprect = gc.get_clip_rectangle()
        clippath, clippath_trans = gc.get_clip_path()
        if clippath is not None:
            clippath_trans = self._make_flip_transform(clippath_trans)
            dictkey = (id(clippath), str(clippath_trans))
        elif cliprect is not None:
            x, y, w, h = cliprect.bounds
            y = self.height-(y+h)
            dictkey = (x, y, w, h)
        else:
            return None

        clip = self._clipd.get(dictkey)
        if clip is None:
            oid = self._make_id('p', dictkey)
            if clippath is not None:
                self._clipd[dictkey] = ((clippath, clippath_trans), oid)
            else:
                self._clipd[dictkey] = (dictkey, oid)
        else:
            clip, oid = clip
        return oid

    def _write_clips(self):
        if not len(self._clipd):
            return
        writer = self.writer
        writer.start('defs')
        for clip, oid in six.itervalues(self._clipd):
            writer.start('clipPath', id=oid)
            if len(clip) == 2:
                clippath, clippath_trans = clip
                path_data = self._convert_path(clippath, clippath_trans, simplify=False)
                writer.element('path', d=path_data)
            else:
                x, y, w, h = clip
                writer.element(
                    'rect',
                    x=short_float_fmt(x),
                    y=short_float_fmt(y),
                    width=short_float_fmt(w),
                    height=short_float_fmt(h))
            writer.end('clipPath')
        writer.end('defs')

    def _write_svgfonts(self):
        if not rcParams['svg.fonttype'] == 'svgfont':
            return

        writer = self.writer
        writer.start('defs')
        for font_fname, chars in six.iteritems(self._fonts):
            font = get_font(font_fname)
            font.set_size(72, 72)
            sfnt = font.get_sfnt()
            writer.start('font', id=sfnt[1, 0, 0, 4].decode("mac_roman"))
            writer.element(
                'font-face',
                attrib={
                    'font-family': font.family_name,
                    'font-style': font.style_name.lower(),
                    'units-per-em': '72',
                    'bbox': ' '.join(
                        short_float_fmt(x / 64.0) for x in font.bbox)})
            for char in chars:
                glyph = font.load_char(char, flags=LOAD_NO_HINTING)
                verts, codes = font.get_path()
                path = Path(verts, codes)
                path_data = self._convert_path(path)
                # name = font.get_glyph_name(char)
                writer.element(
                    'glyph',
                    d=path_data,
                    attrib={
                        # 'glyph-name': name,
                        'unicode': unichr(char),
                        'horiz-adv-x':
                        short_float_fmt(glyph.linearHoriAdvance / 65536.0)})
            writer.end('font')
        writer.end('defs')

    def open_group(self, s, gid=None):
        """
        Open a grouping element with label *s*. If *gid* is given, use
        *gid* as the id of the group.
        """
        if gid:
            self.writer.start('g', id=gid)
        else:
            self._groupd[s] = self._groupd.get(s, 0) + 1
            self.writer.start('g', id="%s_%d" % (s, self._groupd[s]))

    def close_group(self, s):
        self.writer.end('g')

    def option_image_nocomposite(self):
        """
        return whether to generate a composite image from multiple images on
        a set of axes
        """
        return not rcParams['image.composite_image']

    def _convert_path(self, path, transform=None, clip=None, simplify=None,
                      sketch=None):
        if clip:
            clip = (0.0, 0.0, self.width, self.height)
        else:
            clip = None
        return _path.convert_to_string(
            path, transform, clip, simplify, sketch, 6,
            [b'M', b'L', b'Q', b'C', b'z'], False).decode('ascii')

    def draw_path(self, gc, path, transform, rgbFace=None):
        trans_and_flip = self._make_flip_transform(transform)
        clip = (rgbFace is None and gc.get_hatch_path() is None)
        simplify = path.should_simplify and clip
        path_data = self._convert_path(
            path, trans_and_flip, clip=clip, simplify=simplify,
            sketch=gc.get_sketch_params())

        attrib = {}
        attrib['style'] = self._get_style(gc, rgbFace)

        clipid = self._get_clip(gc)
        if clipid is not None:
            attrib['clip-path'] = 'url(#%s)' % clipid

        if gc.get_url() is not None:
            self.writer.start('a', {'xlink:href': gc.get_url()})
        self.writer.element('path', d=path_data, attrib=attrib)
        if gc.get_url() is not None:
            self.writer.end('a')

    def draw_markers(self, gc, marker_path, marker_trans, path, trans, rgbFace=None):
        if not len(path.vertices):
            return

        writer = self.writer
        path_data = self._convert_path(
            marker_path,
            marker_trans + Affine2D().scale(1.0, -1.0),
            simplify=False)
        style = self._get_style_dict(gc, rgbFace)
        dictkey = (path_data, generate_css(style))
        oid = self._markers.get(dictkey)
        style = generate_css({k: v for k, v in six.iteritems(style)
                              if k.startswith('stroke')})

        if oid is None:
            oid = self._make_id('m', dictkey)
            writer.start('defs')
            writer.element('path', id=oid, d=path_data, style=style)
            writer.end('defs')
            self._markers[dictkey] = oid

        attrib = {}
        clipid = self._get_clip(gc)
        if clipid is not None:
            attrib['clip-path'] = 'url(#%s)' % clipid
        writer.start('g', attrib=attrib)

        trans_and_flip = self._make_flip_transform(trans)
        attrib = {'xlink:href': '#%s' % oid}
        clip = (0, 0, self.width*72, self.height*72)
        for vertices, code in path.iter_segments(
                trans_and_flip, clip=clip, simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                attrib['x'] = short_float_fmt(x)
                attrib['y'] = short_float_fmt(y)
                attrib['style'] = self._get_style(gc, rgbFace)
                writer.element('use', attrib=attrib)
        writer.end('g')

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offsetTrans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position):
        # Is the optimization worth it? Rough calculation:
        # cost of emitting a path in-line is
        #    (len_path + 5) * uses_per_path
        # cost of definition+use is
        #    (len_path + 3) + 9 * uses_per_path
        len_path = len(paths[0].vertices) if len(paths) > 0 else 0
        uses_per_path = self._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        should_do_optimization = \
            len_path + 9 * uses_per_path + 3 < (len_path + 5) * uses_per_path
        if not should_do_optimization:
            return RendererBase.draw_path_collection(
                self, gc, master_transform, paths, all_transforms,
                offsets, offsetTrans, facecolors, edgecolors,
                linewidths, linestyles, antialiaseds, urls,
                offset_position)

        writer = self.writer
        path_codes = []
        writer.start('defs')
        for i, (path, transform) in enumerate(self._iter_collection_raw_paths(
            master_transform, paths, all_transforms)):
            transform = Affine2D(transform.get_matrix()).scale(1.0, -1.0)
            d = self._convert_path(path, transform, simplify=False)
            oid = 'C%x_%x_%s' % (self._path_collection_id, i,
                                  self._make_id('', d))
            writer.element('path', id=oid, d=d)
            path_codes.append(oid)
        writer.end('defs')

        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
            gc, master_transform, all_transforms, path_codes, offsets,
            offsetTrans, facecolors, edgecolors, linewidths, linestyles,
            antialiaseds, urls, offset_position):
            clipid = self._get_clip(gc0)
            url = gc0.get_url()
            if url is not None:
                writer.start('a', attrib={'xlink:href': url})
            if clipid is not None:
                writer.start('g', attrib={'clip-path': 'url(#%s)' % clipid})
            attrib = {
                'xlink:href': '#%s' % path_id,
                'x': short_float_fmt(xo),
                'y': short_float_fmt(self.height - yo),
                'style': self._get_style(gc0, rgbFace)
                }
            writer.element('use', attrib=attrib)
            if clipid is not None:
                writer.end('g')
            if url is not None:
                writer.end('a')

        self._path_collection_id += 1

    def draw_gouraud_triangle(self, gc, points, colors, trans):
        # This uses a method described here:
        #
        #   http://www.svgopen.org/2005/papers/Converting3DFaceToSVG/index.html
        #
        # that uses three overlapping linear gradients to simulate a
        # Gouraud triangle.  Each gradient goes from fully opaque in
        # one corner to fully transparent along the opposite edge.
        # The line between the stop points is perpendicular to the
        # opposite edge.  Underlying these three gradients is a solid
        # triangle whose color is the average of all three points.

        writer = self.writer
        if not self._has_gouraud:
            self._has_gouraud = True
            writer.start(
                'filter',
                id='colorAdd')
            writer.element(
                'feComposite',
                attrib={'in': 'SourceGraphic'},
                in2='BackgroundImage',
                operator='arithmetic',
                k2="1", k3="1")
            writer.end('filter')

        avg_color = np.sum(colors[:, :], axis=0) / 3.0
        # Just skip fully-transparent triangles
        if avg_color[-1] == 0.0:
            return

        trans_and_flip = self._make_flip_transform(trans)
        tpoints = trans_and_flip.transform(points)

        writer.start('defs')
        for i in range(3):
            x1, y1 = tpoints[i]
            x2, y2 = tpoints[(i + 1) % 3]
            x3, y3 = tpoints[(i + 2) % 3]
            c = colors[i][:]

            if x2 == x3:
                xb = x2
                yb = y1
            elif y2 == y3:
                xb = x1
                yb = y2
            else:
                m1 = (y2 - y3) / (x2 - x3)
                b1 = y2 - (m1 * x2)
                m2 = -(1.0 / m1)
                b2 = y1 - (m2 * x1)
                xb = (-b1 + b2) / (m1 - m2)
                yb = m2 * xb + b2

            writer.start(
                'linearGradient',
                id="GR%x_%d" % (self._n_gradients, i),
                x1=short_float_fmt(x1), y1=short_float_fmt(y1),
                x2=short_float_fmt(xb), y2=short_float_fmt(yb))
            writer.element(
                'stop',
                offset='0',
                style=generate_css({'stop-color': rgb2hex(c),
                                    'stop-opacity': short_float_fmt(c[-1])}))
            writer.element(
                'stop',
                offset='1',
                style=generate_css({'stop-color': rgb2hex(c),
                                    'stop-opacity': "0"}))
            writer.end('linearGradient')

        writer.element(
            'polygon',
            id='GT%x' % self._n_gradients,
            points=" ".join([short_float_fmt(x)
                             for x in (x1, y1, x2, y2, x3, y3)]))
        writer.end('defs')

        avg_color = np.sum(colors[:, :], axis=0) / 3.0
        href = '#GT%x' % self._n_gradients
        writer.element(
            'use',
            attrib={'xlink:href': href,
                    'fill': rgb2hex(avg_color),
                    'fill-opacity': short_float_fmt(avg_color[-1])})
        for i in range(3):
            writer.element(
                'use',
                attrib={'xlink:href': href,
                        'fill': 'url(#GR%x_%d)' % (self._n_gradients, i),
                        'fill-opacity': '1',
                        'filter': 'url(#colorAdd)'})

        self._n_gradients += 1

    def draw_gouraud_triangles(self, gc, triangles_array, colors_array,
                               transform):
        attrib = {}
        clipid = self._get_clip(gc)
        if clipid is not None:
            attrib['clip-path'] = 'url(#%s)' % clipid

        self.writer.start('g', attrib=attrib)

        transform = transform.frozen()
        for tri, col in zip(triangles_array, colors_array):
            self.draw_gouraud_triangle(gc, tri, col, transform)

        self.writer.end('g')

    def option_scale_image(self):
        return True

    def get_image_magnification(self):
        return self.image_dpi / 72.0

    def draw_image(self, gc, x, y, im, transform=None):
        h, w = im.shape[:2]

        if w == 0 or h == 0:
            return

        attrib = {}
        clipid = self._get_clip(gc)
        if clipid is not None:
            # Can't apply clip-path directly to the image because the
            # image has a transformation, which would also be applied
            # to the clip-path
            self.writer.start('g', attrib={'clip-path': 'url(#%s)' % clipid})

        oid = gc.get_gid()
        url = gc.get_url()
        if url is not None:
            self.writer.start('a', attrib={'xlink:href': url})
        if rcParams['svg.image_inline']:
            bytesio = io.BytesIO()
            _png.write_png(im, bytesio)
            oid = oid or self._make_id('image', bytesio.getvalue())
            attrib['xlink:href'] = (
                "data:image/png;base64,\n" +
                base64.b64encode(bytesio.getvalue()).decode('ascii'))
        else:
            self._imaged[self.basename] = self._imaged.get(self.basename, 0) + 1
            filename = '%s.image%d.png'%(self.basename, self._imaged[self.basename])
            _log.info('Writing image file for inclusion: %s', filename)
            _png.write_png(im, filename)
            oid = oid or 'Im_' + self._make_id('image', filename)
            attrib['xlink:href'] = filename

        attrib['id'] = oid

        if transform is None:
            w = 72.0 * w / self.image_dpi
            h = 72.0 * h / self.image_dpi

            self.writer.element(
                'image',
                transform=generate_transform([
                    ('scale', (1, -1)), ('translate', (0, -h))]),
                x=short_float_fmt(x),
                y=short_float_fmt(-(self.height - y - h)),
                width=short_float_fmt(w), height=short_float_fmt(h),
                attrib=attrib)
        else:
            alpha = gc.get_alpha()
            if alpha != 1.0:
                attrib['opacity'] = short_float_fmt(alpha)

            flipped = (
                Affine2D().scale(1.0 / w, 1.0 / h) +
                transform +
                Affine2D()
                .translate(x, y)
                .scale(1.0, -1.0)
                .translate(0.0, self.height))

            attrib['transform'] = generate_transform(
                [('matrix', flipped.frozen())])
            self.writer.element(
                'image',
                width=short_float_fmt(w), height=short_float_fmt(h),
                attrib=attrib)

        if url is not None:
            self.writer.end('a')
        if clipid is not None:
            self.writer.end('g')

    def _adjust_char_id(self, char_id):
        return char_id.replace("%20", "_")

    def _draw_text_as_path(self, gc, x, y, s, prop, angle, ismath, mtext=None):
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
        writer = self.writer

        writer.comment(s)

        glyph_map=self._glyph_map

        text2path = self._text2path
        color = rgb2hex(gc.get_rgb())
        fontsize = prop.get_size_in_points()

        style = {}
        if color != '#000000':
            style['fill'] = color
        if gc.get_alpha() != 1.0:
            style['opacity'] = short_float_fmt(gc.get_alpha())

        if not ismath:
            font = text2path._get_font(prop)
            _glyphs = text2path.get_glyphs_with_font(
                font, s, glyph_map=glyph_map, return_new_glyphs_only=True)
            glyph_info, glyph_map_new, rects = _glyphs

            if glyph_map_new:
                writer.start('defs')
                for char_id, glyph_path in six.iteritems(glyph_map_new):
                    path = Path(*glyph_path)
                    path_data = self._convert_path(path, simplify=False)
                    writer.element('path', id=char_id, d=path_data)
                writer.end('defs')

                glyph_map.update(glyph_map_new)

            attrib = {}
            attrib['style'] = generate_css(style)
            font_scale = fontsize / text2path.FONT_SCALE
            attrib['transform'] = generate_transform([
                ('translate', (x, y)),
                ('rotate', (-angle,)),
                ('scale', (font_scale, -font_scale))])

            writer.start('g', attrib=attrib)
            for glyph_id, xposition, yposition, scale in glyph_info:
                attrib={'xlink:href': '#%s' % glyph_id}
                if xposition != 0.0:
                    attrib['x'] = short_float_fmt(xposition)
                if yposition != 0.0:
                    attrib['y'] = short_float_fmt(yposition)
                writer.element(
                    'use',
                    attrib=attrib)

            writer.end('g')
        else:
            if ismath == "TeX":
                _glyphs = text2path.get_glyphs_tex(prop, s, glyph_map=glyph_map,
                                                   return_new_glyphs_only=True)
            else:
                _glyphs = text2path.get_glyphs_mathtext(prop, s, glyph_map=glyph_map,
                                                        return_new_glyphs_only=True)

            glyph_info, glyph_map_new, rects = _glyphs

            # we store the character glyphs w/o flipping. Instead, the
            # coordinate will be flipped when this characters are
            # used.
            if glyph_map_new:
                writer.start('defs')
                for char_id, glyph_path in six.iteritems(glyph_map_new):
                    char_id = self._adjust_char_id(char_id)
                    # Some characters are blank
                    if not len(glyph_path[0]):
                        path_data = ""
                    else:
                        path = Path(*glyph_path)
                        path_data = self._convert_path(path, simplify=False)
                    writer.element('path', id=char_id, d=path_data)
                writer.end('defs')

                glyph_map.update(glyph_map_new)

            attrib = {}
            font_scale = fontsize / text2path.FONT_SCALE
            attrib['style'] = generate_css(style)
            attrib['transform'] = generate_transform([
                ('translate', (x, y)),
                ('rotate', (-angle,)),
                ('scale', (font_scale, -font_scale))])

            writer.start('g', attrib=attrib)
            for char_id, xposition, yposition, scale in glyph_info:
                char_id = self._adjust_char_id(char_id)

                writer.element(
                    'use',
                    transform=generate_transform([
                        ('translate', (xposition, yposition)),
                        ('scale', (scale,)),
                        ]),
                    attrib={'xlink:href': '#%s' % char_id})

            for verts, codes in rects:
                path = Path(verts, codes)
                path_data = self._convert_path(path, simplify=False)
                writer.element('path', d=path_data)

            writer.end('g')

    def _draw_text_as_text(self, gc, x, y, s, prop, angle, ismath, mtext=None):
        writer = self.writer

        color = rgb2hex(gc.get_rgb())
        style = {}
        if color != '#000000':
            style['fill'] = color
        if gc.get_alpha() != 1.0:
            style['opacity'] = short_float_fmt(gc.get_alpha())

        if not ismath:
            font = self._get_font(prop)
            font.set_text(s, 0.0, flags=LOAD_NO_HINTING)

            fontsize = prop.get_size_in_points()

            fontfamily = font.family_name
            fontstyle = prop.get_style()

            attrib = {}
            # Must add "px" to workaround a Firefox bug
            style['font-size'] = short_float_fmt(fontsize) + 'px'
            style['font-family'] = six.text_type(fontfamily)
            style['font-style'] = prop.get_style().lower()
            style['font-weight'] = six.text_type(prop.get_weight()).lower()
            attrib['style'] = generate_css(style)

            if mtext and (angle == 0 or mtext.get_rotation_mode() == "anchor"):
                # If text anchoring can be supported, get the original
                # coordinates and add alignment information.

                # Get anchor coordinates.
                transform = mtext.get_transform()
                ax, ay = transform.transform_point(mtext.get_position())
                ay = self.height - ay

                # Don't do vertical anchor alignment. Most applications do not
                # support 'alignment-baseline' yet. Apply the vertical layout
                # to the anchor point manually for now.
                angle_rad = np.deg2rad(angle)
                dir_vert = np.array([np.sin(angle_rad), np.cos(angle_rad)])
                v_offset = np.dot(dir_vert, [(x - ax), (y - ay)])
                ax = ax + v_offset * dir_vert[0]
                ay = ay + v_offset * dir_vert[1]

                ha_mpl_to_svg = {'left': 'start', 'right': 'end',
                                 'center': 'middle'}
                style['text-anchor'] = ha_mpl_to_svg[mtext.get_ha()]

                attrib['x'] = short_float_fmt(ax)
                attrib['y'] = short_float_fmt(ay)
                attrib['style'] = generate_css(style)
                attrib['transform'] = "rotate(%s, %s, %s)" % (
                    short_float_fmt(-angle),
                    short_float_fmt(ax),
                    short_float_fmt(ay))
                writer.element('text', s, attrib=attrib)
            else:
                attrib['transform'] = generate_transform([
                    ('translate', (x, y)),
                    ('rotate', (-angle,))])

                writer.element('text', s, attrib=attrib)

            if rcParams['svg.fonttype'] == 'svgfont':
                fontset = self._fonts.setdefault(font.fname, set())
                for c in s:
                    fontset.add(ord(c))
        else:
            writer.comment(s)

            width, height, descent, svg_elements, used_characters = \
                   self.mathtext_parser.parse(s, 72, prop)
            svg_glyphs = svg_elements.svg_glyphs
            svg_rects = svg_elements.svg_rects

            attrib = {}
            attrib['style'] = generate_css(style)
            attrib['transform'] = generate_transform([
                ('translate', (x, y)),
                ('rotate', (-angle,))])

            # Apply attributes to 'g', not 'text', because we likely
            # have some rectangles as well with the same style and
            # transformation
            writer.start('g', attrib=attrib)

            writer.start('text')

            # Sort the characters by font, and output one tspan for
            # each
            spans = OrderedDict()
            for font, fontsize, thetext, new_x, new_y, metrics in svg_glyphs:
                style = generate_css({
                    'font-size': short_float_fmt(fontsize) + 'px',
                    'font-family': font.family_name,
                    'font-style': font.style_name.lower(),
                    'font-weight': font.style_name.lower()})
                if thetext == 32:
                    thetext = 0xa0 # non-breaking space
                spans.setdefault(style, []).append((new_x, -new_y, thetext))

            if rcParams['svg.fonttype'] == 'svgfont':
                for font, fontsize, thetext, new_x, new_y, metrics in svg_glyphs:
                    fontset = self._fonts.setdefault(font.fname, set())
                    fontset.add(thetext)

            for style, chars in six.iteritems(spans):
                chars.sort()

                same_y = True
                if len(chars) > 1:
                    last_y = chars[0][1]
                    for i in xrange(1, len(chars)):
                        if chars[i][1] != last_y:
                            same_y = False
                            break
                if same_y:
                    ys = six.text_type(chars[0][1])
                else:
                    ys = ' '.join(six.text_type(c[1]) for c in chars)

                attrib = {
                    'style': style,
                    'x': ' '.join(short_float_fmt(c[0]) for c in chars),
                    'y': ys
                    }

                writer.element(
                    'tspan',
                    ''.join(unichr(c[2]) for c in chars),
                    attrib=attrib)

            writer.end('text')

            if len(svg_rects):
                for x, y, width, height in svg_rects:
                    writer.element(
                        'rect',
                        x=short_float_fmt(x),
                        y=short_float_fmt(-y + height),
                        width=short_float_fmt(width),
                        height=short_float_fmt(height)
                        )

            writer.end('g')

    def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
        self._draw_text_as_path(gc, x, y, s, prop, angle, ismath="TeX")

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        clipid = self._get_clip(gc)
        if clipid is not None:
            # Cannot apply clip-path directly to the text, because
            # is has a transformation
            self.writer.start(
                'g', attrib={'clip-path': 'url(#%s)' % clipid})

        if gc.get_url() is not None:
            self.writer.start('a', {'xlink:href': gc.get_url()})

        if rcParams['svg.fonttype'] == 'path':
            self._draw_text_as_path(gc, x, y, s, prop, angle, ismath, mtext)
        else:
            self._draw_text_as_text(gc, x, y, s, prop, angle, ismath, mtext)

        if gc.get_url() is not None:
            self.writer.end('a')

        if clipid is not None:
            self.writer.end('g')

    def flipy(self):
        return True

    def get_canvas_width_height(self):
        return self.width, self.height

    def get_text_width_height_descent(self, s, prop, ismath):
        return self._text2path.get_text_width_height_descent(s, prop, ismath)


class FigureCanvasSVG(FigureCanvasBase):
    filetypes = {'svg': 'Scalable Vector Graphics',
                 'svgz': 'Scalable Vector Graphics'}

    fixed_dpi = 72

    def print_svg(self, filename, *args, **kwargs):
        with cbook.open_file_cm(filename, "w", encoding="utf-8") as fh:

            filename = getattr(fh, 'name', '')
            if not isinstance(filename, six.string_types):
                filename = ''

            if cbook.file_requires_unicode(fh):
                detach = False
            else:
                if six.PY3:
                    fh = io.TextIOWrapper(fh, 'utf-8')
                else:
                    fh = codecs.getwriter('utf-8')(fh)
                detach = True

            result = self._print_svg(filename, fh, **kwargs)

            # Detach underlying stream from wrapper so that it remains open in
            # the caller.
            if detach:
                if six.PY3:
                    fh.detach()
                else:
                    fh.reset()
                    fh.stream = io.BytesIO()

        return result

    def print_svgz(self, filename, *args, **kwargs):
        with cbook.open_file_cm(filename, "wb") as fh, \
                gzip.GzipFile(mode='w', fileobj=fh) as gzipwriter:
            return self.print_svg(gzipwriter)

    def _print_svg(self, filename, fh, **kwargs):
        image_dpi = kwargs.pop("dpi", 72)
        self.figure.set_dpi(72.0)
        width, height = self.figure.get_size_inches()
        w, h = width * 72, height * 72

        _bbox_inches_restore = kwargs.pop("bbox_inches_restore", None)
        renderer = MixedModeRenderer(
            self.figure, width, height, image_dpi,
            RendererSVG(w, h, fh, filename, image_dpi),
            bbox_inches_restore=_bbox_inches_restore)

        self.figure.draw(renderer)
        renderer.finalize()

    def get_default_filetype(self):
        return 'svg'

class FigureManagerSVG(FigureManagerBase):
    pass


svgProlog = """\
<?xml version="1.0" encoding="utf-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<!-- Created with matplotlib (http://matplotlib.org/) -->
"""


@_Backend.export
class _BackendSVG(_Backend):
    FigureCanvas = FigureCanvasSVG
    FigureManager = FigureManagerSVG
