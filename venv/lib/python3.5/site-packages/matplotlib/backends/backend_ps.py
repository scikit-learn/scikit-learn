"""
A PostScript backend, which can produce both PostScript .ps and .eps
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import StringIO

import glob, os, shutil, sys, time, datetime
import io
import logging

from tempfile import mkstemp
from matplotlib import cbook, __version__, rcParams, checkdep_ghostscript
from matplotlib.afm import AFM
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, GraphicsContextBase,
    RendererBase)

from matplotlib.cbook import (get_realpath_and_stat, is_writable_file_like,
                              maxdict, file_requires_unicode)
from matplotlib.compat.subprocess import subprocess

from matplotlib.font_manager import findfont, is_opentype_cff_font, get_font
from matplotlib.ft2font import KERNING_DEFAULT, LOAD_NO_HINTING
from matplotlib.ttconv import convert_ttf_to_ps
from matplotlib.mathtext import MathTextParser
from matplotlib._mathtext_data import uni2type1
from matplotlib.path import Path
from matplotlib import _path
from matplotlib.transforms import Affine2D

from matplotlib.backends.backend_mixed import MixedModeRenderer


import numpy as np
import binascii
import re

_log = logging.getLogger(__name__)

backend_version = 'Level II'

debugPS = 0


class PsBackendHelper(object):

    def __init__(self):
        self._cached = {}

    @property
    def gs_exe(self):
        """
        executable name of ghostscript.
        """
        try:
            return self._cached["gs_exe"]
        except KeyError:
            pass

        gs_exe, gs_version = checkdep_ghostscript()
        if gs_exe is None:
            gs_exe = 'gs'

        self._cached["gs_exe"] = str(gs_exe)
        return str(gs_exe)

    @property
    def gs_version(self):
        """
        version of ghostscript.
        """
        try:
            return self._cached["gs_version"]
        except KeyError:
            pass

        from matplotlib.compat.subprocess import Popen, PIPE
        s = Popen([self.gs_exe, "--version"], stdout=PIPE)
        pipe, stderr = s.communicate()
        if six.PY3:
            ver = pipe.decode('ascii')
        else:
            ver = pipe
        try:
            gs_version = tuple(map(int, ver.strip().split(".")))
        except ValueError:
            # if something went wrong parsing return null version number
            gs_version = (0, 0)
        self._cached["gs_version"] = gs_version
        return gs_version

    @property
    def supports_ps2write(self):
        """
        True if the installed ghostscript supports ps2write device.
        """
        return self.gs_version[0] >= 9

ps_backend_helper = PsBackendHelper()

papersize = {'letter': (8.5,11),
             'legal': (8.5,14),
             'ledger': (11,17),
             'a0': (33.11,46.81),
             'a1': (23.39,33.11),
             'a2': (16.54,23.39),
             'a3': (11.69,16.54),
             'a4': (8.27,11.69),
             'a5': (5.83,8.27),
             'a6': (4.13,5.83),
             'a7': (2.91,4.13),
             'a8': (2.07,2.91),
             'a9': (1.457,2.05),
             'a10': (1.02,1.457),
             'b0': (40.55,57.32),
             'b1': (28.66,40.55),
             'b2': (20.27,28.66),
             'b3': (14.33,20.27),
             'b4': (10.11,14.33),
             'b5': (7.16,10.11),
             'b6': (5.04,7.16),
             'b7': (3.58,5.04),
             'b8': (2.51,3.58),
             'b9': (1.76,2.51),
             'b10': (1.26,1.76)}

def _get_papertype(w, h):
    keys = list(six.iterkeys(papersize))
    keys.sort()
    keys.reverse()
    for key in keys:
        if key.startswith('l'): continue
        pw, ph = papersize[key]
        if (w < pw) and (h < ph): return key
    return 'a0'

def _num_to_str(val):
    if isinstance(val, six.string_types): return val

    ival = int(val)
    if val == ival: return str(ival)

    s = "%1.3f"%val
    s = s.rstrip("0")
    s = s.rstrip(".")
    return s

def _nums_to_str(*args):
    return ' '.join(map(_num_to_str,args))


def quote_ps_string(s):
    "Quote dangerous characters of S for use in a PostScript string constant."
    s = s.replace(b"\\", b"\\\\")
    s = s.replace(b"(", b"\\(")
    s = s.replace(b")", b"\\)")
    s = s.replace(b"'", b"\\251")
    s = s.replace(b"`", b"\\301")
    s = re.sub(br"[^ -~\n]", lambda x: br"\%03o" % ord(x.group()), s)
    return s.decode('ascii')


def _move_path_to_path_or_stream(src, dst):
    """Move the contents of file at *src* to path-or-filelike *dst*.

    If *dst* is a path, the metadata of *src* are *not* copied.
    """
    if is_writable_file_like(dst):
        fh = (io.open(src, 'r', encoding='latin-1')
              if file_requires_unicode(dst)
              else io.open(src, 'rb'))
        with fh:
            shutil.copyfileobj(fh, dst)
    else:
        # Py3: shutil.move(src, dst, copy_function=shutil.copyfile)
        open(dst, 'w').close()
        mode = os.stat(dst).st_mode
        shutil.move(src, dst)
        os.chmod(dst, mode)


class RendererPS(RendererBase):
    """
    The renderer handles all the drawing primitives using a graphics
    context instance that controls the colors/styles.
    """

    afmfontd = maxdict(50)

    def __init__(self, width, height, pswriter, imagedpi=72):
        """
        Although postscript itself is dpi independent, we need to
        imform the image code about a requested dpi to generate high
        res images and them scale them before embeddin them
        """
        RendererBase.__init__(self)
        self.width = width
        self.height = height
        self._pswriter = pswriter
        if rcParams['text.usetex']:
            self.textcnt = 0
            self.psfrag = []
        self.imagedpi = imagedpi

        # current renderer state (None=uninitialised)
        self.color = None
        self.linewidth = None
        self.linejoin = None
        self.linecap = None
        self.linedash = None
        self.fontname = None
        self.fontsize = None
        self._hatches = {}
        self.image_magnification = imagedpi/72.0
        self._clip_paths = {}
        self._path_collection_id = 0

        self.used_characters = {}
        self.mathtext_parser = MathTextParser("PS")

        self._afm_font_dir = os.path.join(
            rcParams['datapath'], 'fonts', 'afm')

    def track_characters(self, font, s):
        """Keeps track of which characters are required from
        each font."""
        realpath, stat_key = get_realpath_and_stat(font.fname)
        used_characters = self.used_characters.setdefault(
            stat_key, (realpath, set()))
        used_characters[1].update([ord(x) for x in s])

    def merge_used_characters(self, other):
        for stat_key, (realpath, charset) in six.iteritems(other):
            used_characters = self.used_characters.setdefault(
                stat_key, (realpath, set()))
            used_characters[1].update(charset)

    def set_color(self, r, g, b, store=1):
        if (r,g,b) != self.color:
            if r==g and r==b:
                self._pswriter.write("%1.3f setgray\n"%r)
            else:
                self._pswriter.write("%1.3f %1.3f %1.3f setrgbcolor\n"%(r,g,b))
            if store: self.color = (r,g,b)

    def set_linewidth(self, linewidth, store=1):
        linewidth = float(linewidth)
        if linewidth != self.linewidth:
            self._pswriter.write("%1.3f setlinewidth\n"%linewidth)
            if store: self.linewidth = linewidth

    def set_linejoin(self, linejoin, store=1):
        if linejoin != self.linejoin:
            self._pswriter.write("%d setlinejoin\n"%linejoin)
            if store: self.linejoin = linejoin

    def set_linecap(self, linecap, store=1):
        if linecap != self.linecap:
            self._pswriter.write("%d setlinecap\n"%linecap)
            if store: self.linecap = linecap

    def set_linedash(self, offset, seq, store=1):
        if self.linedash is not None:
            oldo, oldseq = self.linedash
            if np.array_equal(seq, oldseq) and oldo == offset:
                return

        if seq is not None and len(seq):
            s="[%s] %d setdash\n"%(_nums_to_str(*seq), offset)
            self._pswriter.write(s)
        else:
            self._pswriter.write("[] 0 setdash\n")
        if store:
            self.linedash = (offset, seq)

    def set_font(self, fontname, fontsize, store=1):
        if rcParams['ps.useafm']: return
        if (fontname,fontsize) != (self.fontname,self.fontsize):
            out = ("/%s findfont\n"
                   "%1.3f scalefont\n"
                   "setfont\n" % (fontname, fontsize))

            self._pswriter.write(out)
            if store: self.fontname = fontname
            if store: self.fontsize = fontsize

    def create_hatch(self, hatch):
        sidelen = 72
        if hatch in self._hatches:
            return self._hatches[hatch]
        name = 'H%d' % len(self._hatches)
        linewidth = rcParams['hatch.linewidth']
        pageheight = self.height * 72
        self._pswriter.write("""\
  << /PatternType 1
     /PaintType 2
     /TilingType 2
     /BBox[0 0 %(sidelen)d %(sidelen)d]
     /XStep %(sidelen)d
     /YStep %(sidelen)d

     /PaintProc {
        pop
        %(linewidth)f setlinewidth
""" % locals())
        self._pswriter.write(
            self._convert_path(Path.hatch(hatch), Affine2D().scale(sidelen),
                               simplify=False))
        self._pswriter.write("""\
        fill
        stroke
     } bind
   >>
   matrix
   0.0 %(pageheight)f translate
   makepattern
   /%(name)s exch def
""" % locals())
        self._hatches[hatch] = name
        return name

    def get_canvas_width_height(self):
        'return the canvas width and height in display coords'
        return self.width * 72.0, self.height * 72.0

    def get_text_width_height_descent(self, s, prop, ismath):
        """
        get the width and height in display coords of the string s
        with FontPropertry prop

        """
        if rcParams['text.usetex']:
            texmanager = self.get_texmanager()
            fontsize = prop.get_size_in_points()
            w, h, d = texmanager.get_text_width_height_descent(s, fontsize,
                                                               renderer=self)
            return w, h, d

        if ismath:
            width, height, descent, pswriter, used_characters = \
                self.mathtext_parser.parse(s, 72, prop)
            return width, height, descent

        if rcParams['ps.useafm']:
            if ismath: s = s[1:-1]
            font = self._get_font_afm(prop)
            l,b,w,h,d = font.get_str_bbox_and_descent(s)

            fontsize = prop.get_size_in_points()
            scale = 0.001*fontsize
            w *= scale
            h *= scale
            d *= scale
            return w, h, d

        font = self._get_font_ttf(prop)
        font.set_text(s, 0.0, flags=LOAD_NO_HINTING)
        w, h = font.get_width_height()
        w /= 64.0  # convert from subpixels
        h /= 64.0
        d = font.get_descent()
        d /= 64.0
        return w, h, d

    def flipy(self):
        'return true if small y numbers are top for renderer'
        return False

    def _get_font_afm(self, prop):
        key = hash(prop)
        font = self.afmfontd.get(key)
        if font is None:
            fname = findfont(prop, fontext='afm', directory=self._afm_font_dir)
            if fname is None:
                fname = findfont(
                    "Helvetica", fontext='afm', directory=self._afm_font_dir)
            font = self.afmfontd.get(fname)
            if font is None:
                with io.open(fname, 'rb') as fh:
                    font = AFM(fh)
                self.afmfontd[fname] = font
            self.afmfontd[key] = font
        return font

    def _get_font_ttf(self, prop):
        fname = findfont(prop)
        font = get_font(fname)
        font.clear()
        size = prop.get_size_in_points()
        font.set_size(size, 72.0)
        return font

    def _rgb(self, rgba):
        h, w = rgba.shape[:2]
        rgb = rgba[::-1, :, :3]
        return h, w, rgb.tostring()

    def _hex_lines(self, s, chars_per_line=128):
        s = binascii.b2a_hex(s)
        nhex = len(s)
        lines = []
        for i in range(0,nhex,chars_per_line):
            limit = min(i+chars_per_line, nhex)
            lines.append(s[i:limit])
        return lines

    def get_image_magnification(self):
        """
        Get the factor by which to magnify images passed to draw_image.
        Allows a backend to have images at a different resolution to other
        artists.
        """
        return self.image_magnification

    def option_scale_image(self):
        """
        ps backend support arbitrary scaling of image.
        """
        return True

    def option_image_nocomposite(self):
        """
        return whether to generate a composite image from multiple images on
        a set of axes
        """
        return not rcParams['image.composite_image']

    def _get_image_h_w_bits_command(self, im):
        h, w, bits = self._rgb(im)
        imagecmd = "false 3 colorimage"

        return h, w, bits, imagecmd

    def draw_image(self, gc, x, y, im, transform=None):
        """
        Draw the Image instance into the current axes; x is the
        distance in pixels from the left hand side of the canvas and y
        is the distance from bottom
        """

        h, w, bits, imagecmd = self._get_image_h_w_bits_command(im)
        hexlines = b'\n'.join(self._hex_lines(bits)).decode('ascii')

        if transform is None:
            matrix = "1 0 0 1 0 0"
            xscale = w / self.image_magnification
            yscale = h / self.image_magnification
        else:
            matrix = " ".join(map(str, transform.frozen().to_values()))
            xscale = 1.0
            yscale = 1.0

        figh = self.height * 72

        bbox = gc.get_clip_rectangle()
        clippath, clippath_trans = gc.get_clip_path()

        clip = []
        if bbox is not None:
            clipx,clipy,clipw,cliph = bbox.bounds
            clip.append('%s clipbox' % _nums_to_str(clipw, cliph, clipx, clipy))
        if clippath is not None:
            id = self._get_clip_path(clippath, clippath_trans)
            clip.append('%s' % id)
        clip = '\n'.join(clip)

        ps = """gsave
%(clip)s
%(x)s %(y)s translate
[%(matrix)s] concat
%(xscale)s %(yscale)s scale
/DataString %(w)s string def
%(w)s %(h)s 8 [ %(w)s 0 0 -%(h)s 0 %(h)s ]
{
currentfile DataString readhexstring pop
} bind %(imagecmd)s
%(hexlines)s
grestore
""" % locals()
        self._pswriter.write(ps)

    def _convert_path(self, path, transform, clip=False, simplify=None):
        if clip:
            clip = (0.0, 0.0, self.width * 72.0,
                    self.height * 72.0)
        else:
            clip = None
        return _path.convert_to_string(
            path, transform, clip, simplify, None,
            6, [b'm', b'l', b'', b'c', b'cl'], True).decode('ascii')

    def _get_clip_path(self, clippath, clippath_transform):
        key = (clippath, id(clippath_transform))
        pid = self._clip_paths.get(key)
        if pid is None:
            pid = 'c%x' % len(self._clip_paths)
            ps_cmd = ['/%s {' % pid]
            ps_cmd.append(self._convert_path(clippath, clippath_transform,
                                             simplify=False))
            ps_cmd.extend(['clip', 'newpath', '} bind def\n'])
            self._pswriter.write('\n'.join(ps_cmd))
            self._clip_paths[key] = pid
        return pid

    def draw_path(self, gc, path, transform, rgbFace=None):
        """
        Draws a Path instance using the given affine transform.
        """
        clip = rgbFace is None and gc.get_hatch_path() is None
        simplify = path.should_simplify and clip
        ps = self._convert_path(path, transform, clip=clip, simplify=simplify)
        self._draw_ps(ps, gc, rgbFace)

    def draw_markers(
            self, gc, marker_path, marker_trans, path, trans, rgbFace=None):
        """
        Draw the markers defined by path at each of the positions in x
        and y.  path coordinates are points, x and y coords will be
        transformed by the transform
        """
        if debugPS: self._pswriter.write('% draw_markers \n')

        if rgbFace:
            if len(rgbFace) == 4 and rgbFace[3] == 0:
                return
            if rgbFace[0] == rgbFace[1] == rgbFace[2]:
                ps_color = '%1.3f setgray' % rgbFace[0]
            else:
                ps_color = '%1.3f %1.3f %1.3f setrgbcolor' % rgbFace[:3]

        # construct the generic marker command:
        ps_cmd = ['/o {', 'gsave', 'newpath', 'translate'] # don't want the translate to be global

        lw = gc.get_linewidth()
        stroke = lw != 0.0
        if stroke:
            ps_cmd.append('%.1f setlinewidth' % lw)
            jint = gc.get_joinstyle()
            ps_cmd.append('%d setlinejoin' % jint)
            cint = gc.get_capstyle()
            ps_cmd.append('%d setlinecap' % cint)

        ps_cmd.append(self._convert_path(marker_path, marker_trans,
                                         simplify=False))

        if rgbFace:
            if stroke:
                ps_cmd.append('gsave')
            ps_cmd.extend([ps_color, 'fill'])
            if stroke:
                ps_cmd.append('grestore')

        if stroke:
            ps_cmd.append('stroke')
        ps_cmd.extend(['grestore', '} bind def'])

        for vertices, code in path.iter_segments(
                trans,
                clip=(0, 0, self.width*72, self.height*72),
                simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                ps_cmd.append("%g %g o" % (x, y))

        ps = '\n'.join(ps_cmd)
        self._draw_ps(ps, gc, rgbFace, fill=False, stroke=False)

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offsetTrans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position):
        # Is the optimization worth it? Rough calculation:
        # cost of emitting a path in-line is
        #     (len_path + 2) * uses_per_path
        # cost of definition+use is
        #     (len_path + 3) + 3 * uses_per_path
        len_path = len(paths[0].vertices) if len(paths) > 0 else 0
        uses_per_path = self._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        should_do_optimization = \
            len_path + 3 * uses_per_path + 3 < (len_path + 2) * uses_per_path
        if not should_do_optimization:
            return RendererBase.draw_path_collection(
                self, gc, master_transform, paths, all_transforms,
                offsets, offsetTrans, facecolors, edgecolors,
                linewidths, linestyles, antialiaseds, urls,
                offset_position)

        write = self._pswriter.write

        path_codes = []
        for i, (path, transform) in enumerate(self._iter_collection_raw_paths(
            master_transform, paths, all_transforms)):
            name = 'p%x_%x' % (self._path_collection_id, i)
            ps_cmd = ['/%s {' % name,
                      'newpath', 'translate']
            ps_cmd.append(self._convert_path(path, transform, simplify=False))
            ps_cmd.extend(['} bind def\n'])
            write('\n'.join(ps_cmd))
            path_codes.append(name)

        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
            gc, master_transform, all_transforms, path_codes, offsets,
            offsetTrans, facecolors, edgecolors, linewidths, linestyles,
            antialiaseds, urls, offset_position):
            ps = "%g %g %s" % (xo, yo, path_id)
            self._draw_ps(ps, gc0, rgbFace)

        self._path_collection_id += 1

    def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
        """
        draw a Text instance
        """
        w, h, bl = self.get_text_width_height_descent(s, prop, ismath)
        fontsize = prop.get_size_in_points()
        thetext = 'psmarker%d' % self.textcnt
        color = '%1.3f,%1.3f,%1.3f'% gc.get_rgb()[:3]
        fontcmd = {'sans-serif' : r'{\sffamily %s}',
               'monospace'  : r'{\ttfamily %s}'}.get(
                rcParams['font.family'][0], r'{\rmfamily %s}')
        s = fontcmd % s
        tex = r'\color[rgb]{%s} %s' % (color, s)

        corr = 0#w/2*(fontsize-10)/10
        if rcParams['text.latex.preview']:
            # use baseline alignment!
            pos = _nums_to_str(x-corr, y)
            self.psfrag.append(r'\psfrag{%s}[Bl][Bl][1][%f]{\fontsize{%f}{%f}%s}'%(thetext, angle, fontsize, fontsize*1.25, tex))
        else:
            # stick to the bottom alignment, but this may give incorrect baseline some times.
            pos = _nums_to_str(x-corr, y-bl)
            self.psfrag.append(r'\psfrag{%s}[bl][bl][1][%f]{\fontsize{%f}{%f}%s}'%(thetext, angle, fontsize, fontsize*1.25, tex))

        ps = """\
gsave
%(pos)s moveto
(%(thetext)s)
show
grestore
    """ % locals()

        self._pswriter.write(ps)
        self.textcnt += 1

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        """
        Draw a Text instance.
        """
        # local to avoid repeated attribute lookups
        write = self._pswriter.write
        if debugPS:
            write("% text\n")

        if len(gc.get_rgb()) == 4 and gc.get_rgb()[3] == 0:
            return  # Special handling for fully transparent.

        if ismath=='TeX':
            return self.draw_tex(gc, x, y, s, prop, angle)

        elif ismath:
            return self.draw_mathtext(gc, x, y, s, prop, angle)

        elif rcParams['ps.useafm']:
            self.set_color(*gc.get_rgb())

            font = self._get_font_afm(prop)
            fontname = font.get_fontname()
            fontsize = prop.get_size_in_points()
            scale = 0.001*fontsize

            thisx = 0
            thisy = font.get_str_bbox_and_descent(s)[4] * scale
            last_name = None
            lines = []
            for c in s:
                name = uni2type1.get(ord(c), 'question')
                try:
                    width = font.get_width_from_char_name(name)
                except KeyError:
                    name = 'question'
                    width = font.get_width_char('?')
                if last_name is not None:
                    kern = font.get_kern_dist_from_name(last_name, name)
                else:
                    kern = 0
                last_name = name
                thisx += kern * scale

                lines.append('%f %f m /%s glyphshow'%(thisx, thisy, name))

                thisx += width * scale

            thetext = "\n".join(lines)
            ps = """\
gsave
/%(fontname)s findfont
%(fontsize)s scalefont
setfont
%(x)f %(y)f translate
%(angle)f rotate
%(thetext)s
grestore
    """ % locals()
            self._pswriter.write(ps)

        else:
            font = self._get_font_ttf(prop)
            font.set_text(s, 0, flags=LOAD_NO_HINTING)
            self.track_characters(font, s)

            self.set_color(*gc.get_rgb())
            sfnt = font.get_sfnt()
            try:
                ps_name = sfnt[1, 0, 0, 6].decode('mac_roman')
            except KeyError:
                ps_name = sfnt[3, 1, 0x0409, 6].decode('utf-16be')
            ps_name = ps_name.encode('ascii', 'replace').decode('ascii')
            self.set_font(ps_name, prop.get_size_in_points())

            lastgind = None
            lines = []
            thisx = 0
            thisy = 0
            for c in s:
                ccode = ord(c)
                gind = font.get_char_index(ccode)
                if gind is None:
                    ccode = ord('?')
                    name = '.notdef'
                    gind = 0
                else:
                    name = font.get_glyph_name(gind)
                glyph = font.load_char(ccode, flags=LOAD_NO_HINTING)

                if lastgind is not None:
                    kern = font.get_kerning(lastgind, gind, KERNING_DEFAULT)
                else:
                    kern = 0
                lastgind = gind
                thisx += kern/64.0

                lines.append('%f %f m /%s glyphshow'%(thisx, thisy, name))
                thisx += glyph.linearHoriAdvance/65536.0


            thetext = '\n'.join(lines)
            ps = """gsave
%(x)f %(y)f translate
%(angle)f rotate
%(thetext)s
grestore
""" % locals()
            self._pswriter.write(ps)

    def new_gc(self):
        return GraphicsContextPS()

    def draw_mathtext(self, gc,
        x, y, s, prop, angle):
        """
        Draw the math text using matplotlib.mathtext
        """
        if debugPS:
            self._pswriter.write("% mathtext\n")

        width, height, descent, pswriter, used_characters = \
            self.mathtext_parser.parse(s, 72, prop)
        self.merge_used_characters(used_characters)
        self.set_color(*gc.get_rgb())
        thetext = pswriter.getvalue()
        ps = """gsave
%(x)f %(y)f translate
%(angle)f rotate
%(thetext)s
grestore
""" % locals()
        self._pswriter.write(ps)

    def draw_gouraud_triangle(self, gc, points, colors, trans):
        self.draw_gouraud_triangles(gc, points.reshape((1, 3, 2)),
                                    colors.reshape((1, 3, 4)), trans)

    def draw_gouraud_triangles(self, gc, points, colors, trans):
        assert len(points) == len(colors)
        assert points.ndim == 3
        assert points.shape[1] == 3
        assert points.shape[2] == 2
        assert colors.ndim == 3
        assert colors.shape[1] == 3
        assert colors.shape[2] == 4

        shape = points.shape
        flat_points = points.reshape((shape[0] * shape[1], 2))
        flat_points = trans.transform(flat_points)
        flat_colors = colors.reshape((shape[0] * shape[1], 4))
        points_min = np.min(flat_points, axis=0) - (1 << 12)
        points_max = np.max(flat_points, axis=0) + (1 << 12)
        factor = np.ceil((2 ** 32 - 1) / (points_max - points_min))

        xmin, ymin = points_min
        xmax, ymax = points_max

        streamarr = np.empty(
            (shape[0] * shape[1],),
            dtype=[('flags', 'u1'),
                   ('points', '>u4', (2,)),
                   ('colors', 'u1', (3,))])
        streamarr['flags'] = 0
        streamarr['points'] = (flat_points - points_min) * factor
        streamarr['colors'] = flat_colors[:, :3] * 255.0

        stream = quote_ps_string(streamarr.tostring())

        self._pswriter.write("""
gsave
<< /ShadingType 4
   /ColorSpace [/DeviceRGB]
   /BitsPerCoordinate 32
   /BitsPerComponent 8
   /BitsPerFlag 8
   /AntiAlias true
   /Decode [ %(xmin)f %(xmax)f %(ymin)f %(ymax)f 0 1 0 1 0 1 ]
   /DataSource (%(stream)s)
>>
shfill
grestore
""" % locals())

    def _draw_ps(self, ps, gc, rgbFace, fill=True, stroke=True, command=None):
        """
        Emit the PostScript sniplet 'ps' with all the attributes from 'gc'
        applied.  'ps' must consist of PostScript commands to construct a path.

        The fill and/or stroke kwargs can be set to False if the
        'ps' string already includes filling and/or stroking, in
        which case _draw_ps is just supplying properties and
        clipping.
        """
        # local variable eliminates all repeated attribute lookups
        write = self._pswriter.write
        if debugPS and command:
            write("% "+command+"\n")
        mightstroke = gc.shouldstroke()
        stroke = stroke and mightstroke
        fill = (fill and rgbFace is not None and
                (len(rgbFace) <= 3 or rgbFace[3] != 0.0))
        hatch = gc.get_hatch()

        if mightstroke:
            self.set_linewidth(gc.get_linewidth())
            jint = gc.get_joinstyle()
            self.set_linejoin(jint)
            cint = gc.get_capstyle()
            self.set_linecap(cint)
            self.set_linedash(*gc.get_dashes())
            self.set_color(*gc.get_rgb()[:3])
        write('gsave\n')

        cliprect = gc.get_clip_rectangle()
        if cliprect:
            x,y,w,h=cliprect.bounds
            write('%1.4g %1.4g %1.4g %1.4g clipbox\n' % (w,h,x,y))
        clippath, clippath_trans = gc.get_clip_path()
        if clippath:
            id = self._get_clip_path(clippath, clippath_trans)
            write('%s\n' % id)

        # Jochen, is the strip necessary? - this could be a honking big string
        write(ps.strip())
        write("\n")

        if fill:
            if stroke or hatch:
                write("gsave\n")
            self.set_color(store=0, *rgbFace[:3])
            write("fill\n")
            if stroke or hatch:
                write("grestore\n")

        if hatch:
            hatch_name = self.create_hatch(hatch)
            write("gsave\n")
            write("%f %f %f " % gc.get_hatch_color()[:3])
            write("%s setpattern fill grestore\n" % hatch_name)

        if stroke:
            write("stroke\n")

        write("grestore\n")


class GraphicsContextPS(GraphicsContextBase):
    def get_capstyle(self):
        return {'butt':0,
                'round':1,
                'projecting':2}[GraphicsContextBase.get_capstyle(self)]

    def get_joinstyle(self):
        return {'miter':0,
                'round':1,
                'bevel':2}[GraphicsContextBase.get_joinstyle(self)]

    def shouldstroke(self):
        return (self.get_linewidth() > 0.0 and
                (len(self.get_rgb()) <= 3 or self.get_rgb()[3] != 0.0))


class FigureCanvasPS(FigureCanvasBase):
    _renderer_class = RendererPS

    fixed_dpi = 72

    def draw(self):
        pass

    filetypes = {'ps'  : 'Postscript',
                 'eps' : 'Encapsulated Postscript'}

    def get_default_filetype(self):
        return 'ps'

    def print_ps(self, outfile, *args, **kwargs):
        return self._print_ps(outfile, 'ps', *args, **kwargs)

    def print_eps(self, outfile, *args, **kwargs):
        return self._print_ps(outfile, 'eps', *args, **kwargs)

    def _print_ps(self, outfile, format, *args, **kwargs):
        papertype = kwargs.pop("papertype", rcParams['ps.papersize'])
        papertype = papertype.lower()
        if papertype == 'auto':
            pass
        elif papertype not in papersize:
            raise RuntimeError('%s is not a valid papertype. Use one of %s' %
                               (papertype, ', '.join(papersize)))

        orientation = kwargs.pop("orientation", "portrait").lower()
        if orientation == 'landscape': isLandscape = True
        elif orientation == 'portrait': isLandscape = False
        else: raise RuntimeError('Orientation must be "portrait" or "landscape"')

        self.figure.set_dpi(72) # Override the dpi kwarg
        imagedpi = kwargs.pop("dpi", 72)
        facecolor = kwargs.pop("facecolor", "w")
        edgecolor = kwargs.pop("edgecolor", "w")

        if rcParams['text.usetex']:
            self._print_figure_tex(outfile, format, imagedpi, facecolor, edgecolor,
                                   orientation, isLandscape, papertype,
                                   **kwargs)
        else:
            self._print_figure(outfile, format, imagedpi, facecolor, edgecolor,
                               orientation, isLandscape, papertype,
                               **kwargs)

    def _print_figure(self, outfile, format, dpi=72, facecolor='w', edgecolor='w',
                      orientation='portrait', isLandscape=False, papertype=None,
                      metadata=None, **kwargs):
        """
        Render the figure to hardcopy.  Set the figure patch face and
        edge colors.  This is useful because some of the GUIs have a
        gray figure face color background and you'll probably want to
        override this on hardcopy

        If outfile is a string, it is interpreted as a file name.
        If the extension matches .ep* write encapsulated postscript,
        otherwise write a stand-alone PostScript file.

        If outfile is a file object, a stand-alone PostScript file is
        written into this file object.

        metadata must be a dictionary. Currently, only the value for
        the key 'Creator' is used.
        """
        isEPSF = format == 'eps'
        if isinstance(outfile,
                      (six.string_types, getattr(os, "PathLike", ()),)):
            outfile = title = getattr(os, "fspath", lambda obj: obj)(outfile)
            passed_in_file_object = False
        elif is_writable_file_like(outfile):
            title = None
            passed_in_file_object = True
        else:
            raise ValueError("outfile must be a path or a file-like object")

        # find the appropriate papertype
        width, height = self.figure.get_size_inches()
        if papertype == 'auto':
            if isLandscape: papertype = _get_papertype(height, width)
            else: papertype = _get_papertype(width, height)

        if isLandscape: paperHeight, paperWidth = papersize[papertype]
        else: paperWidth, paperHeight = papersize[papertype]

        if rcParams['ps.usedistiller'] and not papertype == 'auto':
            # distillers will improperly clip eps files if the pagesize is
            # too small
            if width>paperWidth or height>paperHeight:
                if isLandscape:
                    papertype = _get_papertype(height, width)
                    paperHeight, paperWidth = papersize[papertype]
                else:
                    papertype = _get_papertype(width, height)
                    paperWidth, paperHeight = papersize[papertype]

        # center the figure on the paper
        xo = 72*0.5*(paperWidth - width)
        yo = 72*0.5*(paperHeight - height)

        l, b, w, h = self.figure.bbox.bounds
        llx = xo
        lly = yo
        urx = llx + w
        ury = lly + h
        rotation = 0
        if isLandscape:
            llx, lly, urx, ury = lly, llx, ury, urx
            xo, yo = 72*paperHeight - yo, xo
            rotation = 90
        bbox = (llx, lly, urx, ury)

        # generate PostScript code for the figure and store it in a string
        origfacecolor = self.figure.get_facecolor()
        origedgecolor = self.figure.get_edgecolor()
        self.figure.set_facecolor(facecolor)
        self.figure.set_edgecolor(edgecolor)


        dryrun = kwargs.get("dryrun", False)
        if dryrun:
            class NullWriter(object):
                def write(self, *kl, **kwargs):
                    pass

            self._pswriter = NullWriter()
        else:
            self._pswriter = io.StringIO()


        # mixed mode rendering
        _bbox_inches_restore = kwargs.pop("bbox_inches_restore", None)
        ps_renderer = self._renderer_class(width, height, self._pswriter,
                                           imagedpi=dpi)
        renderer = MixedModeRenderer(self.figure,
            width, height, dpi, ps_renderer,
            bbox_inches_restore=_bbox_inches_restore)

        self.figure.draw(renderer)

        if dryrun: # return immediately if dryrun (tightbbox=True)
            return

        self.figure.set_facecolor(origfacecolor)
        self.figure.set_edgecolor(origedgecolor)

        # check for custom metadata
        if metadata is not None and 'Creator' in metadata:
            creator_str = metadata['Creator']
        else:
            creator_str = "matplotlib version " + __version__ + \
                ", http://matplotlib.org/"

        def print_figure_impl(fh):
            # write the PostScript headers
            if isEPSF:
                print("%!PS-Adobe-3.0 EPSF-3.0", file=fh)
            else:
                print("%!PS-Adobe-3.0", file=fh)
            if title:
                print("%%Title: "+title, file=fh)
            print("%%Creator: " + creator_str, file=fh)
            # get source date from SOURCE_DATE_EPOCH, if set
            # See https://reproducible-builds.org/specs/source-date-epoch/
            source_date_epoch = os.getenv("SOURCE_DATE_EPOCH")
            if source_date_epoch:
                source_date = datetime.datetime.utcfromtimestamp(
                    int(source_date_epoch)).strftime("%a %b %d %H:%M:%S %Y")
            else:
                source_date = time.ctime()
            print("%%CreationDate: "+source_date, file=fh)
            print("%%Orientation: " + orientation, file=fh)
            if not isEPSF:
                print("%%DocumentPaperSizes: "+papertype, file=fh)
            print("%%%%BoundingBox: %d %d %d %d" % bbox, file=fh)
            if not isEPSF:
                print("%%Pages: 1", file=fh)
            print("%%EndComments", file=fh)

            Ndict = len(psDefs)
            print("%%BeginProlog", file=fh)
            if not rcParams['ps.useafm']:
                Ndict += len(ps_renderer.used_characters)
            print("/mpldict %d dict def" % Ndict, file=fh)
            print("mpldict begin", file=fh)
            for d in psDefs:
                d = d.strip()
                for l in d.split('\n'):
                    print(l.strip(), file=fh)
            if not rcParams['ps.useafm']:
                for font_filename, chars in six.itervalues(
                        ps_renderer.used_characters):
                    if len(chars):
                        font = get_font(font_filename)
                        glyph_ids = []
                        for c in chars:
                            gind = font.get_char_index(c)
                            glyph_ids.append(gind)

                        fonttype = rcParams['ps.fonttype']

                        # Can not use more than 255 characters from a
                        # single font for Type 3
                        if len(glyph_ids) > 255:
                            fonttype = 42

                        # The ttf to ps (subsetting) support doesn't work for
                        # OpenType fonts that are Postscript inside (like the
                        # STIX fonts).  This will simply turn that off to avoid
                        # errors.
                        if is_opentype_cff_font(font_filename):
                            raise RuntimeError(
                                "OpenType CFF fonts can not be saved using "
                                "the internal Postscript backend at this "
                                "time; consider using the Cairo backend")
                        else:
                            fh.flush()
                            convert_ttf_to_ps(
                                font_filename.encode(
                                    sys.getfilesystemencoding()),
                                fh, fonttype, glyph_ids)
            print("end", file=fh)
            print("%%EndProlog", file=fh)

            if not isEPSF:
                print("%%Page: 1 1", file=fh)
            print("mpldict begin", file=fh)

            print("%s translate" % _nums_to_str(xo, yo), file=fh)
            if rotation:
                print("%d rotate" % rotation, file=fh)
            print("%s clipbox" % _nums_to_str(width*72, height*72, 0, 0),
                  file=fh)

            # write the figure
            content = self._pswriter.getvalue()
            if not isinstance(content, six.text_type):
                content = content.decode('ascii')
            print(content, file=fh)

            # write the trailer
            print("end", file=fh)
            print("showpage", file=fh)
            if not isEPSF:
                print("%%EOF", file=fh)
            fh.flush()

        if rcParams['ps.usedistiller']:
            # We are going to use an external program to process the output.
            # Write to a temporary file.
            fd, tmpfile = mkstemp()
            try:
                with io.open(fd, 'w', encoding='latin-1') as fh:
                    print_figure_impl(fh)
                if rcParams['ps.usedistiller'] == 'ghostscript':
                    gs_distill(tmpfile, isEPSF, ptype=papertype, bbox=bbox)
                elif rcParams['ps.usedistiller'] == 'xpdf':
                    xpdf_distill(tmpfile, isEPSF, ptype=papertype, bbox=bbox)

                _move_path_to_path_or_stream(tmpfile, outfile)
            finally:
                if os.path.isfile(tmpfile):
                    os.unlink(tmpfile)

        else:
            # Write directly to outfile.
            if passed_in_file_object:
                requires_unicode = file_requires_unicode(outfile)

                if (not requires_unicode and
                        (six.PY3 or not isinstance(outfile, StringIO))):
                    fh = io.TextIOWrapper(outfile, encoding="latin-1")

                    # Prevent the io.TextIOWrapper from closing the
                    # underlying file
                    def do_nothing():
                        pass
                    fh.close = do_nothing
                else:
                    fh = outfile

                print_figure_impl(fh)
            else:
                with io.open(outfile, 'w', encoding='latin-1') as fh:
                    print_figure_impl(fh)

    def _print_figure_tex(self, outfile, format, dpi, facecolor, edgecolor,
                          orientation, isLandscape, papertype, metadata=None,
                          **kwargs):
        """
        If text.usetex is True in rc, a temporary pair of tex/eps files
        are created to allow tex to manage the text layout via the PSFrags
        package. These files are processed to yield the final ps or eps file.

        metadata must be a dictionary. Currently, only the value for
        the key 'Creator' is used.
        """
        isEPSF = format == 'eps'
        if isinstance(outfile, six.string_types):
            title = outfile
        elif is_writable_file_like(outfile):
            title = None
        else:
            raise ValueError("outfile must be a path or a file-like object")

        self.figure.dpi = 72  # ignore the dpi kwarg
        width, height = self.figure.get_size_inches()
        xo = 0
        yo = 0

        l, b, w, h = self.figure.bbox.bounds
        llx = xo
        lly = yo
        urx = llx + w
        ury = lly + h
        bbox = (llx, lly, urx, ury)

        # generate PostScript code for the figure and store it in a string
        origfacecolor = self.figure.get_facecolor()
        origedgecolor = self.figure.get_edgecolor()
        self.figure.set_facecolor(facecolor)
        self.figure.set_edgecolor(edgecolor)

        dryrun = kwargs.get("dryrun", False)
        if dryrun:
            class NullWriter(object):
                def write(self, *kl, **kwargs):
                    pass

            self._pswriter = NullWriter()
        else:
            self._pswriter = io.StringIO()

        # mixed mode rendering
        _bbox_inches_restore = kwargs.pop("bbox_inches_restore", None)
        ps_renderer = self._renderer_class(width, height,
                                           self._pswriter, imagedpi=dpi)
        renderer = MixedModeRenderer(self.figure,
                                     width, height, dpi, ps_renderer,
                                     bbox_inches_restore=_bbox_inches_restore)

        self.figure.draw(renderer)

        if dryrun:  # return immediately if dryrun (tightbbox=True)
            return

        self.figure.set_facecolor(origfacecolor)
        self.figure.set_edgecolor(origedgecolor)

        # check for custom metadata
        if metadata is not None and 'Creator' in metadata:
            creator_str = metadata['Creator']
        else:
            creator_str = "matplotlib version " + __version__ + \
                ", http://matplotlib.org/"

        # write to a temp file, we'll move it to outfile when done

        fd, tmpfile = mkstemp()
        try:
            with io.open(fd, 'w', encoding='latin-1') as fh:
                # write the Encapsulated PostScript headers
                print("%!PS-Adobe-3.0 EPSF-3.0", file=fh)
                if title:
                    print("%%Title: "+title, file=fh)
                print("%%Creator: " + creator_str, file=fh)
                # get source date from SOURCE_DATE_EPOCH, if set
                # See https://reproducible-builds.org/specs/source-date-epoch/
                source_date_epoch = os.getenv("SOURCE_DATE_EPOCH")
                if source_date_epoch:
                    source_date = datetime.datetime.utcfromtimestamp(
                        int(source_date_epoch)).strftime(
                            "%a %b %d %H:%M:%S %Y")
                else:
                    source_date = time.ctime()
                print("%%CreationDate: "+source_date, file=fh)
                print("%%%%BoundingBox: %d %d %d %d" % bbox, file=fh)
                print("%%EndComments", file=fh)

                Ndict = len(psDefs)
                print("%%BeginProlog", file=fh)
                print("/mpldict %d dict def" % Ndict, file=fh)
                print("mpldict begin", file=fh)
                for d in psDefs:
                    d = d.strip()
                    for l in d.split('\n'):
                        print(l.strip(), file=fh)
                print("end", file=fh)
                print("%%EndProlog", file=fh)

                print("mpldict begin", file=fh)
                print("%s translate" % _nums_to_str(xo, yo), file=fh)
                print("%s clipbox" % _nums_to_str(width*72, height*72, 0, 0),
                      file=fh)

                # write the figure
                print(self._pswriter.getvalue(), file=fh)

                # write the trailer
                print("end", file=fh)
                print("showpage", file=fh)
                fh.flush()

            if isLandscape:  # now we are ready to rotate
                isLandscape = True
                width, height = height, width
                bbox = (lly, llx, ury, urx)

            # set the paper size to the figure size if isEPSF. The
            # resulting ps file has the given size with correct bounding
            # box so that there is no need to call 'pstoeps'
            if isEPSF:
                paperWidth, paperHeight = self.figure.get_size_inches()
                if isLandscape:
                    paperWidth, paperHeight = paperHeight, paperWidth
            else:
                temp_papertype = _get_papertype(width, height)
                if papertype == 'auto':
                    papertype = temp_papertype
                    paperWidth, paperHeight = papersize[temp_papertype]
                else:
                    paperWidth, paperHeight = papersize[papertype]
                    if (width > paperWidth or height > paperHeight) and isEPSF:
                        paperWidth, paperHeight = papersize[temp_papertype]
                        _log.info('Your figure is too big to fit on %s paper. '
                                  '%s paper will be used to prevent clipping.',
                                  papertype, temp_papertype)

            texmanager = ps_renderer.get_texmanager()
            font_preamble = texmanager.get_font_preamble()
            custom_preamble = texmanager.get_custom_preamble()

            psfrag_rotated = convert_psfrags(tmpfile, ps_renderer.psfrag,
                                             font_preamble,
                                             custom_preamble, paperWidth,
                                             paperHeight,
                                             orientation)

            if (rcParams['ps.usedistiller'] == 'ghostscript'
                    or rcParams['text.usetex']):
                gs_distill(tmpfile, isEPSF, ptype=papertype, bbox=bbox,
                           rotated=psfrag_rotated)
            elif rcParams['ps.usedistiller'] == 'xpdf':
                xpdf_distill(tmpfile, isEPSF, ptype=papertype, bbox=bbox,
                             rotated=psfrag_rotated)

            _move_path_to_path_or_stream(tmpfile, outfile)
        finally:
            if os.path.isfile(tmpfile):
                os.unlink(tmpfile)


def convert_psfrags(tmpfile, psfrags, font_preamble, custom_preamble,
                    paperWidth, paperHeight, orientation):
    """
    When we want to use the LaTeX backend with postscript, we write PSFrag tags
    to a temporary postscript file, each one marking a position for LaTeX to
    render some text. convert_psfrags generates a LaTeX document containing the
    commands to convert those tags to text. LaTeX/dvips produces the postscript
    file that includes the actual text.
    """
    tmpdir = os.path.split(tmpfile)[0]
    epsfile = tmpfile+'.eps'
    shutil.move(tmpfile, epsfile)
    latexfile = tmpfile+'.tex'
    dvifile = tmpfile+'.dvi'
    psfile = tmpfile+'.ps'

    if orientation == 'landscape':
        angle = 90
    else:
        angle = 0

    if rcParams['text.latex.unicode']:
        unicode_preamble = """\\usepackage{ucs}
\\usepackage[utf8x]{inputenc}"""
    else:
        unicode_preamble = ''

    s = """\\documentclass{article}
%s
%s
%s
\\usepackage[dvips, papersize={%sin,%sin}, body={%sin,%sin}, margin={0in,0in}]{geometry}
\\usepackage{psfrag}
\\usepackage[dvips]{graphicx}
\\usepackage{color}
\\pagestyle{empty}
\\begin{document}
\\begin{figure}
\\centering
\\leavevmode
%s
\\includegraphics*[angle=%s]{%s}
\\end{figure}
\\end{document}
"""% (font_preamble, unicode_preamble, custom_preamble, paperWidth, paperHeight,
      paperWidth, paperHeight,
      '\n'.join(psfrags), angle, os.path.split(epsfile)[-1])

    with io.open(latexfile, 'wb') as latexh:
        if rcParams['text.latex.unicode']:
            latexh.write(s.encode('utf8'))
        else:
            try:
                latexh.write(s.encode('ascii'))
            except UnicodeEncodeError:
                _log.info("You are using unicode and latex, but have "
                          "not enabled the matplotlib 'text.latex.unicode' "
                          "rcParam.")
                raise

    # Replace \\ for / so latex does not think there is a function call
    latexfile = latexfile.replace("\\", "/")
    # Replace ~ so Latex does not think it is line break
    latexfile = latexfile.replace("~", "\\string~")
    command = [str("latex"), "-interaction=nonstopmode",
               '"%s"' % latexfile]
    _log.debug('%s', command)
    try:
        report = subprocess.check_output(command, cwd=tmpdir,
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            ('LaTeX was not able to process the following '
             'file:\n%s\n\n'
             'Here is the full report generated by LaTeX:\n%s '
             '\n\n' % (latexfile,
                       exc.output.decode("utf-8"))))
    _log.debug(report)

    command = [str('dvips'), '-q', '-R0', '-o', os.path.basename(psfile),
               os.path.basename(dvifile)]
    _log.debug(command)
    try:
        report = subprocess.check_output(command, cwd=tmpdir,
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            ('dvips was not able to process the following '
             'file:\n%s\n\n'
             'Here is the full report generated by dvips:\n%s '
             '\n\n' % (dvifile,
                       exc.output.decode("utf-8"))))
    _log.debug(report)
    os.remove(epsfile)
    shutil.move(psfile, tmpfile)

    # check if the dvips created a ps in landscape paper.  Somehow,
    # above latex+dvips results in a ps file in a landscape mode for a
    # certain figure sizes (e.g., 8.3in,5.8in which is a5). And the
    # bounding box of the final output got messed up. We check see if
    # the generated ps file is in landscape and return this
    # information. The return value is used in pstoeps step to recover
    # the correct bounding box. 2010-06-05 JJL
    with io.open(tmpfile) as fh:
        if "Landscape" in fh.read(1000):
            psfrag_rotated = True
        else:
            psfrag_rotated = False

    if not debugPS:
        for fname in glob.glob(tmpfile+'.*'):
            os.remove(fname)

    return psfrag_rotated


def gs_distill(tmpfile, eps=False, ptype='letter', bbox=None, rotated=False):
    """
    Use ghostscript's pswrite or epswrite device to distill a file.
    This yields smaller files without illegal encapsulated postscript
    operators. The output is low-level, converting text to outlines.
    """

    if eps:
        paper_option = "-dEPSCrop"
    else:
        paper_option = "-sPAPERSIZE=%s" % ptype

    psfile = tmpfile + '.ps'
    dpi = rcParams['ps.distiller.res']

    gs_exe = ps_backend_helper.gs_exe
    if ps_backend_helper.supports_ps2write: # gs version >= 9
        device_name = "ps2write"
    else:
        device_name = "pswrite"

    command = [str(gs_exe), "-dBATCH", "-dNOPAUSE", "-r%d" % dpi,
               "-sDEVICE=%s" % device_name, paper_option,
               "-sOutputFile=%s" % psfile, tmpfile]
    _log.debug(command)
    try:
        report = subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            ('ghostscript was not able to process your image.\n'
             'Here is the full report generated by ghostscript:\n%s '
             '\n\n' % exc.output.decode("utf-8")))
    _log.debug(report)
    os.remove(tmpfile)
    shutil.move(psfile, tmpfile)


    # While it is best if above steps preserve the original bounding
    # box, there seem to be cases when it is not. For those cases,
    # the original bbox can be restored during the pstoeps step.

    if eps:
        # For some versions of gs, above steps result in an ps file
        # where the original bbox is no more correct. Do not adjust
        # bbox for now.
        if ps_backend_helper.supports_ps2write:
            # fo gs version >= 9 w/ ps2write device
            pstoeps(tmpfile, bbox, rotated=rotated)
        else:
            pstoeps(tmpfile)


def xpdf_distill(tmpfile, eps=False, ptype='letter', bbox=None, rotated=False):
    """
    Use ghostscript's ps2pdf and xpdf's/poppler's pdftops to distill a file.
    This yields smaller files without illegal encapsulated postscript
    operators. This distiller is preferred, generating high-level postscript
    output that treats text as text.
    """
    pdffile = tmpfile + '.pdf'
    psfile = tmpfile + '.ps'

    # Pass options as `-foo#bar` instead of `-foo=bar` to keep Windows happy
    # (https://www.ghostscript.com/doc/9.22/Use.htm#MS_Windows).
    command = [str("ps2pdf"),
               "-dAutoFilterColorImages#false",
               "-dAutoFilterGrayImages#false",
               "-dAutoRotatePages#false",
               "-sGrayImageFilter#FlateEncode",
               "-sColorImageFilter#FlateEncode",
               "-dEPSCrop" if eps else "-sPAPERSIZE#%s" % ptype,
               tmpfile, pdffile]
    _log.debug(command)

    try:
        report = subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            ('ps2pdf was not able to process your image.\n'
             'Here is the full report generated by ps2pdf:\n%s '
             '\n\n' % exc.output.decode("utf-8")))
    _log.debug(report)

    command = [str("pdftops"), "-paper", "match", "-level2", pdffile, psfile]
    _log.debug(command)
    try:
        report = subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            ('pdftops was not able to process your image.\n'
             'Here is the full report generated by pdftops:\n%s '
             '\n\n' % exc.output.decode("utf-8")))
    _log.debug(report)
    os.remove(tmpfile)
    shutil.move(psfile, tmpfile)

    if eps:
        pstoeps(tmpfile)

    for fname in glob.glob(tmpfile+'.*'):
        os.remove(fname)


def get_bbox_header(lbrt, rotated=False):
    """
    return a postscript header stringfor the given bbox lbrt=(l, b, r, t).
    Optionally, return rotate command.
    """

    l, b, r, t = lbrt
    if rotated:
        rotate = "%.2f %.2f  translate\n90 rotate" % (l+r, 0)
    else:
        rotate = ""
    bbox_info = '%%%%BoundingBox: %d %d %d %d' % (l, b, np.ceil(r), np.ceil(t))
    hires_bbox_info = '%%%%HiResBoundingBox: %.6f %.6f %.6f %.6f' % (l, b, r, t)

    return '\n'.join([bbox_info, hires_bbox_info]), rotate


# get_bbox is deprecated. I don't see any reason to use ghostscript to
# find the bounding box, as the required bounding box is alread known.
def get_bbox(tmpfile, bbox):
    """
    Use ghostscript's bbox device to find the center of the bounding box.
    Return an appropriately sized bbox centered around that point. A bit of a
    hack.
    """

    gs_exe = ps_backend_helper.gs_exe
    command = [gs_exe, "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", "%s" % tmpfile]
    _log.debug(command)
    p = subprocess.Popen(command, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         close_fds=True)
    (stdout, stderr) = (p.stdout, p.stderr)
    _log.debug(stdout.read())
    bbox_info = stderr.read()
    _log.info(bbox_info)
    bbox_found = re.search('%%HiResBoundingBox: .*', bbox_info)
    if bbox_found:
        bbox_info = bbox_found.group()
    else:
        raise RuntimeError('Ghostscript was not able to extract a bounding box.\
Here is the Ghostscript output:\n\n%s' % bbox_info)
    l, b, r, t = [float(i) for i in bbox_info.split()[-4:]]

    # this is a hack to deal with the fact that ghostscript does not return the
    # intended bbox, but a tight bbox. For now, we just center the ink in the
    # intended bbox. This is not ideal, users may intend the ink to not be
    # centered.
    if bbox is None:
        l, b, r, t = (l-1, b-1, r+1, t+1)
    else:
        x = (l+r)/2
        y = (b+t)/2
        dx = (bbox[2]-bbox[0])/2
        dy = (bbox[3]-bbox[1])/2
        l,b,r,t = (x-dx, y-dy, x+dx, y+dy)

    bbox_info = '%%%%BoundingBox: %d %d %d %d' % (l, b, np.ceil(r), np.ceil(t))
    hires_bbox_info = '%%%%HiResBoundingBox: %.6f %.6f %.6f %.6f' % (l, b, r, t)

    return '\n'.join([bbox_info, hires_bbox_info])


def pstoeps(tmpfile, bbox=None, rotated=False):
    """
    Convert the postscript to encapsulated postscript.  The bbox of
    the eps file will be replaced with the given *bbox* argument. If
    None, original bbox will be used.
    """

    # if rotated==True, the output eps file need to be rotated
    if bbox:
        bbox_info, rotate = get_bbox_header(bbox, rotated=rotated)
    else:
        bbox_info, rotate = None, None

    epsfile = tmpfile + '.eps'
    with io.open(epsfile, 'wb') as epsh, io.open(tmpfile, 'rb') as tmph:
        write = epsh.write
        # Modify the header:
        for line in tmph:
            if line.startswith(b'%!PS'):
                write(b"%!PS-Adobe-3.0 EPSF-3.0\n")
                if bbox:
                    write(bbox_info.encode('ascii') + b'\n')
            elif line.startswith(b'%%EndComments'):
                write(line)
                write(b'%%BeginProlog\n'
                      b'save\n'
                      b'countdictstack\n'
                      b'mark\n'
                      b'newpath\n'
                      b'/showpage {} def\n'
                      b'/setpagedevice {pop} def\n'
                      b'%%EndProlog\n'
                      b'%%Page 1 1\n')
                if rotate:
                    write(rotate.encode('ascii') + b'\n')
                break
            elif bbox and line.startswith((b'%%Bound', b'%%HiResBound',
                                           b'%%DocumentMedia', b'%%Pages')):
                pass
            else:
                write(line)
        # Now rewrite the rest of the file, and modify the trailer.
        # This is done in a second loop such that the header of the embedded
        # eps file is not modified.
        for line in tmph:
            if line.startswith(b'%%EOF'):
                write(b'cleartomark\n'
                      b'countdictstack\n'
                      b'exch sub { end } repeat\n'
                      b'restore\n'
                      b'showpage\n'
                      b'%%EOF\n')
            elif line.startswith(b'%%PageBoundingBox'):
                pass
            else:
                write(line)

    os.remove(tmpfile)
    shutil.move(epsfile, tmpfile)


class FigureManagerPS(FigureManagerBase):
    pass


# The following Python dictionary psDefs contains the entries for the
# PostScript dictionary mpldict.  This dictionary implements most of
# the matplotlib primitives and some abbreviations.
#
# References:
# http://www.adobe.com/products/postscript/pdfs/PLRM.pdf
# http://www.mactech.com/articles/mactech/Vol.09/09.04/PostscriptTutorial/
# http://www.math.ubc.ca/people/faculty/cass/graphics/text/www/
#

# The usage comments use the notation of the operator summary
# in the PostScript Language reference manual.
psDefs = [
    # x y  *m*  -
    "/m { moveto } bind def",
    # x y  *l*  -
    "/l { lineto } bind def",
    # x y  *r*  -
    "/r { rlineto } bind def",
    # x1 y1 x2 y2 x y *c*  -
    "/c { curveto } bind def",
    # *closepath*  -
    "/cl { closepath } bind def",
    # w h x y  *box*  -
    """/box {
      m
      1 index 0 r
      0 exch r
      neg 0 r
      cl
    } bind def""",
    # w h x y  *clipbox*  -
    """/clipbox {
      box
      clip
      newpath
    } bind def""",
]


@_Backend.export
class _BackendPS(_Backend):
    FigureCanvas = FigureCanvasPS
    FigureManager = FigureManagerPS
