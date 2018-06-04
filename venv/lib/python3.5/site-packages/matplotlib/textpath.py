# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict

import six
from six.moves import zip

import warnings

import numpy as np

from matplotlib.path import Path
from matplotlib import rcParams
import matplotlib.font_manager as font_manager
from matplotlib.ft2font import KERNING_DEFAULT, LOAD_NO_HINTING
from matplotlib.ft2font import LOAD_TARGET_LIGHT
from matplotlib.mathtext import MathTextParser
import matplotlib.dviread as dviread
from matplotlib.font_manager import FontProperties, get_font
from matplotlib.transforms import Affine2D
from six.moves.urllib.parse import quote as urllib_quote


class TextToPath(object):
    """
    A class that convert a given text to a path using ttf fonts.
    """

    FONT_SCALE = 100.
    DPI = 72

    def __init__(self):
        """
        Initialization
        """
        self.mathtext_parser = MathTextParser('path')
        self.tex_font_map = None

        from matplotlib.cbook import maxdict
        self._ps_fontd = maxdict(50)

        self._texmanager = None

        self._adobe_standard_encoding = None

    def _get_adobe_standard_encoding(self):
        enc_name = dviread.find_tex_file('8a.enc')
        enc = dviread.Encoding(enc_name)
        return {c: i for i, c in enumerate(enc.encoding)}

    def _get_font(self, prop):
        """
        find a ttf font.
        """
        fname = font_manager.findfont(prop)
        font = get_font(fname)
        font.set_size(self.FONT_SCALE, self.DPI)

        return font

    def _get_hinting_flag(self):
        return LOAD_NO_HINTING

    def _get_char_id(self, font, ccode):
        """
        Return a unique id for the given font and character-code set.
        """
        sfnt = font.get_sfnt()
        try:
            ps_name = sfnt[1, 0, 0, 6].decode('mac_roman')
        except KeyError:
            ps_name = sfnt[3, 1, 0x0409, 6].decode('utf-16be')
        char_id = urllib_quote('%s-%x' % (ps_name, ccode))
        return char_id

    def _get_char_id_ps(self, font, ccode):
        """
        Return a unique id for the given font and character-code set (for tex).
        """
        ps_name = font.get_ps_font_info()[2]
        char_id = urllib_quote('%s-%d' % (ps_name, ccode))
        return char_id

    def glyph_to_path(self, font, currx=0.):
        """
        convert the ft2font glyph to vertices and codes.
        """
        verts, codes = font.get_path()
        if currx != 0.0:
            verts[:, 0] += currx
        return verts, codes

    def get_text_width_height_descent(self, s, prop, ismath):
        if rcParams['text.usetex']:
            texmanager = self.get_texmanager()
            fontsize = prop.get_size_in_points()
            w, h, d = texmanager.get_text_width_height_descent(s, fontsize,
                                                               renderer=None)
            return w, h, d

        fontsize = prop.get_size_in_points()
        scale = fontsize / self.FONT_SCALE

        if ismath:
            prop = prop.copy()
            prop.set_size(self.FONT_SCALE)

            width, height, descent, trash, used_characters = \
                self.mathtext_parser.parse(s, 72, prop)
            return width * scale, height * scale, descent * scale

        font = self._get_font(prop)
        font.set_text(s, 0.0, flags=LOAD_NO_HINTING)
        w, h = font.get_width_height()
        w /= 64.0  # convert from subpixels
        h /= 64.0
        d = font.get_descent()
        d /= 64.0
        return w * scale, h * scale, d * scale

    def get_text_path(self, prop, s, ismath=False, usetex=False):
        """
        convert text *s* to path (a tuple of vertices and codes for
        matplotlib.path.Path).

        *prop*
          font property

        *s*
          text to be converted

        *usetex*
          If True, use matplotlib usetex mode.

        *ismath*
          If True, use mathtext parser. Effective only if usetex == False.


        """
        if not usetex:
            if not ismath:
                font = self._get_font(prop)
                glyph_info, glyph_map, rects = self.get_glyphs_with_font(
                                                    font, s)
            else:
                glyph_info, glyph_map, rects = self.get_glyphs_mathtext(
                                                    prop, s)
        else:
            glyph_info, glyph_map, rects = self.get_glyphs_tex(prop, s)

        verts, codes = [], []

        for glyph_id, xposition, yposition, scale in glyph_info:
            verts1, codes1 = glyph_map[glyph_id]
            if len(verts1):
                verts1 = np.array(verts1) * scale + [xposition, yposition]
                verts.extend(verts1)
                codes.extend(codes1)

        for verts1, codes1 in rects:
            verts.extend(verts1)
            codes.extend(codes1)

        return verts, codes

    def get_glyphs_with_font(self, font, s, glyph_map=None,
                             return_new_glyphs_only=False):
        """
        convert the string *s* to vertices and codes using the
        provided ttf font.
        """

        # Mostly copied from backend_svg.py.

        lastgind = None

        currx = 0
        xpositions = []
        glyph_ids = []

        if glyph_map is None:
            glyph_map = OrderedDict()

        if return_new_glyphs_only:
            glyph_map_new = OrderedDict()
        else:
            glyph_map_new = glyph_map

        # I'm not sure if I get kernings right. Needs to be verified. -JJL

        for c in s:
            ccode = ord(c)
            gind = font.get_char_index(ccode)
            if gind is None:
                ccode = ord('?')
                gind = 0

            if lastgind is not None:
                kern = font.get_kerning(lastgind, gind, KERNING_DEFAULT)
            else:
                kern = 0

            glyph = font.load_char(ccode, flags=LOAD_NO_HINTING)
            horiz_advance = (glyph.linearHoriAdvance / 65536.0)

            char_id = self._get_char_id(font, ccode)
            if char_id not in glyph_map:
                glyph_map_new[char_id] = self.glyph_to_path(font)

            currx += (kern / 64.0)

            xpositions.append(currx)
            glyph_ids.append(char_id)

            currx += horiz_advance

            lastgind = gind

        ypositions = [0] * len(xpositions)
        sizes = [1.] * len(xpositions)

        rects = []

        return (list(zip(glyph_ids, xpositions, ypositions, sizes)),
                     glyph_map_new, rects)

    def get_glyphs_mathtext(self, prop, s, glyph_map=None,
                            return_new_glyphs_only=False):
        """
        convert the string *s* to vertices and codes by parsing it with
        mathtext.
        """

        prop = prop.copy()
        prop.set_size(self.FONT_SCALE)

        width, height, descent, glyphs, rects = self.mathtext_parser.parse(
            s, self.DPI, prop)

        if not glyph_map:
            glyph_map = OrderedDict()

        if return_new_glyphs_only:
            glyph_map_new = OrderedDict()
        else:
            glyph_map_new = glyph_map

        xpositions = []
        ypositions = []
        glyph_ids = []
        sizes = []

        currx, curry = 0, 0
        for font, fontsize, ccode, ox, oy in glyphs:
            char_id = self._get_char_id(font, ccode)
            if char_id not in glyph_map:
                font.clear()
                font.set_size(self.FONT_SCALE, self.DPI)
                glyph = font.load_char(ccode, flags=LOAD_NO_HINTING)
                glyph_map_new[char_id] = self.glyph_to_path(font)

            xpositions.append(ox)
            ypositions.append(oy)
            glyph_ids.append(char_id)
            size = fontsize / self.FONT_SCALE
            sizes.append(size)

        myrects = []
        for ox, oy, w, h in rects:
            vert1 = [(ox, oy), (ox, oy + h), (ox + w, oy + h),
                     (ox + w, oy), (ox, oy), (0, 0)]
            code1 = [Path.MOVETO,
                     Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                     Path.CLOSEPOLY]
            myrects.append((vert1, code1))

        return (list(zip(glyph_ids, xpositions, ypositions, sizes)),
                glyph_map_new, myrects)

    def get_texmanager(self):
        """
        return the :class:`matplotlib.texmanager.TexManager` instance
        """
        if self._texmanager is None:
            from matplotlib.texmanager import TexManager
            self._texmanager = TexManager()
        return self._texmanager

    def get_glyphs_tex(self, prop, s, glyph_map=None,
                       return_new_glyphs_only=False):
        """
        convert the string *s* to vertices and codes using matplotlib's usetex
        mode.
        """

        # codes are modstly borrowed from pdf backend.

        texmanager = self.get_texmanager()

        if self.tex_font_map is None:
            self.tex_font_map = dviread.PsfontsMap(
                                    dviread.find_tex_file('pdftex.map'))

        if self._adobe_standard_encoding is None:
            self._adobe_standard_encoding = self._get_adobe_standard_encoding()

        fontsize = prop.get_size_in_points()
        if hasattr(texmanager, "get_dvi"):
            dvifilelike = texmanager.get_dvi(s, self.FONT_SCALE)
            dvi = dviread.DviFromFileLike(dvifilelike, self.DPI)
        else:
            dvifile = texmanager.make_dvi(s, self.FONT_SCALE)
            dvi = dviread.Dvi(dvifile, self.DPI)
        with dvi:
            page = next(iter(dvi))

        if glyph_map is None:
            glyph_map = OrderedDict()

        if return_new_glyphs_only:
            glyph_map_new = OrderedDict()
        else:
            glyph_map_new = glyph_map

        glyph_ids, xpositions, ypositions, sizes = [], [], [], []

        # Gather font information and do some setup for combining
        # characters into strings.
        # oldfont, seq = None, []
        for x1, y1, dvifont, glyph, width in page.text:
            font_and_encoding = self._ps_fontd.get(dvifont.texname)
            font_bunch = self.tex_font_map[dvifont.texname]

            if font_and_encoding is None:
                if font_bunch.filename is None:
                    raise ValueError(
                        ("No usable font file found for %s (%s). "
                         "The font may lack a Type-1 version.")
                        % (font_bunch.psname, dvifont.texname))

                font = get_font(font_bunch.filename)

                for charmap_name, charmap_code in [("ADOBE_CUSTOM",
                                                    1094992451),
                                                   ("ADOBE_STANDARD",
                                                    1094995778)]:
                    try:
                        font.select_charmap(charmap_code)
                    except (ValueError, RuntimeError):
                        pass
                    else:
                        break
                else:
                    charmap_name = ""
                    warnings.warn("No supported encoding in font (%s)." %
                                  font_bunch.filename)

                if charmap_name == "ADOBE_STANDARD" and font_bunch.encoding:
                    enc0 = dviread.Encoding(font_bunch.encoding)
                    enc = {i: self._adobe_standard_encoding.get(c, None)
                           for i, c in enumerate(enc0.encoding)}
                else:
                    enc = {}
                self._ps_fontd[dvifont.texname] = font, enc

            else:
                font, enc = font_and_encoding

            ft2font_flag = LOAD_TARGET_LIGHT

            char_id = self._get_char_id_ps(font, glyph)

            if char_id not in glyph_map:
                font.clear()
                font.set_size(self.FONT_SCALE, self.DPI)
                if enc:
                    charcode = enc.get(glyph, None)
                else:
                    charcode = glyph

                if charcode is not None:
                    glyph0 = font.load_char(charcode, flags=ft2font_flag)
                else:
                    warnings.warn("The glyph (%d) of font (%s) cannot be "
                                  "converted with the encoding. Glyph may "
                                  "be wrong" % (glyph, font_bunch.filename))

                    glyph0 = font.load_char(glyph, flags=ft2font_flag)

                glyph_map_new[char_id] = self.glyph_to_path(font)

            glyph_ids.append(char_id)
            xpositions.append(x1)
            ypositions.append(y1)
            sizes.append(dvifont.size / self.FONT_SCALE)

        myrects = []

        for ox, oy, h, w in page.boxes:
            vert1 = [(ox, oy), (ox + w, oy), (ox + w, oy + h),
                     (ox, oy + h), (ox, oy), (0, 0)]
            code1 = [Path.MOVETO,
                     Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                     Path.CLOSEPOLY]
            myrects.append((vert1, code1))

        return (list(zip(glyph_ids, xpositions, ypositions, sizes)),
                glyph_map_new, myrects)


text_to_path = TextToPath()


class TextPath(Path):
    """
    Create a path from the text.
    """

    def __init__(self, xy, s, size=None, prop=None,
                 _interpolation_steps=1, usetex=False,
                 *kl, **kwargs):
        """
        Create a path from the text. No support for TeX yet. Note that
        it simply is a path, not an artist. You need to use the
        PathPatch (or other artists) to draw this path onto the
        canvas.

        xy : position of the text.
        s : text
        size : font size
        prop : font property
        """

        if prop is None:
            prop = FontProperties()

        if size is None:
            size = prop.get_size_in_points()

        self._xy = xy
        self.set_size(size)

        self._cached_vertices = None

        self._vertices, self._codes = self.text_get_vertices_codes(
                                            prop, s,
                                            usetex=usetex)

        self._should_simplify = False
        self._simplify_threshold = rcParams['path.simplify_threshold']
        self._has_nonfinite = False
        self._interpolation_steps = _interpolation_steps

    def set_size(self, size):
        """
        set the size of the text
        """
        self._size = size
        self._invalid = True

    def get_size(self):
        """
        get the size of the text
        """
        return self._size

    def _get_vertices(self):
        """
        Return the cached path after updating it if necessary.
        """
        self._revalidate_path()
        return self._cached_vertices

    def _get_codes(self):
        """
        Return the codes
        """
        return self._codes

    vertices = property(_get_vertices)
    codes = property(_get_codes)

    def _revalidate_path(self):
        """
        update the path if necessary.

        The path for the text is initially create with the font size
        of FONT_SCALE, and this path is rescaled to other size when
        necessary.

        """
        if (self._invalid or
            (self._cached_vertices is None)):
            tr = Affine2D().scale(
                    self._size / text_to_path.FONT_SCALE,
                    self._size / text_to_path.FONT_SCALE).translate(*self._xy)
            self._cached_vertices = tr.transform(self._vertices)
            self._invalid = False

    def is_math_text(self, s):
        """
        Returns True if the given string *s* contains any mathtext.
        """
        # copied from Text.is_math_text -JJL

        # Did we find an even number of non-escaped dollar signs?
        # If so, treat is as math text.
        dollar_count = s.count(r'$') - s.count(r'\$')
        even_dollars = (dollar_count > 0 and dollar_count % 2 == 0)

        if rcParams['text.usetex']:
            return s, 'TeX'

        if even_dollars:
            return s, True
        else:
            return s.replace(r'\$', '$'), False

    def text_get_vertices_codes(self, prop, s, usetex):
        """
        convert the string *s* to vertices and codes using the
        provided font property *prop*. Mostly copied from
        backend_svg.py.
        """

        if usetex:
            verts, codes = text_to_path.get_text_path(prop, s, usetex=True)
        else:
            clean_line, ismath = self.is_math_text(s)
            verts, codes = text_to_path.get_text_path(prop, clean_line,
                                                      ismath=ismath)

        return verts, codes
