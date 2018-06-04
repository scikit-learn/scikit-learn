# -*- coding: utf-8 -*-

"""
A PDF matplotlib backend
Author: Jouni K Seppänen <jks@iki.fi>
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six import unichr

import codecs
import collections
from datetime import datetime
from functools import total_ordering
from io import BytesIO
import logging
from math import ceil, cos, floor, pi, sin
import os
import re
import struct
import sys
import time
import warnings
import zlib

import numpy as np

from matplotlib import cbook, __version__, rcParams
from matplotlib._pylab_helpers import Gcf
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, GraphicsContextBase,
    RendererBase)
from matplotlib.backends.backend_mixed import MixedModeRenderer
from matplotlib.cbook import (Bunch, get_realpath_and_stat,
                              is_writable_file_like, maxdict)
from matplotlib.figure import Figure
from matplotlib.font_manager import findfont, is_opentype_cff_font, get_font
from matplotlib.afm import AFM
import matplotlib.type1font as type1font
import matplotlib.dviread as dviread
from matplotlib.ft2font import (FIXED_WIDTH, ITALIC, LOAD_NO_SCALE,
                                LOAD_NO_HINTING, KERNING_UNFITTED)
from matplotlib.mathtext import MathTextParser
from matplotlib.transforms import Affine2D, BboxBase
from matplotlib.path import Path
from matplotlib.dates import UTC
from matplotlib import _path
from matplotlib import _png
from matplotlib import ttconv

_log = logging.getLogger(__name__)

# Overview
#
# The low-level knowledge about pdf syntax lies mainly in the pdfRepr
# function and the classes Reference, Name, Operator, and Stream.  The
# PdfFile class knows about the overall structure of pdf documents.
# It provides a "write" method for writing arbitrary strings in the
# file, and an "output" method that passes objects through the pdfRepr
# function before writing them in the file.  The output method is
# called by the RendererPdf class, which contains the various draw_foo
# methods.  RendererPdf contains a GraphicsContextPdf instance, and
# each draw_foo calls self.check_gc before outputting commands.  This
# method checks whether the pdf graphics state needs to be modified
# and outputs the necessary commands.  GraphicsContextPdf represents
# the graphics state, and its "delta" method returns the commands that
# modify the state.

# Add "pdf.use14corefonts: True" in your configuration file to use only
# the 14 PDF core fonts. These fonts do not need to be embedded; every
# PDF viewing application is required to have them. This results in very
# light PDF files you can use directly in LaTeX or ConTeXt documents
# generated with pdfTeX, without any conversion.

# These fonts are: Helvetica, Helvetica-Bold, Helvetica-Oblique,
# Helvetica-BoldOblique, Courier, Courier-Bold, Courier-Oblique,
# Courier-BoldOblique, Times-Roman, Times-Bold, Times-Italic,
# Times-BoldItalic, Symbol, ZapfDingbats.
#
# Some tricky points:
#
# 1. The clip path can only be widened by popping from the state
# stack.  Thus the state must be pushed onto the stack before narrowing
# the clip path.  This is taken care of by GraphicsContextPdf.
#
# 2. Sometimes it is necessary to refer to something (e.g., font,
# image, or extended graphics state, which contains the alpha value)
# in the page stream by a name that needs to be defined outside the
# stream.  PdfFile provides the methods fontName, imageObject, and
# alphaState for this purpose.  The implementations of these methods
# should perhaps be generalized.

# TODOs:
#
# * encoding of fonts, including mathtext fonts and unicode support
# * TTF support has lots of small TODOs, e.g., how do you know if a font
#   is serif/sans-serif, or symbolic/non-symbolic?
# * draw_markers, draw_line_collection, etc.


def fill(strings, linelen=75):
    """Make one string from sequence of strings, with whitespace
    in between. The whitespace is chosen to form lines of at most
    linelen characters, if possible."""
    currpos = 0
    lasti = 0
    result = []
    for i, s in enumerate(strings):
        length = len(s)
        if currpos + length < linelen:
            currpos += length + 1
        else:
            result.append(b' '.join(strings[lasti:i]))
            lasti = i
            currpos = length
    result.append(b' '.join(strings[lasti:]))
    return b'\n'.join(result)

# PDF strings are supposed to be able to include any eight-bit data,
# except that unbalanced parens and backslashes must be escaped by a
# backslash. However, sf bug #2708559 shows that the carriage return
# character may get read as a newline; these characters correspond to
# \gamma and \Omega in TeX's math font encoding. Escaping them fixes
# the bug.
_string_escape_regex = re.compile(br'([\\()\r\n])')


def _string_escape(match):
    m = match.group(0)
    if m in br'\()':
        return b'\\' + m
    elif m == b'\n':
        return br'\n'
    elif m == b'\r':
        return br'\r'
    assert False


def pdfRepr(obj):
    """Map Python objects to PDF syntax."""

    # Some objects defined later have their own pdfRepr method.
    if hasattr(obj, 'pdfRepr'):
        return obj.pdfRepr()

    # Floats. PDF does not have exponential notation (1.0e-10) so we
    # need to use %f with some precision.  Perhaps the precision
    # should adapt to the magnitude of the number?
    elif isinstance(obj, (float, np.floating)):
        if not np.isfinite(obj):
            raise ValueError("Can only output finite numbers in PDF")
        r = ("%.10f" % obj).encode('ascii')
        return r.rstrip(b'0').rstrip(b'.')

    # Booleans. Needs to be tested before integers since
    # isinstance(True, int) is true.
    elif isinstance(obj, bool):
        return [b'false', b'true'][obj]

    # Integers are written as such.
    elif isinstance(obj, (six.integer_types, np.integer)):
        return ("%d" % obj).encode('ascii')

    # Unicode strings are encoded in UTF-16BE with byte-order mark.
    elif isinstance(obj, six.text_type):
        try:
            # But maybe it's really ASCII?
            s = obj.encode('ASCII')
            return pdfRepr(s)
        except UnicodeEncodeError:
            s = codecs.BOM_UTF16_BE + obj.encode('UTF-16BE')
            return pdfRepr(s)

    # Strings are written in parentheses, with backslashes and parens
    # escaped. Actually balanced parens are allowed, but it is
    # simpler to escape them all. TODO: cut long strings into lines;
    # I believe there is some maximum line length in PDF.
    elif isinstance(obj, bytes):
        return b'(' + _string_escape_regex.sub(_string_escape, obj) + b')'

    # Dictionaries. The keys must be PDF names, so if we find strings
    # there, we make Name objects from them. The values may be
    # anything, so the caller must ensure that PDF names are
    # represented as Name objects.
    elif isinstance(obj, dict):
        r = [b"<<"]
        r.extend([Name(key).pdfRepr() + b" " + pdfRepr(obj[key])
                  for key in sorted(obj)])
        r.append(b">>")
        return fill(r)

    # Lists.
    elif isinstance(obj, (list, tuple)):
        r = [b"["]
        r.extend([pdfRepr(val) for val in obj])
        r.append(b"]")
        return fill(r)

    # The null keyword.
    elif obj is None:
        return b'null'

    # A date.
    elif isinstance(obj, datetime):
        r = obj.strftime('D:%Y%m%d%H%M%S')
        z = obj.utcoffset()
        if z is not None:
            z = z.seconds
        else:
            if time.daylight:
                z = time.altzone
            else:
                z = time.timezone
        if z == 0:
            r += 'Z'
        elif z < 0:
            r += "+%02d'%02d'" % ((-z) // 3600, (-z) % 3600)
        else:
            r += "-%02d'%02d'" % (z // 3600, z % 3600)
        return pdfRepr(r)

    # A bounding box
    elif isinstance(obj, BboxBase):
        return fill([pdfRepr(val) for val in obj.bounds])

    else:
        raise TypeError("Don't know a PDF representation for {} objects"
                        .format(type(obj)))


class Reference(object):
    """PDF reference object.
    Use PdfFile.reserveObject() to create References.
    """

    def __init__(self, id):
        self.id = id

    def __repr__(self):
        return "<Reference %d>" % self.id

    def pdfRepr(self):
        return ("%d 0 R" % self.id).encode('ascii')

    def write(self, contents, file):
        write = file.write
        write(("%d 0 obj\n" % self.id).encode('ascii'))
        write(pdfRepr(contents))
        write(b"\nendobj\n")


@total_ordering
class Name(object):
    """PDF name object."""
    __slots__ = ('name',)
    _regex = re.compile(r'[^!-~]')

    def __init__(self, name):
        if isinstance(name, Name):
            self.name = name.name
        else:
            if isinstance(name, bytes):
                name = name.decode('ascii')
            self.name = self._regex.sub(Name.hexify, name).encode('ascii')

    def __repr__(self):
        return "<Name %s>" % self.name

    def __str__(self):
        return '/' + six.text_type(self.name)

    def __eq__(self, other):
        return isinstance(other, Name) and self.name == other.name

    def __lt__(self, other):
        return isinstance(other, Name) and self.name < other.name

    def __hash__(self):
        return hash(self.name)

    @staticmethod
    def hexify(match):
        return '#%02x' % ord(match.group())

    def pdfRepr(self):
        return b'/' + self.name


class Operator(object):
    """PDF operator object."""
    __slots__ = ('op',)

    def __init__(self, op):
        self.op = op

    def __repr__(self):
        return '<Operator %s>' % self.op

    def pdfRepr(self):
        return self.op


class Verbatim(object):
    """Store verbatim PDF command content for later inclusion in the
    stream."""
    def __init__(self, x):
        self._x = x

    def pdfRepr(self):
        return self._x


# PDF operators (not an exhaustive list)
_pdfops = dict(
    close_fill_stroke=b'b', fill_stroke=b'B', fill=b'f', closepath=b'h',
    close_stroke=b's', stroke=b'S', endpath=b'n', begin_text=b'BT',
    end_text=b'ET', curveto=b'c', rectangle=b're', lineto=b'l', moveto=b'm',
    concat_matrix=b'cm', use_xobject=b'Do', setgray_stroke=b'G',
    setgray_nonstroke=b'g', setrgb_stroke=b'RG', setrgb_nonstroke=b'rg',
    setcolorspace_stroke=b'CS', setcolorspace_nonstroke=b'cs',
    setcolor_stroke=b'SCN', setcolor_nonstroke=b'scn', setdash=b'd',
    setlinejoin=b'j', setlinecap=b'J', setgstate=b'gs', gsave=b'q',
    grestore=b'Q', textpos=b'Td', selectfont=b'Tf', textmatrix=b'Tm',
    show=b'Tj', showkern=b'TJ', setlinewidth=b'w', clip=b'W', shading=b'sh')

Op = Bunch(**{name: Operator(value) for name, value in six.iteritems(_pdfops)})


def _paint_path(fill, stroke):
    """Return the PDF operator to paint a path in the following way:
    fill:   fill the path with the fill color
    stroke: stroke the outline of the path with the line color"""
    if stroke:
        if fill:
            return Op.fill_stroke
        else:
            return Op.stroke
    else:
        if fill:
            return Op.fill
        else:
            return Op.endpath
Op.paint_path = _paint_path


class Stream(object):
    """PDF stream object.

    This has no pdfRepr method. Instead, call begin(), then output the
    contents of the stream by calling write(), and finally call end().
    """
    __slots__ = ('id', 'len', 'pdfFile', 'file', 'compressobj', 'extra', 'pos')

    def __init__(self, id, len, file, extra=None, png=None):
        """id: object id of stream; len: an unused Reference object for the
        length of the stream, or None (to use a memory buffer); file:
        a PdfFile; extra: a dictionary of extra key-value pairs to
        include in the stream header; png: if the data is already
        png compressed, the decode parameters"""
        self.id = id            # object id
        self.len = len          # id of length object
        self.pdfFile = file
        self.file = file.fh      # file to which the stream is written
        self.compressobj = None  # compression object
        if extra is None:
            self.extra = dict()
        else:
            self.extra = extra.copy()
        if png is not None:
            self.extra.update({'Filter':      Name('FlateDecode'),
                               'DecodeParms': png})

        self.pdfFile.recordXref(self.id)
        if rcParams['pdf.compression'] and not png:
            self.compressobj = zlib.compressobj(rcParams['pdf.compression'])
        if self.len is None:
            self.file = BytesIO()
        else:
            self._writeHeader()
            self.pos = self.file.tell()

    def _writeHeader(self):
        write = self.file.write
        write(("%d 0 obj\n" % self.id).encode('ascii'))
        dict = self.extra
        dict['Length'] = self.len
        if rcParams['pdf.compression']:
            dict['Filter'] = Name('FlateDecode')

        write(pdfRepr(dict))
        write(b"\nstream\n")

    def end(self):
        """Finalize stream."""

        self._flush()
        if self.len is None:
            contents = self.file.getvalue()
            self.len = len(contents)
            self.file = self.pdfFile.fh
            self._writeHeader()
            self.file.write(contents)
            self.file.write(b"\nendstream\nendobj\n")
        else:
            length = self.file.tell() - self.pos
            self.file.write(b"\nendstream\nendobj\n")
            self.pdfFile.writeObject(self.len, length)

    def write(self, data):
        """Write some data on the stream."""

        if self.compressobj is None:
            self.file.write(data)
        else:
            compressed = self.compressobj.compress(data)
            self.file.write(compressed)

    def _flush(self):
        """Flush the compression object."""

        if self.compressobj is not None:
            compressed = self.compressobj.flush()
            self.file.write(compressed)
            self.compressobj = None


class PdfFile(object):
    """PDF file object."""

    def __init__(self, filename, metadata=None):
        self.nextObject = 1     # next free object id
        self.xrefTable = [[0, 65535, 'the zero object']]
        self.passed_in_file_object = False
        self.original_file_like = None
        self.tell_base = 0
        fh, opened = cbook.to_filehandle(filename, "wb", return_opened=True)
        if not opened:
            try:
                self.tell_base = filename.tell()
            except IOError:
                fh = BytesIO()
                self.original_file_like = filename
            else:
                fh = filename
                self.passed_in_file_object = True

        self._core14fontdir = os.path.join(
            rcParams['datapath'], 'fonts', 'pdfcorefonts')
        self.fh = fh
        self.currentstream = None  # stream object to write to, if any
        fh.write(b"%PDF-1.4\n")    # 1.4 is the first version to have alpha
        # Output some eight-bit chars as a comment so various utilities
        # recognize the file as binary by looking at the first few
        # lines (see note in section 3.4.1 of the PDF reference).
        fh.write(b"%\254\334 \253\272\n")

        self.rootObject = self.reserveObject('root')
        self.pagesObject = self.reserveObject('pages')
        self.pageList = []
        self.fontObject = self.reserveObject('fonts')
        self.alphaStateObject = self.reserveObject('extended graphics states')
        self.hatchObject = self.reserveObject('tiling patterns')
        self.gouraudObject = self.reserveObject('Gouraud triangles')
        self.XObjectObject = self.reserveObject('external objects')
        self.resourceObject = self.reserveObject('resources')

        root = {'Type': Name('Catalog'),
                'Pages': self.pagesObject}
        self.writeObject(self.rootObject, root)

        # get source date from SOURCE_DATE_EPOCH, if set
        # See https://reproducible-builds.org/specs/source-date-epoch/
        source_date_epoch = os.getenv("SOURCE_DATE_EPOCH")
        if source_date_epoch:
            source_date = datetime.utcfromtimestamp(int(source_date_epoch))
            source_date = source_date.replace(tzinfo=UTC)
        else:
            source_date = datetime.today()

        self.infoDict = {
            'Creator': 'matplotlib %s, http://matplotlib.org' % __version__,
            'Producer': 'matplotlib pdf backend %s' % __version__,
            'CreationDate': source_date
        }
        if metadata is not None:
            self.infoDict.update(metadata)
        self.infoDict = {k: v for (k, v) in self.infoDict.items()
                         if v is not None}

        self.fontNames = {}     # maps filenames to internal font names
        self.nextFont = 1       # next free internal font name
        self.dviFontInfo = {}   # maps dvi font names to embedding information
        self._texFontMap = None  # maps TeX font names to PostScript fonts
        # differently encoded Type-1 fonts may share the same descriptor
        self.type1Descriptors = {}
        self.used_characters = {}

        self.alphaStates = {}   # maps alpha values to graphics state objects
        self.nextAlphaState = 1
        # reproducible writeHatches needs an ordered dict:
        self.hatchPatterns = collections.OrderedDict()
        self.nextHatch = 1
        self.gouraudTriangles = []

        self._images = collections.OrderedDict()   # reproducible writeImages
        self.nextImage = 1

        self.markers = collections.OrderedDict()   # reproducible writeMarkers
        self.multi_byte_charprocs = {}

        self.paths = []

        self.pageAnnotations = []  # A list of annotations for the
                                   # current page

        # The PDF spec recommends to include every procset
        procsets = [Name(x)
                    for x in "PDF Text ImageB ImageC ImageI".split()]

        # Write resource dictionary.
        # Possibly TODO: more general ExtGState (graphics state dictionaries)
        #                ColorSpace Pattern Shading Properties
        resources = {'Font': self.fontObject,
                     'XObject': self.XObjectObject,
                     'ExtGState': self.alphaStateObject,
                     'Pattern': self.hatchObject,
                     'Shading': self.gouraudObject,
                     'ProcSet': procsets}
        self.writeObject(self.resourceObject, resources)

    def newPage(self, width, height):
        self.endStream()

        self.width, self.height = width, height
        contentObject = self.reserveObject('page contents')
        thePage = {'Type': Name('Page'),
                   'Parent': self.pagesObject,
                   'Resources': self.resourceObject,
                   'MediaBox': [0, 0, 72 * width, 72 * height],
                   'Contents': contentObject,
                   'Group': {'Type': Name('Group'),
                             'S': Name('Transparency'),
                             'CS': Name('DeviceRGB')},
                   'Annots': self.pageAnnotations,
                   }
        pageObject = self.reserveObject('page')
        self.writeObject(pageObject, thePage)
        self.pageList.append(pageObject)

        self.beginStream(contentObject.id,
                         self.reserveObject('length of content stream'))
        # Initialize the pdf graphics state to match the default mpl
        # graphics context: currently only the join style needs to be set
        self.output(GraphicsContextPdf.joinstyles['round'], Op.setlinejoin)

        # Clear the list of annotations for the next page
        self.pageAnnotations = []

    def newTextnote(self, text, positionRect=[-100, -100, 0, 0]):
        # Create a new annotation of type text
        theNote = {'Type': Name('Annot'),
                   'Subtype': Name('Text'),
                   'Contents': text,
                   'Rect': positionRect,
                   }
        annotObject = self.reserveObject('annotation')
        self.writeObject(annotObject, theNote)
        self.pageAnnotations.append(annotObject)

    def finalize(self):
        "Write out the various deferred objects and the pdf end matter."

        self.endStream()
        self.writeFonts()
        self.writeObject(
            self.alphaStateObject,
            {val[0]: val[1] for val in six.itervalues(self.alphaStates)})
        self.writeHatches()
        self.writeGouraudTriangles()
        xobjects = {
            name: ob for image, name, ob in six.itervalues(self._images)}
        for tup in six.itervalues(self.markers):
            xobjects[tup[0]] = tup[1]
        for name, value in six.iteritems(self.multi_byte_charprocs):
            xobjects[name] = value
        for name, path, trans, ob, join, cap, padding, filled, stroked \
                in self.paths:
            xobjects[name] = ob
        self.writeObject(self.XObjectObject, xobjects)
        self.writeImages()
        self.writeMarkers()
        self.writePathCollectionTemplates()
        self.writeObject(self.pagesObject,
                         {'Type': Name('Pages'),
                          'Kids': self.pageList,
                          'Count': len(self.pageList)})
        self.writeInfoDict()

        # Finalize the file
        self.writeXref()
        self.writeTrailer()

    def close(self):
        "Flush all buffers and free all resources."

        self.endStream()
        if self.passed_in_file_object:
            self.fh.flush()
        else:
            if self.original_file_like is not None:
                self.original_file_like.write(self.fh.getvalue())
            self.fh.close()

    def write(self, data):
        if self.currentstream is None:
            self.fh.write(data)
        else:
            self.currentstream.write(data)

    def output(self, *data):
        self.write(fill([pdfRepr(x) for x in data]))
        self.write(b'\n')

    def beginStream(self, id, len, extra=None, png=None):
        assert self.currentstream is None
        self.currentstream = Stream(id, len, self, extra, png)

    def endStream(self):
        if self.currentstream is not None:
            self.currentstream.end()
            self.currentstream = None

    def fontName(self, fontprop):
        """
        Select a font based on fontprop and return a name suitable for
        Op.selectfont. If fontprop is a string, it will be interpreted
        as the filename of the font.
        """

        if isinstance(fontprop, six.string_types):
            filename = fontprop
        elif rcParams['pdf.use14corefonts']:
            filename = findfont(
                fontprop, fontext='afm', directory=self._core14fontdir)
            if filename is None:
                filename = findfont(
                    "Helvetica", fontext='afm', directory=self._core14fontdir)
        else:
            filename = findfont(fontprop)

        Fx = self.fontNames.get(filename)
        if Fx is None:
            Fx = Name('F%d' % self.nextFont)
            self.fontNames[filename] = Fx
            self.nextFont += 1
            _log.debug('Assigning font %s = %r', Fx, filename)

        return Fx

    @property
    def texFontMap(self):
        # lazy-load texFontMap, it takes a while to parse
        # and usetex is a relatively rare use case
        if self._texFontMap is None:
            self._texFontMap = dviread.PsfontsMap(
                dviread.find_tex_file('pdftex.map'))

        return self._texFontMap

    def dviFontName(self, dvifont):
        """
        Given a dvi font object, return a name suitable for Op.selectfont.
        This registers the font information in self.dviFontInfo if not yet
        registered.
        """

        dvi_info = self.dviFontInfo.get(dvifont.texname)
        if dvi_info is not None:
            return dvi_info.pdfname

        psfont = self.texFontMap[dvifont.texname]
        if psfont.filename is None:
            raise ValueError(
                "No usable font file found for {} (TeX: {}); "
                "the font may lack a Type-1 version"
                .format(psfont.psname, dvifont.texname))

        pdfname = Name('F%d' % self.nextFont)
        self.nextFont += 1
        _log.debug('Assigning font %s = %s (dvi)', pdfname, dvifont.texname)
        self.dviFontInfo[dvifont.texname] = Bunch(
            dvifont=dvifont,
            pdfname=pdfname,
            fontfile=psfont.filename,
            basefont=psfont.psname,
            encodingfile=psfont.encoding,
            effects=psfont.effects)
        return pdfname

    def writeFonts(self):
        fonts = {}
        for dviname, info in sorted(self.dviFontInfo.items()):
            Fx = info.pdfname
            _log.debug('Embedding Type-1 font %s from dvi.', dviname)
            fonts[Fx] = self._embedTeXFont(info)
        for filename in sorted(self.fontNames):
            Fx = self.fontNames[filename]
            _log.debug('Embedding font %s.', filename)
            if filename.endswith('.afm'):
                # from pdf.use14corefonts
                _log.debug('Writing AFM font.')
                fonts[Fx] = self._write_afm_font(filename)
            else:
                # a normal TrueType font
                _log.debug('Writing TrueType font.')
                realpath, stat_key = get_realpath_and_stat(filename)
                chars = self.used_characters.get(stat_key)
                if chars is not None and len(chars[1]):
                    fonts[Fx] = self.embedTTF(realpath, chars[1])
        self.writeObject(self.fontObject, fonts)

    def _write_afm_font(self, filename):
        with open(filename, 'rb') as fh:
            font = AFM(fh)
        fontname = font.get_fontname()
        fontdict = {'Type': Name('Font'),
                    'Subtype': Name('Type1'),
                    'BaseFont': Name(fontname),
                    'Encoding': Name('WinAnsiEncoding')}
        fontdictObject = self.reserveObject('font dictionary')
        self.writeObject(fontdictObject, fontdict)
        return fontdictObject

    def _embedTeXFont(self, fontinfo):
        _log.debug('Embedding TeX font %s - fontinfo=%s',
                   fontinfo.dvifont.texname, fontinfo.__dict__)

        # Widths
        widthsObject = self.reserveObject('font widths')
        self.writeObject(widthsObject, fontinfo.dvifont.widths)

        # Font dictionary
        fontdictObject = self.reserveObject('font dictionary')
        fontdict = {
            'Type':      Name('Font'),
            'Subtype':   Name('Type1'),
            'FirstChar': 0,
            'LastChar':  len(fontinfo.dvifont.widths) - 1,
            'Widths':    widthsObject,
            }

        # Encoding (if needed)
        if fontinfo.encodingfile is not None:
            enc = dviread.Encoding(fontinfo.encodingfile)
            differencesArray = [Name(ch) for ch in enc]
            differencesArray = [0] + differencesArray
            fontdict['Encoding'] = \
                {'Type': Name('Encoding'),
                 'Differences': differencesArray}

        # If no file is specified, stop short
        if fontinfo.fontfile is None:
            _log.warning(
                "Because of TeX configuration (pdftex.map, see updmap option "
                "pdftexDownloadBase14) the font %s is not embedded. This is "
                "deprecated as of PDF 1.5 and it may cause the consumer "
                "application to show something that was not intended.",
                fontinfo.basefont)
            fontdict['BaseFont'] = Name(fontinfo.basefont)
            self.writeObject(fontdictObject, fontdict)
            return fontdictObject

        # We have a font file to embed - read it in and apply any effects
        t1font = type1font.Type1Font(fontinfo.fontfile)
        if fontinfo.effects:
            t1font = t1font.transform(fontinfo.effects)
        fontdict['BaseFont'] = Name(t1font.prop['FontName'])

        # Font descriptors may be shared between differently encoded
        # Type-1 fonts, so only create a new descriptor if there is no
        # existing descriptor for this font.
        effects = (fontinfo.effects.get('slant', 0.0),
                   fontinfo.effects.get('extend', 1.0))
        fontdesc = self.type1Descriptors.get((fontinfo.fontfile, effects))
        if fontdesc is None:
            fontdesc = self.createType1Descriptor(t1font, fontinfo.fontfile)
            self.type1Descriptors[(fontinfo.fontfile, effects)] = fontdesc
        fontdict['FontDescriptor'] = fontdesc

        self.writeObject(fontdictObject, fontdict)
        return fontdictObject

    def createType1Descriptor(self, t1font, fontfile):
        # Create and write the font descriptor and the font file
        # of a Type-1 font
        fontdescObject = self.reserveObject('font descriptor')
        fontfileObject = self.reserveObject('font file')

        italic_angle = t1font.prop['ItalicAngle']
        fixed_pitch = t1font.prop['isFixedPitch']

        flags = 0
        # fixed width
        if fixed_pitch:
            flags |= 1 << 0
        # TODO: serif
        if 0:
            flags |= 1 << 1
        # TODO: symbolic (most TeX fonts are)
        if 1:
            flags |= 1 << 2
        # non-symbolic
        else:
            flags |= 1 << 5
        # italic
        if italic_angle:
            flags |= 1 << 6
        # TODO: all caps
        if 0:
            flags |= 1 << 16
        # TODO: small caps
        if 0:
            flags |= 1 << 17
        # TODO: force bold
        if 0:
            flags |= 1 << 18

        ft2font = get_font(fontfile)

        descriptor = {
            'Type':        Name('FontDescriptor'),
            'FontName':    Name(t1font.prop['FontName']),
            'Flags':       flags,
            'FontBBox':    ft2font.bbox,
            'ItalicAngle': italic_angle,
            'Ascent':      ft2font.ascender,
            'Descent':     ft2font.descender,
            'CapHeight':   1000,  # TODO: find this out
            'XHeight':     500,  # TODO: this one too
            'FontFile':    fontfileObject,
            'FontFamily':  t1font.prop['FamilyName'],
            'StemV':       50,  # TODO
            # (see also revision 3874; but not all TeX distros have AFM files!)
            # 'FontWeight': a number where 400 = Regular, 700 = Bold
            }

        self.writeObject(fontdescObject, descriptor)

        self.beginStream(fontfileObject.id, None,
                         {'Length1': len(t1font.parts[0]),
                          'Length2': len(t1font.parts[1]),
                          'Length3': 0})
        self.currentstream.write(t1font.parts[0])
        self.currentstream.write(t1font.parts[1])
        self.endStream()

        return fontdescObject

    def _get_xobject_symbol_name(self, filename, symbol_name):
        return "%s-%s" % (
            os.path.splitext(os.path.basename(filename))[0],
            symbol_name)

    _identityToUnicodeCMap = """/CIDInit /ProcSet findresource begin
12 dict begin
begincmap
/CIDSystemInfo
<< /Registry (Adobe)
   /Ordering (UCS)
   /Supplement 0
>> def
/CMapName /Adobe-Identity-UCS def
/CMapType 2 def
1 begincodespacerange
<0000> <ffff>
endcodespacerange
%d beginbfrange
%s
endbfrange
endcmap
CMapName currentdict /CMap defineresource pop
end
end"""

    def embedTTF(self, filename, characters):
        """Embed the TTF font from the named file into the document."""

        font = get_font(filename)
        fonttype = rcParams['pdf.fonttype']

        def cvt(length, upe=font.units_per_EM, nearest=True):
            "Convert font coordinates to PDF glyph coordinates"
            value = length / upe * 1000
            if nearest:
                return np.round(value)
            # Perhaps best to round away from zero for bounding
            # boxes and the like
            if value < 0:
                return floor(value)
            else:
                return ceil(value)

        def embedTTFType3(font, characters, descriptor):
            """The Type 3-specific part of embedding a Truetype font"""
            widthsObject = self.reserveObject('font widths')
            fontdescObject = self.reserveObject('font descriptor')
            fontdictObject = self.reserveObject('font dictionary')
            charprocsObject = self.reserveObject('character procs')
            differencesArray = []
            firstchar, lastchar = 0, 255
            bbox = [cvt(x, nearest=False) for x in font.bbox]

            fontdict = {
                'Type': Name('Font'),
                'BaseFont': ps_name,
                'FirstChar': firstchar,
                'LastChar': lastchar,
                'FontDescriptor': fontdescObject,
                'Subtype': Name('Type3'),
                'Name': descriptor['FontName'],
                'FontBBox': bbox,
                'FontMatrix': [.001, 0, 0, .001, 0, 0],
                'CharProcs': charprocsObject,
                'Encoding': {
                    'Type': Name('Encoding'),
                    'Differences': differencesArray},
                'Widths': widthsObject
                }

            # Make the "Widths" array
            from encodings import cp1252
            # The "decoding_map" was changed
            # to a "decoding_table" as of Python 2.5.
            if hasattr(cp1252, 'decoding_map'):
                def decode_char(charcode):
                    return cp1252.decoding_map[charcode] or 0
            else:
                def decode_char(charcode):
                    return ord(cp1252.decoding_table[charcode])

            def get_char_width(charcode):
                s = decode_char(charcode)
                width = font.load_char(
                    s, flags=LOAD_NO_SCALE | LOAD_NO_HINTING).horiAdvance
                return cvt(width)

            widths = [get_char_width(charcode)
                      for charcode in range(firstchar, lastchar+1)]
            descriptor['MaxWidth'] = max(widths)

            # Make the "Differences" array, sort the ccodes < 255 from
            # the multi-byte ccodes, and build the whole set of glyph ids
            # that we need from this font.
            glyph_ids = []
            differences = []
            multi_byte_chars = set()
            for c in characters:
                ccode = c
                gind = font.get_char_index(ccode)
                glyph_ids.append(gind)
                glyph_name = font.get_glyph_name(gind)
                if ccode <= 255:
                    differences.append((ccode, glyph_name))
                else:
                    multi_byte_chars.add(glyph_name)
            differences.sort()

            last_c = -2
            for c, name in differences:
                if c != last_c + 1:
                    differencesArray.append(c)
                differencesArray.append(Name(name))
                last_c = c

            # Make the charprocs array (using ttconv to generate the
            # actual outlines)
            rawcharprocs = ttconv.get_pdf_charprocs(
                filename.encode(sys.getfilesystemencoding()), glyph_ids)
            charprocs = {}
            for charname in sorted(rawcharprocs):
                stream = rawcharprocs[charname]
                charprocDict = {'Length': len(stream)}
                # The 2-byte characters are used as XObjects, so they
                # need extra info in their dictionary
                if charname in multi_byte_chars:
                    charprocDict['Type'] = Name('XObject')
                    charprocDict['Subtype'] = Name('Form')
                    charprocDict['BBox'] = bbox
                    # Each glyph includes bounding box information,
                    # but xpdf and ghostscript can't handle it in a
                    # Form XObject (they segfault!!!), so we remove it
                    # from the stream here.  It's not needed anyway,
                    # since the Form XObject includes it in its BBox
                    # value.
                    stream = stream[stream.find(b"d1") + 2:]
                charprocObject = self.reserveObject('charProc')
                self.beginStream(charprocObject.id, None, charprocDict)
                self.currentstream.write(stream)
                self.endStream()

                # Send the glyphs with ccode > 255 to the XObject dictionary,
                # and the others to the font itself
                if charname in multi_byte_chars:
                    name = self._get_xobject_symbol_name(filename, charname)
                    self.multi_byte_charprocs[name] = charprocObject
                else:
                    charprocs[charname] = charprocObject

            # Write everything out
            self.writeObject(fontdictObject, fontdict)
            self.writeObject(fontdescObject, descriptor)
            self.writeObject(widthsObject, widths)
            self.writeObject(charprocsObject, charprocs)

            return fontdictObject

        def embedTTFType42(font, characters, descriptor):
            """The Type 42-specific part of embedding a Truetype font"""
            fontdescObject = self.reserveObject('font descriptor')
            cidFontDictObject = self.reserveObject('CID font dictionary')
            type0FontDictObject = self.reserveObject('Type 0 font dictionary')
            cidToGidMapObject = self.reserveObject('CIDToGIDMap stream')
            fontfileObject = self.reserveObject('font file stream')
            wObject = self.reserveObject('Type 0 widths')
            toUnicodeMapObject = self.reserveObject('ToUnicode map')

            cidFontDict = {
                'Type': Name('Font'),
                'Subtype': Name('CIDFontType2'),
                'BaseFont': ps_name,
                'CIDSystemInfo': {
                    'Registry': 'Adobe',
                    'Ordering': 'Identity',
                    'Supplement': 0},
                'FontDescriptor': fontdescObject,
                'W': wObject,
                'CIDToGIDMap': cidToGidMapObject
                }

            type0FontDict = {
                'Type': Name('Font'),
                'Subtype': Name('Type0'),
                'BaseFont': ps_name,
                'Encoding': Name('Identity-H'),
                'DescendantFonts': [cidFontDictObject],
                'ToUnicode': toUnicodeMapObject
                }

            # Make fontfile stream
            descriptor['FontFile2'] = fontfileObject
            length1Object = self.reserveObject('decoded length of a font')
            self.beginStream(
                fontfileObject.id,
                self.reserveObject('length of font stream'),
                {'Length1': length1Object})
            with open(filename, 'rb') as fontfile:
                length1 = 0
                while True:
                    data = fontfile.read(4096)
                    if not data:
                        break
                    length1 += len(data)
                    self.currentstream.write(data)
            self.endStream()
            self.writeObject(length1Object, length1)

            # Make the 'W' (Widths) array, CidToGidMap and ToUnicode CMap
            # at the same time
            cid_to_gid_map = ['\0'] * 65536
            widths = []
            max_ccode = 0
            for c in characters:
                ccode = c
                gind = font.get_char_index(ccode)
                glyph = font.load_char(ccode,
                                       flags=LOAD_NO_SCALE | LOAD_NO_HINTING)
                widths.append((ccode, cvt(glyph.horiAdvance)))
                if ccode < 65536:
                    cid_to_gid_map[ccode] = unichr(gind)
                max_ccode = max(ccode, max_ccode)
            widths.sort()
            cid_to_gid_map = cid_to_gid_map[:max_ccode + 1]

            last_ccode = -2
            w = []
            max_width = 0
            unicode_groups = []
            for ccode, width in widths:
                if ccode != last_ccode + 1:
                    w.append(ccode)
                    w.append([width])
                    unicode_groups.append([ccode, ccode])
                else:
                    w[-1].append(width)
                    unicode_groups[-1][1] = ccode
                max_width = max(max_width, width)
                last_ccode = ccode

            unicode_bfrange = []
            for start, end in unicode_groups:
                unicode_bfrange.append(
                    "<%04x> <%04x> [%s]" %
                    (start, end,
                     " ".join(["<%04x>" % x for x in range(start, end+1)])))
            unicode_cmap = (self._identityToUnicodeCMap %
                            (len(unicode_groups),
                             "\n".join(unicode_bfrange))).encode('ascii')

            # CIDToGIDMap stream
            cid_to_gid_map = "".join(cid_to_gid_map).encode("utf-16be")
            self.beginStream(cidToGidMapObject.id,
                             None,
                             {'Length':  len(cid_to_gid_map)})
            self.currentstream.write(cid_to_gid_map)
            self.endStream()

            # ToUnicode CMap
            self.beginStream(toUnicodeMapObject.id,
                             None,
                             {'Length': unicode_cmap})
            self.currentstream.write(unicode_cmap)
            self.endStream()

            descriptor['MaxWidth'] = max_width

            # Write everything out
            self.writeObject(cidFontDictObject, cidFontDict)
            self.writeObject(type0FontDictObject, type0FontDict)
            self.writeObject(fontdescObject, descriptor)
            self.writeObject(wObject, w)

            return type0FontDictObject

        # Beginning of main embedTTF function...

        # You are lost in a maze of TrueType tables, all different...
        sfnt = font.get_sfnt()
        try:
            ps_name = sfnt[1, 0, 0, 6].decode('mac_roman')  # Macintosh scheme
        except KeyError:
            # Microsoft scheme:
            ps_name = sfnt[3, 1, 0x0409, 6].decode('utf-16be')
            # (see freetype/ttnameid.h)
        ps_name = ps_name.encode('ascii', 'replace')
        ps_name = Name(ps_name)
        pclt = font.get_sfnt_table('pclt') or {'capHeight': 0, 'xHeight': 0}
        post = font.get_sfnt_table('post') or {'italicAngle': (0, 0)}
        ff = font.face_flags
        sf = font.style_flags

        flags = 0
        symbolic = False  # ps_name.name in ('Cmsy10', 'Cmmi10', 'Cmex10')
        if ff & FIXED_WIDTH:
            flags |= 1 << 0
        if 0:  # TODO: serif
            flags |= 1 << 1
        if symbolic:
            flags |= 1 << 2
        else:
            flags |= 1 << 5
        if sf & ITALIC:
            flags |= 1 << 6
        if 0:  # TODO: all caps
            flags |= 1 << 16
        if 0:  # TODO: small caps
            flags |= 1 << 17
        if 0:  # TODO: force bold
            flags |= 1 << 18

        descriptor = {
            'Type': Name('FontDescriptor'),
            'FontName': ps_name,
            'Flags': flags,
            'FontBBox': [cvt(x, nearest=False) for x in font.bbox],
            'Ascent': cvt(font.ascender, nearest=False),
            'Descent': cvt(font.descender, nearest=False),
            'CapHeight': cvt(pclt['capHeight'], nearest=False),
            'XHeight': cvt(pclt['xHeight']),
            'ItalicAngle': post['italicAngle'][1],  # ???
            'StemV': 0  # ???
            }

        # The font subsetting to a Type 3 font does not work for
        # OpenType (.otf) that embed a Postscript CFF font, so avoid that --
        # save as a (non-subsetted) Type 42 font instead.
        if is_opentype_cff_font(filename):
            fonttype = 42
            _log.warning("%r can not be subsetted into a Type 3 font. The "
                         "entire font will be embedded in the output.",
                         os.path.basename(filename))

        if fonttype == 3:
            return embedTTFType3(font, characters, descriptor)
        elif fonttype == 42:
            return embedTTFType42(font, characters, descriptor)

    def alphaState(self, alpha):
        """Return name of an ExtGState that sets alpha to the given value."""

        state = self.alphaStates.get(alpha, None)
        if state is not None:
            return state[0]

        name = Name('A%d' % self.nextAlphaState)
        self.nextAlphaState += 1
        self.alphaStates[alpha] = \
            (name, {'Type': Name('ExtGState'),
                    'CA': alpha[0], 'ca': alpha[1]})
        return name

    def hatchPattern(self, hatch_style):
        # The colors may come in as numpy arrays, which aren't hashable
        if hatch_style is not None:
            edge, face, hatch = hatch_style
            if edge is not None:
                edge = tuple(edge)
            if face is not None:
                face = tuple(face)
            hatch_style = (edge, face, hatch)

        pattern = self.hatchPatterns.get(hatch_style, None)
        if pattern is not None:
            return pattern

        name = Name('H%d' % self.nextHatch)
        self.nextHatch += 1
        self.hatchPatterns[hatch_style] = name
        return name

    def writeHatches(self):
        hatchDict = dict()
        sidelen = 72.0
        for hatch_style, name in six.iteritems(self.hatchPatterns):
            ob = self.reserveObject('hatch pattern')
            hatchDict[name] = ob
            res = {'Procsets':
                   [Name(x) for x in "PDF Text ImageB ImageC ImageI".split()]}
            self.beginStream(
                ob.id, None,
                {'Type': Name('Pattern'),
                 'PatternType': 1, 'PaintType': 1, 'TilingType': 1,
                 'BBox': [0, 0, sidelen, sidelen],
                 'XStep': sidelen, 'YStep': sidelen,
                 'Resources': res,
                 # Change origin to match Agg at top-left.
                 'Matrix': [1, 0, 0, 1, 0, self.height * 72]})

            stroke_rgb, fill_rgb, path = hatch_style
            self.output(stroke_rgb[0], stroke_rgb[1], stroke_rgb[2],
                        Op.setrgb_stroke)
            if fill_rgb is not None:
                self.output(fill_rgb[0], fill_rgb[1], fill_rgb[2],
                            Op.setrgb_nonstroke,
                            0, 0, sidelen, sidelen, Op.rectangle,
                            Op.fill)

            self.output(rcParams['hatch.linewidth'], Op.setlinewidth)

            self.output(*self.pathOperations(
                Path.hatch(path),
                Affine2D().scale(sidelen),
                simplify=False))
            self.output(Op.fill_stroke)

            self.endStream()
        self.writeObject(self.hatchObject, hatchDict)

    def addGouraudTriangles(self, points, colors):
        name = Name('GT%d' % len(self.gouraudTriangles))
        self.gouraudTriangles.append((name, points, colors))
        return name

    def writeGouraudTriangles(self):
        gouraudDict = dict()
        for name, points, colors in self.gouraudTriangles:
            ob = self.reserveObject('Gouraud triangle')
            gouraudDict[name] = ob
            shape = points.shape
            flat_points = points.reshape((shape[0] * shape[1], 2))
            flat_colors = colors.reshape((shape[0] * shape[1], 4))
            points_min = np.min(flat_points, axis=0) - (1 << 8)
            points_max = np.max(flat_points, axis=0) + (1 << 8)
            factor = 0xffffffff / (points_max - points_min)

            self.beginStream(
                ob.id, None,
                {'ShadingType': 4,
                 'BitsPerCoordinate': 32,
                 'BitsPerComponent': 8,
                 'BitsPerFlag': 8,
                 'ColorSpace': Name('DeviceRGB'),
                 'AntiAlias': True,
                 'Decode': [points_min[0], points_max[0],
                            points_min[1], points_max[1],
                            0, 1, 0, 1, 0, 1]
                 })

            streamarr = np.empty(
                (shape[0] * shape[1],),
                dtype=[(str('flags'), str('u1')),
                       (str('points'), str('>u4'), (2,)),
                       (str('colors'), str('u1'), (3,))])
            streamarr['flags'] = 0
            streamarr['points'] = (flat_points - points_min) * factor
            streamarr['colors'] = flat_colors[:, :3] * 255.0

            self.write(streamarr.tostring())
            self.endStream()
        self.writeObject(self.gouraudObject, gouraudDict)

    def imageObject(self, image):
        """Return name of an image XObject representing the given image."""

        entry = self._images.get(id(image), None)
        if entry is not None:
            return entry[1]

        name = Name('I%d' % self.nextImage)
        ob = self.reserveObject('image %d' % self.nextImage)
        self.nextImage += 1
        self._images[id(image)] = (image, name, ob)
        return name

    def _unpack(self, im):
        """
        Unpack the image object im into height, width, data, alpha,
        where data and alpha are HxWx3 (RGB) or HxWx1 (grayscale or alpha)
        arrays, except alpha is None if the image is fully opaque.
        """
        h, w = im.shape[:2]
        im = im[::-1]
        if im.ndim == 2:
            return h, w, im, None
        else:
            rgb = im[:, :, :3]
            rgb = np.array(rgb, order='C')
            # PDF needs a separate alpha image
            if im.shape[2] == 4:
                alpha = im[:, :, 3][..., None]
                if np.all(alpha == 255):
                    alpha = None
                else:
                    alpha = np.array(alpha, order='C')
            else:
                alpha = None
            return h, w, rgb, alpha

    def _writePng(self, data):
        """
        Write the image *data* into the pdf file using png
        predictors with Flate compression.
        """

        buffer = BytesIO()
        _png.write_png(data, buffer)
        buffer.seek(8)
        written = 0
        header = bytearray(8)
        while True:
            n = buffer.readinto(header)
            assert n == 8
            length, type = struct.unpack(b'!L4s', bytes(header))
            if type == b'IDAT':
                data = bytearray(length)
                n = buffer.readinto(data)
                assert n == length
                self.currentstream.write(bytes(data))
                written += n
            elif type == b'IEND':
                break
            else:
                buffer.seek(length, 1)
            buffer.seek(4, 1)   # skip CRC

    def _writeImg(self, data, height, width, grayscale, id, smask=None):
        """
        Write the image *data* of size *height* x *width*, as grayscale
        if *grayscale* is true and RGB otherwise, as pdf object *id*
        and with the soft mask (alpha channel) *smask*, which should be
        either None or a *height* x *width* x 1 array.
        """

        obj = {'Type':             Name('XObject'),
               'Subtype':          Name('Image'),
               'Width':            width,
               'Height':           height,
               'ColorSpace':       Name('DeviceGray' if grayscale
                                        else 'DeviceRGB'),
               'BitsPerComponent': 8}
        if smask:
            obj['SMask'] = smask
        if rcParams['pdf.compression']:
            png = {'Predictor': 10,
                   'Colors':    1 if grayscale else 3,
                   'Columns':   width}
        else:
            png = None
        self.beginStream(
            id,
            self.reserveObject('length of image stream'),
            obj,
            png=png
            )
        if png:
            self._writePng(data)
        else:
            self.currentstream.write(data.tostring())
        self.endStream()

    def writeImages(self):
        for img, name, ob in six.itervalues(self._images):
            height, width, data, adata = self._unpack(img)
            if adata is not None:
                smaskObject = self.reserveObject("smask")
                self._writeImg(adata, height, width, True, smaskObject.id)
            else:
                smaskObject = None
            self._writeImg(data, height, width, False,
                           ob.id, smaskObject)

    def markerObject(self, path, trans, fill, stroke, lw, joinstyle,
                     capstyle):
        """Return name of a marker XObject representing the given path."""
        # self.markers used by markerObject, writeMarkers, close:
        # mapping from (path operations, fill?, stroke?) to
        #   [name, object reference, bounding box, linewidth]
        # This enables different draw_markers calls to share the XObject
        # if the gc is sufficiently similar: colors etc can vary, but
        # the choices of whether to fill and whether to stroke cannot.
        # We need a bounding box enclosing all of the XObject path,
        # but since line width may vary, we store the maximum of all
        # occurring line widths in self.markers.
        # close() is somewhat tightly coupled in that it expects the
        # first two components of each value in self.markers to be the
        # name and object reference.
        pathops = self.pathOperations(path, trans, simplify=False)
        key = (tuple(pathops), bool(fill), bool(stroke), joinstyle, capstyle)
        result = self.markers.get(key)
        if result is None:
            name = Name('M%d' % len(self.markers))
            ob = self.reserveObject('marker %d' % len(self.markers))
            bbox = path.get_extents(trans)
            self.markers[key] = [name, ob, bbox, lw]
        else:
            if result[-1] < lw:
                result[-1] = lw
            name = result[0]
        return name

    def writeMarkers(self):
        for ((pathops, fill, stroke, joinstyle, capstyle),
             (name, ob, bbox, lw)) in six.iteritems(self.markers):
            bbox = bbox.padded(lw * 0.5)
            self.beginStream(
                ob.id, None,
                {'Type': Name('XObject'), 'Subtype': Name('Form'),
                 'BBox': list(bbox.extents)})
            self.output(GraphicsContextPdf.joinstyles[joinstyle],
                        Op.setlinejoin)
            self.output(GraphicsContextPdf.capstyles[capstyle], Op.setlinecap)
            self.output(*pathops)
            self.output(Op.paint_path(fill, stroke))
            self.endStream()

    def pathCollectionObject(self, gc, path, trans, padding, filled, stroked):
        name = Name('P%d' % len(self.paths))
        ob = self.reserveObject('path %d' % len(self.paths))
        self.paths.append(
            (name, path, trans, ob, gc.get_joinstyle(), gc.get_capstyle(),
             padding, filled, stroked))
        return name

    def writePathCollectionTemplates(self):
        for (name, path, trans, ob, joinstyle, capstyle, padding, filled,
             stroked) in self.paths:
            pathops = self.pathOperations(path, trans, simplify=False)
            bbox = path.get_extents(trans)
            if not np.all(np.isfinite(bbox.extents)):
                extents = [0, 0, 0, 0]
            else:
                bbox = bbox.padded(padding)
                extents = list(bbox.extents)
            self.beginStream(
                ob.id, None,
                {'Type': Name('XObject'), 'Subtype': Name('Form'),
                 'BBox': extents})
            self.output(GraphicsContextPdf.joinstyles[joinstyle],
                        Op.setlinejoin)
            self.output(GraphicsContextPdf.capstyles[capstyle], Op.setlinecap)
            self.output(*pathops)
            self.output(Op.paint_path(filled, stroked))
            self.endStream()

    @staticmethod
    def pathOperations(path, transform, clip=None, simplify=None, sketch=None):
        return [Verbatim(_path.convert_to_string(
            path, transform, clip, simplify, sketch,
            6,
            [Op.moveto.op, Op.lineto.op, b'', Op.curveto.op, Op.closepath.op],
            True))]

    def writePath(self, path, transform, clip=False, sketch=None):
        if clip:
            clip = (0.0, 0.0, self.width * 72, self.height * 72)
            simplify = path.should_simplify
        else:
            clip = None
            simplify = False
        cmds = self.pathOperations(path, transform, clip, simplify=simplify,
                                   sketch=sketch)
        self.output(*cmds)

    def reserveObject(self, name=''):
        """Reserve an ID for an indirect object.
        The name is used for debugging in case we forget to print out
        the object with writeObject.
        """

        id = self.nextObject
        self.nextObject += 1
        self.xrefTable.append([None, 0, name])
        return Reference(id)

    def recordXref(self, id):
        self.xrefTable[id][0] = self.fh.tell() - self.tell_base

    def writeObject(self, object, contents):
        self.recordXref(object.id)
        object.write(contents, self)

    def writeXref(self):
        """Write out the xref table."""

        self.startxref = self.fh.tell() - self.tell_base
        self.write(("xref\n0 %d\n" % self.nextObject).encode('ascii'))
        i = 0
        borken = False
        for offset, generation, name in self.xrefTable:
            if offset is None:
                print('No offset for object %d (%s)' % (i, name),
                      file=sys.stderr)
                borken = True
            else:
                if name == 'the zero object':
                    key = "f"
                else:
                    key = "n"
                text = "%010d %05d %s \n" % (offset, generation, key)
                self.write(text.encode('ascii'))
            i += 1
        if borken:
            raise AssertionError('Indirect object does not exist')

    def writeInfoDict(self):
        """Write out the info dictionary, checking it for good form"""

        def is_string_like(x):
            return isinstance(x, six.string_types)

        def is_date(x):
            return isinstance(x, datetime)

        check_trapped = (lambda x: isinstance(x, Name) and
                         x.name in ('True', 'False', 'Unknown'))

        keywords = {'Title': is_string_like,
                    'Author': is_string_like,
                    'Subject': is_string_like,
                    'Keywords': is_string_like,
                    'Creator': is_string_like,
                    'Producer': is_string_like,
                    'CreationDate': is_date,
                    'ModDate': is_date,
                    'Trapped': check_trapped}
        for k in self.infoDict:
            if k not in keywords:
                warnings.warn('Unknown infodict keyword: %s' % k)
            else:
                if not keywords[k](self.infoDict[k]):
                    warnings.warn('Bad value for infodict keyword %s' % k)

        self.infoObject = self.reserveObject('info')
        self.writeObject(self.infoObject, self.infoDict)

    def writeTrailer(self):
        """Write out the PDF trailer."""

        self.write(b"trailer\n")
        self.write(pdfRepr(
            {'Size': self.nextObject,
             'Root': self.rootObject,
             'Info': self.infoObject}))
        # Could add 'ID'
        self.write(("\nstartxref\n%d\n%%%%EOF\n" %
                    self.startxref).encode('ascii'))


class RendererPdf(RendererBase):
    afm_font_cache = maxdict(50)

    def __init__(self, file, image_dpi, height, width):
        RendererBase.__init__(self)
        self.height = height
        self.width = width
        self.file = file
        self.gc = self.new_gc()
        self.mathtext_parser = MathTextParser("Pdf")
        self.image_dpi = image_dpi

    def finalize(self):
        self.file.output(*self.gc.finalize())

    def check_gc(self, gc, fillcolor=None):
        orig_fill = getattr(gc, '_fillcolor', (0., 0., 0.))
        gc._fillcolor = fillcolor

        orig_alphas = getattr(gc, '_effective_alphas', (1.0, 1.0))

        if gc.get_rgb() is None:
            # it should not matter what color here
            # since linewidth should be 0
            # unless affected by global settings in rcParams
            # hence setting zero alpha just incase
            gc.set_foreground((0, 0, 0, 0), isRGBA=True)

        if gc._forced_alpha:
            gc._effective_alphas = (gc._alpha, gc._alpha)
        elif fillcolor is None or len(fillcolor) < 4:
            gc._effective_alphas = (gc._rgb[3], 1.0)
        else:
            gc._effective_alphas = (gc._rgb[3], fillcolor[3])

        delta = self.gc.delta(gc)
        if delta:
            self.file.output(*delta)

        # Restore gc to avoid unwanted side effects
        gc._fillcolor = orig_fill
        gc._effective_alphas = orig_alphas

    def track_characters(self, font, s):
        """Keeps track of which characters are required from
        each font."""
        if isinstance(font, six.string_types):
            fname = font
        else:
            fname = font.fname
        realpath, stat_key = get_realpath_and_stat(fname)
        used_characters = self.file.used_characters.setdefault(
            stat_key, (realpath, set()))
        used_characters[1].update([ord(x) for x in s])

    def merge_used_characters(self, other):
        for stat_key, (realpath, charset) in six.iteritems(other):
            used_characters = self.file.used_characters.setdefault(
                stat_key, (realpath, set()))
            used_characters[1].update(charset)

    def get_image_magnification(self):
        return self.image_dpi/72.0

    def option_scale_image(self):
        """
        pdf backend support arbitrary scaling of image.
        """
        return True

    def option_image_nocomposite(self):
        """
        return whether to generate a composite image from multiple images on
        a set of axes
        """
        return not rcParams['image.composite_image']

    def draw_image(self, gc, x, y, im, transform=None):
        h, w = im.shape[:2]
        if w == 0 or h == 0:
            return

        if transform is None:
            # If there's no transform, alpha has already been applied
            gc.set_alpha(1.0)

        self.check_gc(gc)

        w = 72.0 * w / self.image_dpi
        h = 72.0 * h / self.image_dpi

        imob = self.file.imageObject(im)

        if transform is None:
            self.file.output(Op.gsave,
                             w, 0, 0, h, x, y, Op.concat_matrix,
                             imob, Op.use_xobject, Op.grestore)
        else:
            tr1, tr2, tr3, tr4, tr5, tr6 = transform.frozen().to_values()

            self.file.output(Op.gsave,
                             1, 0, 0, 1, x, y, Op.concat_matrix,
                             tr1, tr2, tr3, tr4, tr5, tr6, Op.concat_matrix,
                             imob, Op.use_xobject, Op.grestore)

    def draw_path(self, gc, path, transform, rgbFace=None):
        self.check_gc(gc, rgbFace)
        self.file.writePath(
            path, transform,
            rgbFace is None and gc.get_hatch_path() is None,
            gc.get_sketch_params())
        self.file.output(self.gc.paint())

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offsetTrans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position):
        # We can only reuse the objects if the presence of fill and
        # stroke (and the amount of alpha for each) is the same for
        # all of them
        can_do_optimization = True
        facecolors = np.asarray(facecolors)
        edgecolors = np.asarray(edgecolors)

        if not len(facecolors):
            filled = False
            can_do_optimization = not gc.get_hatch()
        else:
            if np.all(facecolors[:, 3] == facecolors[0, 3]):
                filled = facecolors[0, 3] != 0.0
            else:
                can_do_optimization = False

        if not len(edgecolors):
            stroked = False
        else:
            if np.all(np.asarray(linewidths) == 0.0):
                stroked = False
            elif np.all(edgecolors[:, 3] == edgecolors[0, 3]):
                stroked = edgecolors[0, 3] != 0.0
            else:
                can_do_optimization = False

        # Is the optimization worth it? Rough calculation:
        # cost of emitting a path in-line is len_path * uses_per_path
        # cost of XObject is len_path + 5 for the definition,
        #    uses_per_path for the uses
        len_path = len(paths[0].vertices) if len(paths) > 0 else 0
        uses_per_path = self._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        should_do_optimization = \
            len_path + uses_per_path + 5 < len_path * uses_per_path

        if (not can_do_optimization) or (not should_do_optimization):
            return RendererBase.draw_path_collection(
                self, gc, master_transform, paths, all_transforms,
                offsets, offsetTrans, facecolors, edgecolors,
                linewidths, linestyles, antialiaseds, urls,
                offset_position)

        padding = np.max(linewidths)
        path_codes = []
        for i, (path, transform) in enumerate(self._iter_collection_raw_paths(
                master_transform, paths, all_transforms)):
            name = self.file.pathCollectionObject(
                gc, path, transform, padding, filled, stroked)
            path_codes.append(name)

        output = self.file.output
        output(*self.gc.push())
        lastx, lasty = 0, 0
        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
                gc, master_transform, all_transforms, path_codes, offsets,
                offsetTrans, facecolors, edgecolors, linewidths, linestyles,
                antialiaseds, urls, offset_position):

            self.check_gc(gc0, rgbFace)
            dx, dy = xo - lastx, yo - lasty
            output(1, 0, 0, 1, dx, dy, Op.concat_matrix, path_id,
                   Op.use_xobject)
            lastx, lasty = xo, yo
        output(*self.gc.pop())

    def draw_markers(self, gc, marker_path, marker_trans, path, trans,
                     rgbFace=None):
        # Same logic as in draw_path_collection
        len_marker_path = len(marker_path)
        uses = len(path)
        if len_marker_path * uses < len_marker_path + uses + 5:
            RendererBase.draw_markers(self, gc, marker_path, marker_trans,
                                      path, trans, rgbFace)
            return

        self.check_gc(gc, rgbFace)
        fill = gc.fill(rgbFace)
        stroke = gc.stroke()

        output = self.file.output
        marker = self.file.markerObject(
            marker_path, marker_trans, fill, stroke, self.gc._linewidth,
            gc.get_joinstyle(), gc.get_capstyle())

        output(Op.gsave)
        lastx, lasty = 0, 0
        for vertices, code in path.iter_segments(
                trans,
                clip=(0, 0, self.file.width*72, self.file.height*72),
                simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                if (x < 0 or y < 0 or
                        x > self.file.width * 72 or y > self.file.height * 72):
                    continue
                dx, dy = x - lastx, y - lasty
                output(1, 0, 0, 1, dx, dy, Op.concat_matrix,
                       marker, Op.use_xobject)
                lastx, lasty = x, y
        output(Op.grestore)

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
        points = points.reshape((shape[0] * shape[1], 2))
        tpoints = trans.transform(points)
        tpoints = tpoints.reshape(shape)
        name = self.file.addGouraudTriangles(tpoints, colors)
        self.check_gc(gc)
        self.file.output(name, Op.shading)

    def _setup_textpos(self, x, y, angle, oldx=0, oldy=0, oldangle=0):
        if angle == oldangle == 0:
            self.file.output(x - oldx, y - oldy, Op.textpos)
        else:
            angle = angle / 180.0 * pi
            self.file.output(cos(angle), sin(angle),
                             -sin(angle), cos(angle),
                             x, y, Op.textmatrix)
            self.file.output(0, 0, Op.textpos)

    def draw_mathtext(self, gc, x, y, s, prop, angle):
        # TODO: fix positioning and encoding
        width, height, descent, glyphs, rects, used_characters = \
            self.mathtext_parser.parse(s, 72, prop)
        self.merge_used_characters(used_characters)

        # When using Type 3 fonts, we can't use character codes higher
        # than 255, so we use the "Do" command to render those
        # instead.
        global_fonttype = rcParams['pdf.fonttype']

        # Set up a global transformation matrix for the whole math expression
        a = angle / 180.0 * pi
        self.file.output(Op.gsave)
        self.file.output(cos(a), sin(a), -sin(a), cos(a), x, y,
                         Op.concat_matrix)

        self.check_gc(gc, gc._rgb)
        self.file.output(Op.begin_text)
        prev_font = None, None
        oldx, oldy = 0, 0
        for ox, oy, fontname, fontsize, num, symbol_name in glyphs:
            if is_opentype_cff_font(fontname):
                fonttype = 42
            else:
                fonttype = global_fonttype

            if fonttype == 42 or num <= 255:
                self._setup_textpos(ox, oy, 0, oldx, oldy)
                oldx, oldy = ox, oy
                if (fontname, fontsize) != prev_font:
                    self.file.output(self.file.fontName(fontname), fontsize,
                                     Op.selectfont)
                    prev_font = fontname, fontsize
                self.file.output(self.encode_string(unichr(num), fonttype),
                                 Op.show)
        self.file.output(Op.end_text)

        # If using Type 3 fonts, render all of the multi-byte characters
        # as XObjects using the 'Do' command.
        if global_fonttype == 3:
            for ox, oy, fontname, fontsize, num, symbol_name in glyphs:
                if is_opentype_cff_font(fontname):
                    fonttype = 42
                else:
                    fonttype = global_fonttype

                if fonttype == 3 and num > 255:
                    self.file.fontName(fontname)
                    self.file.output(Op.gsave,
                                     0.001 * fontsize, 0,
                                     0, 0.001 * fontsize,
                                     ox, oy, Op.concat_matrix)
                    name = self.file._get_xobject_symbol_name(
                        fontname, symbol_name)
                    self.file.output(Name(name), Op.use_xobject)
                    self.file.output(Op.grestore)

        # Draw any horizontal lines in the math layout
        for ox, oy, width, height in rects:
            self.file.output(Op.gsave, ox, oy, width, height,
                             Op.rectangle, Op.fill, Op.grestore)

        # Pop off the global transformation
        self.file.output(Op.grestore)

    def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
        texmanager = self.get_texmanager()
        fontsize = prop.get_size_in_points()
        dvifile = texmanager.make_dvi(s, fontsize)
        with dviread.Dvi(dvifile, 72) as dvi:
            page = next(iter(dvi))

        # Gather font information and do some setup for combining
        # characters into strings. The variable seq will contain a
        # sequence of font and text entries. A font entry is a list
        # ['font', name, size] where name is a Name object for the
        # font. A text entry is ['text', x, y, glyphs, x+w] where x
        # and y are the starting coordinates, w is the width, and
        # glyphs is a list; in this phase it will always contain just
        # one one-character string, but later it may have longer
        # strings interspersed with kern amounts.
        oldfont, seq = None, []
        for x1, y1, dvifont, glyph, width in page.text:
            if dvifont != oldfont:
                pdfname = self.file.dviFontName(dvifont)
                seq += [['font', pdfname, dvifont.size]]
                oldfont = dvifont
            # We need to convert the glyph numbers to bytes, and the easiest
            # way to do this on both Python 2 and 3 is .encode('latin-1')
            seq += [['text', x1, y1,
                     [six.unichr(glyph).encode('latin-1')], x1+width]]

        # Find consecutive text strings with constant y coordinate and
        # combine into a sequence of strings and kerns, or just one
        # string (if any kerns would be less than 0.1 points).
        i, curx, fontsize = 0, 0, None
        while i < len(seq)-1:
            elt, nxt = seq[i:i+2]
            if elt[0] == 'font':
                fontsize = elt[2]
            elif elt[0] == nxt[0] == 'text' and elt[2] == nxt[2]:
                offset = elt[4] - nxt[1]
                if abs(offset) < 0.1:
                    elt[3][-1] += nxt[3][0]
                    elt[4] += nxt[4]-nxt[1]
                else:
                    elt[3] += [offset*1000.0/fontsize, nxt[3][0]]
                    elt[4] = nxt[4]
                del seq[i+1]
                continue
            i += 1

        # Create a transform to map the dvi contents to the canvas.
        mytrans = Affine2D().rotate_deg(angle).translate(x, y)

        # Output the text.
        self.check_gc(gc, gc._rgb)
        self.file.output(Op.begin_text)
        curx, cury, oldx, oldy = 0, 0, 0, 0
        for elt in seq:
            if elt[0] == 'font':
                self.file.output(elt[1], elt[2], Op.selectfont)
            elif elt[0] == 'text':
                curx, cury = mytrans.transform_point((elt[1], elt[2]))
                self._setup_textpos(curx, cury, angle, oldx, oldy)
                oldx, oldy = curx, cury
                if len(elt[3]) == 1:
                    self.file.output(elt[3][0], Op.show)
                else:
                    self.file.output(elt[3], Op.showkern)
            else:
                assert False
        self.file.output(Op.end_text)

        # Then output the boxes (e.g., variable-length lines of square
        # roots).
        boxgc = self.new_gc()
        boxgc.copy_properties(gc)
        boxgc.set_linewidth(0)
        pathops = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
                   Path.CLOSEPOLY]
        for x1, y1, h, w in page.boxes:
            path = Path([[x1, y1], [x1+w, y1], [x1+w, y1+h], [x1, y1+h],
                         [0, 0]], pathops)
            self.draw_path(boxgc, path, mytrans, gc._rgb)

    def encode_string(self, s, fonttype):
        if fonttype in (1, 3):
            return s.encode('cp1252', 'replace')
        return s.encode('utf-16be', 'replace')

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        # TODO: combine consecutive texts into one BT/ET delimited section

        # This function is rather complex, since there is no way to
        # access characters of a Type 3 font with codes > 255.  (Type
        # 3 fonts can not have a CIDMap).  Therefore, we break the
        # string into chunks, where each chunk contains exclusively
        # 1-byte or exclusively 2-byte characters, and output each
        # chunk a separate command.  1-byte characters use the regular
        # text show command (Tj), whereas 2-byte characters use the
        # use XObject command (Do).  If using Type 42 fonts, all of
        # this complication is avoided, but of course, those fonts can
        # not be subsetted.

        self.check_gc(gc, gc._rgb)
        if ismath:
            return self.draw_mathtext(gc, x, y, s, prop, angle)

        fontsize = prop.get_size_in_points()

        if rcParams['pdf.use14corefonts']:
            font = self._get_font_afm(prop)
            l, b, w, h = font.get_str_bbox(s)
            fonttype = 1
        else:
            font = self._get_font_ttf(prop)
            self.track_characters(font, s)
            font.set_text(s, 0.0, flags=LOAD_NO_HINTING)

            fonttype = rcParams['pdf.fonttype']

            # We can't subset all OpenType fonts, so switch to Type 42
            # in that case.
            if is_opentype_cff_font(font.fname):
                fonttype = 42

        def check_simple_method(s):
            """Determine if we should use the simple or woven method
            to output this text, and chunks the string into 1-byte and
            2-byte sections if necessary."""
            use_simple_method = True
            chunks = []

            if not rcParams['pdf.use14corefonts']:
                if fonttype == 3 and not isinstance(s, bytes) and len(s) != 0:
                    # Break the string into chunks where each chunk is either
                    # a string of chars <= 255, or a single character > 255.
                    s = six.text_type(s)
                    for c in s:
                        if ord(c) <= 255:
                            char_type = 1
                        else:
                            char_type = 2
                        if len(chunks) and chunks[-1][0] == char_type:
                            chunks[-1][1].append(c)
                        else:
                            chunks.append((char_type, [c]))
                    use_simple_method = (len(chunks) == 1 and
                                         chunks[-1][0] == 1)
            return use_simple_method, chunks

        def draw_text_simple():
            """Outputs text using the simple method."""
            self.file.output(Op.begin_text,
                             self.file.fontName(prop),
                             fontsize,
                             Op.selectfont)
            self._setup_textpos(x, y, angle)
            self.file.output(self.encode_string(s, fonttype), Op.show,
                             Op.end_text)

        def draw_text_woven(chunks):
            """Outputs text using the woven method, alternating
            between chunks of 1-byte characters and 2-byte characters.
            Only used for Type 3 fonts."""
            chunks = [(a, ''.join(b)) for a, b in chunks]

            # Do the rotation and global translation as a single matrix
            # concatenation up front
            self.file.output(Op.gsave)
            a = angle / 180.0 * pi
            self.file.output(cos(a), sin(a), -sin(a), cos(a), x, y,
                             Op.concat_matrix)

            # Output all the 1-byte characters in a BT/ET group, then
            # output all the 2-byte characters.
            for mode in (1, 2):
                newx = oldx = 0
                # Output a 1-byte character chunk
                if mode == 1:
                    self.file.output(Op.begin_text,
                                     self.file.fontName(prop),
                                     fontsize,
                                     Op.selectfont)

                for chunk_type, chunk in chunks:
                    if mode == 1 and chunk_type == 1:
                        self._setup_textpos(newx, 0, 0, oldx, 0, 0)
                        self.file.output(self.encode_string(chunk, fonttype),
                                         Op.show)
                        oldx = newx

                    lastgind = None
                    for c in chunk:
                        ccode = ord(c)
                        gind = font.get_char_index(ccode)
                        if gind is not None:
                            if mode == 2 and chunk_type == 2:
                                glyph_name = font.get_glyph_name(gind)
                                self.file.output(Op.gsave)
                                self.file.output(0.001 * fontsize, 0,
                                                 0, 0.001 * fontsize,
                                                 newx, 0, Op.concat_matrix)
                                name = self.file._get_xobject_symbol_name(
                                    font.fname, glyph_name)
                                self.file.output(Name(name), Op.use_xobject)
                                self.file.output(Op.grestore)

                            # Move the pointer based on the character width
                            # and kerning
                            glyph = font.load_char(ccode,
                                                   flags=LOAD_NO_HINTING)
                            if lastgind is not None:
                                kern = font.get_kerning(
                                    lastgind, gind, KERNING_UNFITTED)
                            else:
                                kern = 0
                            lastgind = gind
                            newx += kern/64.0 + glyph.linearHoriAdvance/65536.0

                if mode == 1:
                    self.file.output(Op.end_text)

            self.file.output(Op.grestore)

        use_simple_method, chunks = check_simple_method(s)
        if use_simple_method:
            return draw_text_simple()
        else:
            return draw_text_woven(chunks)

    def get_text_width_height_descent(self, s, prop, ismath):
        if rcParams['text.usetex']:
            texmanager = self.get_texmanager()
            fontsize = prop.get_size_in_points()
            w, h, d = texmanager.get_text_width_height_descent(s, fontsize,
                                                               renderer=self)
            return w, h, d

        if ismath:
            w, h, d, glyphs, rects, used_characters = \
                self.mathtext_parser.parse(s, 72, prop)

        elif rcParams['pdf.use14corefonts']:
            font = self._get_font_afm(prop)
            l, b, w, h, d = font.get_str_bbox_and_descent(s)
            scale = prop.get_size_in_points()
            w *= scale / 1000
            h *= scale / 1000
            d *= scale / 1000
        else:
            font = self._get_font_ttf(prop)
            font.set_text(s, 0.0, flags=LOAD_NO_HINTING)
            w, h = font.get_width_height()
            scale = (1.0 / 64.0)
            w *= scale
            h *= scale
            d = font.get_descent()
            d *= scale
        return w, h, d

    def _get_font_afm(self, prop):
        key = hash(prop)
        font = self.afm_font_cache.get(key)
        if font is None:
            filename = findfont(
                prop, fontext='afm', directory=self.file._core14fontdir)
            if filename is None:
                filename = findfont(
                    "Helvetica", fontext='afm',
                    directory=self.file._core14fontdir)
            font = self.afm_font_cache.get(filename)
            if font is None:
                with open(filename, 'rb') as fh:
                    font = AFM(fh)
                    self.afm_font_cache[filename] = font
            self.afm_font_cache[key] = font
        return font

    def _get_font_ttf(self, prop):
        filename = findfont(prop)
        font = get_font(filename)
        font.clear()
        font.set_size(prop.get_size_in_points(), 72)
        return font

    def flipy(self):
        return False

    def get_canvas_width_height(self):
        return self.file.width * 72.0, self.file.height * 72.0

    def new_gc(self):
        return GraphicsContextPdf(self.file)


class GraphicsContextPdf(GraphicsContextBase):

    def __init__(self, file):
        GraphicsContextBase.__init__(self)
        self._fillcolor = (0.0, 0.0, 0.0)
        self._effective_alphas = (1.0, 1.0)
        self.file = file
        self.parent = None

    def __repr__(self):
        d = dict(self.__dict__)
        del d['file']
        del d['parent']
        return repr(d)

    def stroke(self):
        """
        Predicate: does the path need to be stroked (its outline drawn)?
        This tests for the various conditions that disable stroking
        the path, in which case it would presumably be filled.
        """
        # _linewidth > 0: in pdf a line of width 0 is drawn at minimum
        #   possible device width, but e.g., agg doesn't draw at all
        return (self._linewidth > 0 and self._alpha > 0 and
                (len(self._rgb) <= 3 or self._rgb[3] != 0.0))

    def fill(self, *args):
        """
        Predicate: does the path need to be filled?

        An optional argument can be used to specify an alternative
        _fillcolor, as needed by RendererPdf.draw_markers.
        """
        if len(args):
            _fillcolor = args[0]
        else:
            _fillcolor = self._fillcolor
        return (self._hatch or
                (_fillcolor is not None and
                 (len(_fillcolor) <= 3 or _fillcolor[3] != 0.0)))

    def paint(self):
        """
        Return the appropriate pdf operator to cause the path to be
        stroked, filled, or both.
        """
        return Op.paint_path(self.fill(), self.stroke())

    capstyles = {'butt': 0, 'round': 1, 'projecting': 2}
    joinstyles = {'miter': 0, 'round': 1, 'bevel': 2}

    def capstyle_cmd(self, style):
        return [self.capstyles[style], Op.setlinecap]

    def joinstyle_cmd(self, style):
        return [self.joinstyles[style], Op.setlinejoin]

    def linewidth_cmd(self, width):
        return [width, Op.setlinewidth]

    def dash_cmd(self, dashes):
        offset, dash = dashes
        if dash is None:
            dash = []
            offset = 0
        return [list(dash), offset, Op.setdash]

    def alpha_cmd(self, alpha, forced, effective_alphas):
        name = self.file.alphaState(effective_alphas)
        return [name, Op.setgstate]

    def hatch_cmd(self, hatch, hatch_color):
        if not hatch:
            if self._fillcolor is not None:
                return self.fillcolor_cmd(self._fillcolor)
            else:
                return [Name('DeviceRGB'), Op.setcolorspace_nonstroke]
        else:
            hatch_style = (hatch_color, self._fillcolor, hatch)
            name = self.file.hatchPattern(hatch_style)
            return [Name('Pattern'), Op.setcolorspace_nonstroke,
                    name, Op.setcolor_nonstroke]

    def rgb_cmd(self, rgb):
        if rcParams['pdf.inheritcolor']:
            return []
        if rgb[0] == rgb[1] == rgb[2]:
            return [rgb[0], Op.setgray_stroke]
        else:
            return list(rgb[:3]) + [Op.setrgb_stroke]

    def fillcolor_cmd(self, rgb):
        if rgb is None or rcParams['pdf.inheritcolor']:
            return []
        elif rgb[0] == rgb[1] == rgb[2]:
            return [rgb[0], Op.setgray_nonstroke]
        else:
            return list(rgb[:3]) + [Op.setrgb_nonstroke]

    def push(self):
        parent = GraphicsContextPdf(self.file)
        parent.copy_properties(self)
        parent.parent = self.parent
        self.parent = parent
        return [Op.gsave]

    def pop(self):
        assert self.parent is not None
        self.copy_properties(self.parent)
        self.parent = self.parent.parent
        return [Op.grestore]

    def clip_cmd(self, cliprect, clippath):
        """Set clip rectangle. Calls self.pop() and self.push()."""
        cmds = []
        # Pop graphics state until we hit the right one or the stack is empty
        while ((self._cliprect, self._clippath) != (cliprect, clippath)
                and self.parent is not None):
            cmds.extend(self.pop())
        # Unless we hit the right one, set the clip polygon
        if ((self._cliprect, self._clippath) != (cliprect, clippath) or
                self.parent is None):
            cmds.extend(self.push())
            if self._cliprect != cliprect:
                cmds.extend([cliprect, Op.rectangle, Op.clip, Op.endpath])
            if self._clippath != clippath:
                path, affine = clippath.get_transformed_path_and_affine()
                cmds.extend(
                    PdfFile.pathOperations(path, affine, simplify=False) +
                    [Op.clip, Op.endpath])
        return cmds

    commands = (
        # must come first since may pop
        (('_cliprect', '_clippath'), clip_cmd),
        (('_alpha', '_forced_alpha', '_effective_alphas'), alpha_cmd),
        (('_capstyle',), capstyle_cmd),
        (('_fillcolor',), fillcolor_cmd),
        (('_joinstyle',), joinstyle_cmd),
        (('_linewidth',), linewidth_cmd),
        (('_dashes',), dash_cmd),
        (('_rgb',), rgb_cmd),
        # must come after fillcolor and rgb
        (('_hatch', '_hatch_color'), hatch_cmd),
        )

    def delta(self, other):
        """
        Copy properties of other into self and return PDF commands
        needed to transform self into other.
        """
        cmds = []
        fill_performed = False
        for params, cmd in self.commands:
            different = False
            for p in params:
                ours = getattr(self, p)
                theirs = getattr(other, p)
                try:
                    if ours is None or theirs is None:
                        different = ours is not theirs
                    else:
                        different = bool(ours != theirs)
                except ValueError:
                    ours = np.asarray(ours)
                    theirs = np.asarray(theirs)
                    different = (ours.shape != theirs.shape or
                                 np.any(ours != theirs))
                if different:
                    break

            # Need to update hatching if we also updated fillcolor
            if params == ('_hatch', '_hatch_color') and fill_performed:
                different = True

            if different:
                if params == ('_fillcolor',):
                    fill_performed = True
                theirs = [getattr(other, p) for p in params]
                cmds.extend(cmd(self, *theirs))
                for p in params:
                    setattr(self, p, getattr(other, p))
        return cmds

    def copy_properties(self, other):
        """
        Copy properties of other into self.
        """
        GraphicsContextBase.copy_properties(self, other)
        fillcolor = getattr(other, '_fillcolor', self._fillcolor)
        effective_alphas = getattr(other, '_effective_alphas',
                                   self._effective_alphas)
        self._fillcolor = fillcolor
        self._effective_alphas = effective_alphas

    def finalize(self):
        """
        Make sure every pushed graphics state is popped.
        """
        cmds = []
        while self.parent is not None:
            cmds.extend(self.pop())
        return cmds

########################################################################
#
# The following functions and classes are for pylab and implement
# window/figure managers, etc...
#
########################################################################


class PdfPages(object):
    """
    A multi-page PDF file.

    Examples
    --------

    >>> import matplotlib.pyplot as plt
    >>> # Initialize:
    >>> with PdfPages('foo.pdf') as pdf:
    ...     # As many times as you like, create a figure fig and save it:
    ...     fig = plt.figure()
    ...     pdf.savefig(fig)
    ...     # When no figure is specified the current figure is saved
    ...     pdf.savefig()

    Notes
    -----

    In reality :class:`PdfPages` is a thin wrapper around :class:`PdfFile`, in
    order to avoid confusion when using :func:`~matplotlib.pyplot.savefig` and
    forgetting the format argument.
    """
    __slots__ = ('_file', 'keep_empty')

    def __init__(self, filename, keep_empty=True, metadata=None):
        """
        Create a new PdfPages object.

        Parameters
        ----------

        filename : str
            Plots using :meth:`PdfPages.savefig` will be written to a file at
            this location. The file is opened at once and any older file with
            the same name is overwritten.
        keep_empty : bool, optional
            If set to False, then empty pdf files will be deleted automatically
            when closed.
        metadata : dictionary, optional
            Information dictionary object (see PDF reference section 10.2.1
            'Document Information Dictionary'), e.g.:
            `{'Creator': 'My software', 'Author': 'Me',
            'Title': 'Awesome fig'}`

            The standard keys are `'Title'`, `'Author'`, `'Subject'`,
            `'Keywords'`, `'Creator'`, `'Producer'`, `'CreationDate'`,
            `'ModDate'`, and `'Trapped'`. Values have been predefined
            for `'Creator'`, `'Producer'` and `'CreationDate'`. They
            can be removed by setting them to `None`.

        """
        self._file = PdfFile(filename, metadata=metadata)
        self.keep_empty = keep_empty

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """
        Finalize this object, making the underlying file a complete
        PDF file.
        """
        self._file.finalize()
        self._file.close()
        if (self.get_pagecount() == 0 and not self.keep_empty and
                not self._file.passed_in_file_object):
            os.remove(self._file.fh.name)
        self._file = None

    def infodict(self):
        """
        Return a modifiable information dictionary object
        (see PDF reference section 10.2.1 'Document Information
        Dictionary').
        """
        return self._file.infoDict

    def savefig(self, figure=None, **kwargs):
        """
        Saves a :class:`~matplotlib.figure.Figure` to this file as a new page.

        Any other keyword arguments are passed to
        :meth:`~matplotlib.figure.Figure.savefig`.

        Parameters
        ----------

        figure : :class:`~matplotlib.figure.Figure` or int, optional
            Specifies what figure is saved to file. If not specified, the
            active figure is saved. If a :class:`~matplotlib.figure.Figure`
            instance is provided, this figure is saved. If an int is specified,
            the figure instance to save is looked up by number.
        """
        if not isinstance(figure, Figure):
            if figure is None:
                manager = Gcf.get_active()
            else:
                manager = Gcf.get_fig_manager(figure)
            if manager is None:
                raise ValueError("No figure {}".format(figure))
            figure = manager.canvas.figure
        # Force use of pdf backend, as PdfPages is tightly coupled with it.
        try:
            orig_canvas = figure.canvas
            figure.canvas = FigureCanvasPdf(figure)
            figure.savefig(self, format="pdf", **kwargs)
        finally:
            figure.canvas = orig_canvas

    def get_pagecount(self):
        """
        Returns the current number of pages in the multipage pdf file.
        """
        return len(self._file.pageList)

    def attach_note(self, text, positionRect=[-100, -100, 0, 0]):
        """
        Add a new text note to the page to be saved next. The optional
        positionRect specifies the position of the new note on the
        page. It is outside the page per default to make sure it is
        invisible on printouts.
        """
        self._file.newTextnote(text, positionRect)


class FigureCanvasPdf(FigureCanvasBase):
    """
    The canvas the figure renders into.  Calls the draw and print fig
    methods, creates the renderers, etc...

    Attributes
    ----------
    figure : `matplotlib.figure.Figure`
        A high-level Figure instance

    """

    fixed_dpi = 72

    def draw(self):
        pass

    filetypes = {'pdf': 'Portable Document Format'}

    def get_default_filetype(self):
        return 'pdf'

    def print_pdf(self, filename, **kwargs):
        image_dpi = kwargs.get('dpi', 72)  # dpi to use for images
        self.figure.set_dpi(72)            # there are 72 pdf points to an inch
        width, height = self.figure.get_size_inches()
        if isinstance(filename, PdfPages):
            file = filename._file
        else:
            file = PdfFile(filename, metadata=kwargs.pop("metadata", None))
        try:
            file.newPage(width, height)
            _bbox_inches_restore = kwargs.pop("bbox_inches_restore", None)
            renderer = MixedModeRenderer(
                self.figure, width, height, image_dpi,
                RendererPdf(file, image_dpi, height, width),
                bbox_inches_restore=_bbox_inches_restore)
            self.figure.draw(renderer)
            renderer.finalize()
            if not isinstance(filename, PdfPages):
                file.finalize()
        finally:
            if isinstance(filename, PdfPages):  # finish off this page
                file.endStream()
            else:            # we opened the file above; now finish it off
                file.close()


class FigureManagerPdf(FigureManagerBase):
    pass


@_Backend.export
class _BackendPdf(_Backend):
    FigureCanvas = FigureCanvasPdf
    FigureManager = FigureManagerPdf
