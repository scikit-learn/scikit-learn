"""
A module for reading dvi files output by TeX. Several limitations make
this not (currently) useful as a general-purpose dvi preprocessor, but
it is currently used by the pdf backend for processing usetex text.

Interface::

  with Dvi(filename, 72) as dvi:
      # iterate over pages:
      for page in dvi:
          w, h, d = page.width, page.height, page.descent
          for x,y,font,glyph,width in page.text:
              fontname = font.texname
              pointsize = font.size
              ...
          for x,y,height,width in page.boxes:
              ...

"""
from __future__ import absolute_import, division, print_function

import six
from six.moves import xrange

from collections import namedtuple
from functools import partial, wraps
import logging
import numpy as np
import os
import re
import struct
import sys
import textwrap

from matplotlib import cbook, rcParams
from matplotlib.compat import subprocess

try:
    from functools import lru_cache
except ImportError:  # Py2
    from backports.functools_lru_cache import lru_cache

if six.PY3:
    def ord(x):
        return x

_log = logging.getLogger(__name__)

# Dvi is a bytecode format documented in
# http://mirrors.ctan.org/systems/knuth/dist/texware/dvitype.web
# http://texdoc.net/texmf-dist/doc/generic/knuth/texware/dvitype.pdf
#
# The file consists of a preamble, some number of pages, a postamble,
# and a finale. Different opcodes are allowed in different contexts,
# so the Dvi object has a parser state:
#
#   pre:       expecting the preamble
#   outer:     between pages (followed by a page or the postamble,
#              also e.g. font definitions are allowed)
#   page:      processing a page
#   post_post: state after the postamble (our current implementation
#              just stops reading)
#   finale:    the finale (unimplemented in our current implementation)

_dvistate = cbook.Bunch(pre=0, outer=1, inpage=2, post_post=3, finale=4)

# The marks on a page consist of text and boxes. A page also has dimensions.
Page = namedtuple('Page', 'text boxes height width descent')
Text = namedtuple('Text', 'x y font glyph width')
Box = namedtuple('Box', 'x y height width')


# Opcode argument parsing
#
# Each of the following functions takes a Dvi object and delta,
# which is the difference between the opcode and the minimum opcode
# with the same meaning. Dvi opcodes often encode the number of
# argument bytes in this delta.

def _arg_raw(dvi, delta):
    """Return *delta* without reading anything more from the dvi file"""
    return delta


def _arg(bytes, signed, dvi, _):
    """Read *bytes* bytes, returning the bytes interpreted as a
    signed integer if *signed* is true, unsigned otherwise."""
    return dvi._arg(bytes, signed)


def _arg_slen(dvi, delta):
    """Signed, length *delta*

    Read *delta* bytes, returning None if *delta* is zero, and
    the bytes interpreted as a signed integer otherwise."""
    if delta == 0:
        return None
    return dvi._arg(delta, True)


def _arg_slen1(dvi, delta):
    """Signed, length *delta*+1

    Read *delta*+1 bytes, returning the bytes interpreted as signed."""
    return dvi._arg(delta+1, True)


def _arg_ulen1(dvi, delta):
    """Unsigned length *delta*+1

    Read *delta*+1 bytes, returning the bytes interpreted as unsigned."""
    return dvi._arg(delta+1, False)


def _arg_olen1(dvi, delta):
    """Optionally signed, length *delta*+1

    Read *delta*+1 bytes, returning the bytes interpreted as
    unsigned integer for 0<=*delta*<3 and signed if *delta*==3."""
    return dvi._arg(delta + 1, delta == 3)


_arg_mapping = dict(raw=_arg_raw,
                    u1=partial(_arg, 1, False),
                    u4=partial(_arg, 4, False),
                    s4=partial(_arg, 4, True),
                    slen=_arg_slen,
                    olen1=_arg_olen1,
                    slen1=_arg_slen1,
                    ulen1=_arg_ulen1)


def _dispatch(table, min, max=None, state=None, args=('raw',)):
    """Decorator for dispatch by opcode. Sets the values in *table*
    from *min* to *max* to this method, adds a check that the Dvi state
    matches *state* if not None, reads arguments from the file according
    to *args*.

    *table*
        the dispatch table to be filled in

    *min*
        minimum opcode for calling this function

    *max*
        maximum opcode for calling this function, None if only *min* is allowed

    *state*
        state of the Dvi object in which these opcodes are allowed

    *args*
        sequence of argument specifications:

        ``'raw'``: opcode minus minimum
        ``'u1'``: read one unsigned byte
        ``'u4'``: read four bytes, treat as an unsigned number
        ``'s4'``: read four bytes, treat as a signed number
        ``'slen'``: read (opcode - minimum) bytes, treat as signed
        ``'slen1'``: read (opcode - minimum + 1) bytes, treat as signed
        ``'ulen1'``: read (opcode - minimum + 1) bytes, treat as unsigned
        ``'olen1'``: read (opcode - minimum + 1) bytes, treat as unsigned
                     if under four bytes, signed if four bytes
    """
    def decorate(method):
        get_args = [_arg_mapping[x] for x in args]

        @wraps(method)
        def wrapper(self, byte):
            if state is not None and self.state != state:
                raise ValueError("state precondition failed")
            return method(self, *[f(self, byte-min) for f in get_args])
        if max is None:
            table[min] = wrapper
        else:
            for i in xrange(min, max+1):
                assert table[i] is None
                table[i] = wrapper
        return wrapper
    return decorate


class Dvi(object):
    """
    A reader for a dvi ("device-independent") file, as produced by TeX.
    The current implementation can only iterate through pages in order,
    and does not even attempt to verify the postamble.

    This class can be used as a context manager to close the underlying
    file upon exit. Pages can be read via iteration. Here is an overly
    simple way to extract text without trying to detect whitespace::

    >>> with matplotlib.dviread.Dvi('input.dvi', 72) as dvi:
    >>>     for page in dvi:
    >>>         print(''.join(unichr(t.glyph) for t in page.text))
    """
    # dispatch table
    _dtable = [None for _ in xrange(256)]
    _dispatch = partial(_dispatch, _dtable)

    def __init__(self, filename, dpi):
        """
        Read the data from the file named *filename* and convert
        TeX's internal units to units of *dpi* per inch.
        *dpi* only sets the units and does not limit the resolution.
        Use None to return TeX's internal units.
        """
        _log.debug('Dvi: %s', filename)
        self.file = open(filename, 'rb')
        self.dpi = dpi
        self.fonts = {}
        self.state = _dvistate.pre
        self.baseline = self._get_baseline(filename)

    def _get_baseline(self, filename):
        if rcParams['text.latex.preview']:
            base, ext = os.path.splitext(filename)
            baseline_filename = base + ".baseline"
            if os.path.exists(baseline_filename):
                with open(baseline_filename, 'rb') as fd:
                    l = fd.read().split()
                height, depth, width = l
                return float(depth)
        return None

    def __enter__(self):
        """
        Context manager enter method, does nothing.
        """
        return self

    def __exit__(self, etype, evalue, etrace):
        """
        Context manager exit method, closes the underlying file if it is open.
        """
        self.close()

    def __iter__(self):
        """
        Iterate through the pages of the file.

        Yields
        ------
        Page
            Details of all the text and box objects on the page.
            The Page tuple contains lists of Text and Box tuples and
            the page dimensions, and the Text and Box tuples contain
            coordinates transformed into a standard Cartesian
            coordinate system at the dpi value given when initializing.
            The coordinates are floating point numbers, but otherwise
            precision is not lost and coordinate values are not clipped to
            integers.
        """
        while True:
            have_page = self._read()
            if have_page:
                yield self._output()
            else:
                break

    def close(self):
        """
        Close the underlying file if it is open.
        """
        if not self.file.closed:
            self.file.close()

    def _output(self):
        """
        Output the text and boxes belonging to the most recent page.
        page = dvi._output()
        """
        minx, miny, maxx, maxy = np.inf, np.inf, -np.inf, -np.inf
        maxy_pure = -np.inf
        for elt in self.text + self.boxes:
            if isinstance(elt, Box):
                x, y, h, w = elt
                e = 0           # zero depth
            else:               # glyph
                x, y, font, g, w = elt
                h, e = font._height_depth_of(g)
            minx = min(minx, x)
            miny = min(miny, y - h)
            maxx = max(maxx, x + w)
            maxy = max(maxy, y + e)
            maxy_pure = max(maxy_pure, y)

        if self.dpi is None:
            # special case for ease of debugging: output raw dvi coordinates
            return Page(text=self.text, boxes=self.boxes,
                        width=maxx-minx, height=maxy_pure-miny,
                        descent=maxy-maxy_pure)

        # convert from TeX's "scaled points" to dpi units
        d = self.dpi / (72.27 * 2**16)
        if self.baseline is None:
            descent = (maxy - maxy_pure) * d
        else:
            descent = self.baseline

        text = [Text((x-minx)*d, (maxy-y)*d - descent, f, g, w*d)
                for (x, y, f, g, w) in self.text]
        boxes = [Box((x-minx)*d, (maxy-y)*d - descent, h*d, w*d)
                 for (x, y, h, w) in self.boxes]

        return Page(text=text, boxes=boxes, width=(maxx-minx)*d,
                    height=(maxy_pure-miny)*d, descent=descent)

    def _read(self):
        """
        Read one page from the file. Return True if successful,
        False if there were no more pages.
        """
        while True:
            byte = ord(self.file.read(1)[0])
            self._dtable[byte](self, byte)
            if byte == 140:                         # end of page
                return True
            if self.state == _dvistate.post_post:   # end of file
                self.close()
                return False

    def _arg(self, nbytes, signed=False):
        """
        Read and return an integer argument *nbytes* long.
        Signedness is determined by the *signed* keyword.
        """
        str = self.file.read(nbytes)
        value = ord(str[0])
        if signed and value >= 0x80:
            value = value - 0x100
        for i in range(1, nbytes):
            value = 0x100*value + ord(str[i])
        return value

    @_dispatch(min=0, max=127, state=_dvistate.inpage)
    def _set_char_immediate(self, char):
        self._put_char_real(char)
        self.h += self.fonts[self.f]._width_of(char)

    @_dispatch(min=128, max=131, state=_dvistate.inpage, args=('olen1',))
    def _set_char(self, char):
        self._put_char_real(char)
        self.h += self.fonts[self.f]._width_of(char)

    @_dispatch(132, state=_dvistate.inpage, args=('s4', 's4'))
    def _set_rule(self, a, b):
        self._put_rule_real(a, b)
        self.h += b

    @_dispatch(min=133, max=136, state=_dvistate.inpage, args=('olen1',))
    def _put_char(self, char):
        self._put_char_real(char)

    def _put_char_real(self, char):
        font = self.fonts[self.f]
        if font._vf is None:
            self.text.append(Text(self.h, self.v, font, char,
                                  font._width_of(char)))
        else:
            scale = font._scale
            for x, y, f, g, w in font._vf[char].text:
                newf = DviFont(scale=_mul2012(scale, f._scale),
                               tfm=f._tfm, texname=f.texname, vf=f._vf)
                self.text.append(Text(self.h + _mul2012(x, scale),
                                      self.v + _mul2012(y, scale),
                                      newf, g, newf._width_of(g)))
            self.boxes.extend([Box(self.h + _mul2012(x, scale),
                                   self.v + _mul2012(y, scale),
                                   _mul2012(a, scale), _mul2012(b, scale))
                               for x, y, a, b in font._vf[char].boxes])

    @_dispatch(137, state=_dvistate.inpage, args=('s4', 's4'))
    def _put_rule(self, a, b):
        self._put_rule_real(a, b)

    def _put_rule_real(self, a, b):
        if a > 0 and b > 0:
            self.boxes.append(Box(self.h, self.v, a, b))

    @_dispatch(138)
    def _nop(self, _):
        pass

    @_dispatch(139, state=_dvistate.outer, args=('s4',)*11)
    def _bop(self, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, p):
        self.state = _dvistate.inpage
        self.h, self.v, self.w, self.x, self.y, self.z = 0, 0, 0, 0, 0, 0
        self.stack = []
        self.text = []          # list of Text objects
        self.boxes = []         # list of Box objects

    @_dispatch(140, state=_dvistate.inpage)
    def _eop(self, _):
        self.state = _dvistate.outer
        del self.h, self.v, self.w, self.x, self.y, self.z, self.stack

    @_dispatch(141, state=_dvistate.inpage)
    def _push(self, _):
        self.stack.append((self.h, self.v, self.w, self.x, self.y, self.z))

    @_dispatch(142, state=_dvistate.inpage)
    def _pop(self, _):
        self.h, self.v, self.w, self.x, self.y, self.z = self.stack.pop()

    @_dispatch(min=143, max=146, state=_dvistate.inpage, args=('slen1',))
    def _right(self, b):
        self.h += b

    @_dispatch(min=147, max=151, state=_dvistate.inpage, args=('slen',))
    def _right_w(self, new_w):
        if new_w is not None:
            self.w = new_w
        self.h += self.w

    @_dispatch(min=152, max=156, state=_dvistate.inpage, args=('slen',))
    def _right_x(self, new_x):
        if new_x is not None:
            self.x = new_x
        self.h += self.x

    @_dispatch(min=157, max=160, state=_dvistate.inpage, args=('slen1',))
    def _down(self, a):
        self.v += a

    @_dispatch(min=161, max=165, state=_dvistate.inpage, args=('slen',))
    def _down_y(self, new_y):
        if new_y is not None:
            self.y = new_y
        self.v += self.y

    @_dispatch(min=166, max=170, state=_dvistate.inpage, args=('slen',))
    def _down_z(self, new_z):
        if new_z is not None:
            self.z = new_z
        self.v += self.z

    @_dispatch(min=171, max=234, state=_dvistate.inpage)
    def _fnt_num_immediate(self, k):
        self.f = k

    @_dispatch(min=235, max=238, state=_dvistate.inpage, args=('olen1',))
    def _fnt_num(self, new_f):
        self.f = new_f

    @_dispatch(min=239, max=242, args=('ulen1',))
    def _xxx(self, datalen):
        special = self.file.read(datalen)
        if six.PY3:
            chr_ = chr
        else:
            def chr_(x):
                return x
        _log.debug(
            'Dvi._xxx: encountered special: %s',
            ''.join([chr_(ch) if 32 <= ord(ch) < 127 else '<%02x>' % ord(ch)
                     for ch in special]))

    @_dispatch(min=243, max=246, args=('olen1', 'u4', 'u4', 'u4', 'u1', 'u1'))
    def _fnt_def(self, k, c, s, d, a, l):
        self._fnt_def_real(k, c, s, d, a, l)

    def _fnt_def_real(self, k, c, s, d, a, l):
        n = self.file.read(a + l)
        fontname = n[-l:].decode('ascii')
        tfm = _tfmfile(fontname)
        if tfm is None:
            if six.PY2:
                error_class = OSError
            else:
                error_class = FileNotFoundError
            raise error_class("missing font metrics file: %s" % fontname)
        if c != 0 and tfm.checksum != 0 and c != tfm.checksum:
            raise ValueError('tfm checksum mismatch: %s' % n)

        vf = _vffile(fontname)

        self.fonts[k] = DviFont(scale=s, tfm=tfm, texname=n, vf=vf)

    @_dispatch(247, state=_dvistate.pre, args=('u1', 'u4', 'u4', 'u4', 'u1'))
    def _pre(self, i, num, den, mag, k):
        comment = self.file.read(k)
        if i != 2:
            raise ValueError("Unknown dvi format %d" % i)
        if num != 25400000 or den != 7227 * 2**16:
            raise ValueError("nonstandard units in dvi file")
            # meaning: TeX always uses those exact values, so it
            # should be enough for us to support those
            # (There are 72.27 pt to an inch so 7227 pt =
            # 7227 * 2**16 sp to 100 in. The numerator is multiplied
            # by 10^5 to get units of 10**-7 meters.)
        if mag != 1000:
            raise ValueError("nonstandard magnification in dvi file")
            # meaning: LaTeX seems to frown on setting \mag, so
            # I think we can assume this is constant
        self.state = _dvistate.outer

    @_dispatch(248, state=_dvistate.outer)
    def _post(self, _):
        self.state = _dvistate.post_post
        # TODO: actually read the postamble and finale?
        # currently post_post just triggers closing the file

    @_dispatch(249)
    def _post_post(self, _):
        raise NotImplementedError

    @_dispatch(min=250, max=255)
    def _malformed(self, offset):
        raise ValueError("unknown command: byte %d", 250 + offset)


class DviFont(object):
    """
    Encapsulation of a font that a DVI file can refer to.

    This class holds a font's texname and size, supports comparison,
    and knows the widths of glyphs in the same units as the AFM file.
    There are also internal attributes (for use by dviread.py) that
    are *not* used for comparison.

    The size is in Adobe points (converted from TeX points).

    Parameters
    ----------

    scale : float
        Factor by which the font is scaled from its natural size.
    tfm : Tfm
        TeX font metrics for this font
    texname : bytes
       Name of the font as used internally by TeX and friends, as an
       ASCII bytestring. This is usually very different from any external
       font names, and :class:`dviread.PsfontsMap` can be used to find
       the external name of the font.
    vf : Vf
       A TeX "virtual font" file, or None if this font is not virtual.

    Attributes
    ----------

    texname : bytes
    size : float
       Size of the font in Adobe points, converted from the slightly
       smaller TeX points.
    widths : list
       Widths of glyphs in glyph-space units, typically 1/1000ths of
       the point size.

    """
    __slots__ = ('texname', 'size', 'widths', '_scale', '_vf', '_tfm')

    def __init__(self, scale, tfm, texname, vf):
        if not isinstance(texname, bytes):
            raise ValueError("texname must be a bytestring, got %s"
                             % type(texname))
        self._scale, self._tfm, self.texname, self._vf = \
            scale, tfm, texname, vf
        self.size = scale * (72.0 / (72.27 * 2**16))
        try:
            nchars = max(tfm.width) + 1
        except ValueError:
            nchars = 0
        self.widths = [(1000*tfm.width.get(char, 0)) >> 20
                       for char in xrange(nchars)]

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
            self.texname == other.texname and self.size == other.size

    def __ne__(self, other):
        return not self.__eq__(other)

    def _width_of(self, char):
        """
        Width of char in dvi units. For internal use by dviread.py.
        """

        width = self._tfm.width.get(char, None)
        if width is not None:
            return _mul2012(width, self._scale)
        _log.debug('No width for char %d in font %s.', char, self.texname)
        return 0

    def _height_depth_of(self, char):
        """
        Height and depth of char in dvi units. For internal use by dviread.py.
        """

        result = []
        for metric, name in ((self._tfm.height, "height"),
                             (self._tfm.depth, "depth")):
            value = metric.get(char, None)
            if value is None:
                _log.debug('No %s for char %d in font %s',
                           name, char, self.texname)
                result.append(0)
            else:
                result.append(_mul2012(value, self._scale))
        return result


class Vf(Dvi):
    """
    A virtual font (\\*.vf file) containing subroutines for dvi files.

    Usage::

      vf = Vf(filename)
      glyph = vf[code]
      glyph.text, glyph.boxes, glyph.width

    Parameters
    ----------

    filename : string or bytestring

    Notes
    -----

    The virtual font format is a derivative of dvi:
    http://mirrors.ctan.org/info/knuth/virtual-fonts
    This class reuses some of the machinery of `Dvi`
    but replaces the `_read` loop and dispatch mechanism.
    """

    def __init__(self, filename):
        Dvi.__init__(self, filename, 0)
        try:
            self._first_font = None
            self._chars = {}
            self._read()
        finally:
            self.close()

    def __getitem__(self, code):
        return self._chars[code]

    def _read(self):
        """
        Read one page from the file. Return True if successful,
        False if there were no more pages.
        """
        packet_char, packet_ends = None, None
        packet_len, packet_width = None, None
        while True:
            byte = ord(self.file.read(1)[0])
            # If we are in a packet, execute the dvi instructions
            if self.state == _dvistate.inpage:
                byte_at = self.file.tell()-1
                if byte_at == packet_ends:
                    self._finalize_packet(packet_char, packet_width)
                    packet_len, packet_char, packet_width = None, None, None
                    # fall through to out-of-packet code
                elif byte_at > packet_ends:
                    raise ValueError("Packet length mismatch in vf file")
                else:
                    if byte in (139, 140) or byte >= 243:
                        raise ValueError(
                            "Inappropriate opcode %d in vf file" % byte)
                    Dvi._dtable[byte](self, byte)
                    continue

            # We are outside a packet
            if byte < 242:          # a short packet (length given by byte)
                packet_len = byte
                packet_char, packet_width = self._arg(1), self._arg(3)
                packet_ends = self._init_packet(byte)
                self.state = _dvistate.inpage
            elif byte == 242:       # a long packet
                packet_len, packet_char, packet_width = \
                            [self._arg(x) for x in (4, 4, 4)]
                self._init_packet(packet_len)
            elif 243 <= byte <= 246:
                k = self._arg(byte - 242, byte == 246)
                c, s, d, a, l = [self._arg(x) for x in (4, 4, 4, 1, 1)]
                self._fnt_def_real(k, c, s, d, a, l)
                if self._first_font is None:
                    self._first_font = k
            elif byte == 247:       # preamble
                i, k = self._arg(1), self._arg(1)
                x = self.file.read(k)
                cs, ds = self._arg(4), self._arg(4)
                self._pre(i, x, cs, ds)
            elif byte == 248:       # postamble (just some number of 248s)
                break
            else:
                raise ValueError("unknown vf opcode %d" % byte)

    def _init_packet(self, pl):
        if self.state != _dvistate.outer:
            raise ValueError("Misplaced packet in vf file")
        self.h, self.v, self.w, self.x, self.y, self.z = 0, 0, 0, 0, 0, 0
        self.stack, self.text, self.boxes = [], [], []
        self.f = self._first_font
        return self.file.tell() + pl

    def _finalize_packet(self, packet_char, packet_width):
        self._chars[packet_char] = Page(
            text=self.text, boxes=self.boxes, width=packet_width,
            height=None, descent=None)
        self.state = _dvistate.outer

    def _pre(self, i, x, cs, ds):
        if self.state != _dvistate.pre:
            raise ValueError("pre command in middle of vf file")
        if i != 202:
            raise ValueError("Unknown vf format %d" % i)
        if len(x):
            _log.debug('vf file comment: %s', x)
        self.state = _dvistate.outer
        # cs = checksum, ds = design size


def _fix2comp(num):
    """
    Convert from two's complement to negative.
    """
    assert 0 <= num < 2**32
    if num & 2**31:
        return num - 2**32
    else:
        return num


def _mul2012(num1, num2):
    """
    Multiply two numbers in 20.12 fixed point format.
    """
    # Separated into a function because >> has surprising precedence
    return (num1*num2) >> 20


class Tfm(object):
    """
    A TeX Font Metric file.

    This implementation covers only the bare minimum needed by the Dvi class.

    Parameters
    ----------
    filename : string or bytestring

    Attributes
    ----------
    checksum : int
       Used for verifying against the dvi file.
    design_size : int
       Design size of the font (unknown units)
    width, height, depth : dict
       Dimensions of each character, need to be scaled by the factor
       specified in the dvi file. These are dicts because indexing may
       not start from 0.
    """
    __slots__ = ('checksum', 'design_size', 'width', 'height', 'depth')

    def __init__(self, filename):
        _log.debug('opening tfm file %s', filename)
        with open(filename, 'rb') as file:
            header1 = file.read(24)
            lh, bc, ec, nw, nh, nd = \
                struct.unpack('!6H', header1[2:14])
            _log.debug('lh=%d, bc=%d, ec=%d, nw=%d, nh=%d, nd=%d',
                       lh, bc, ec, nw, nh, nd)
            header2 = file.read(4*lh)
            self.checksum, self.design_size = \
                struct.unpack('!2I', header2[:8])
            # there is also encoding information etc.
            char_info = file.read(4*(ec-bc+1))
            widths = file.read(4*nw)
            heights = file.read(4*nh)
            depths = file.read(4*nd)

        self.width, self.height, self.depth = {}, {}, {}
        widths, heights, depths = \
            [struct.unpack('!%dI' % (len(x)/4), x)
             for x in (widths, heights, depths)]
        for idx, char in enumerate(xrange(bc, ec+1)):
            byte0 = ord(char_info[4*idx])
            byte1 = ord(char_info[4*idx+1])
            self.width[char] = _fix2comp(widths[byte0])
            self.height[char] = _fix2comp(heights[byte1 >> 4])
            self.depth[char] = _fix2comp(depths[byte1 & 0xf])


PsFont = namedtuple('Font', 'texname psname effects encoding filename')


class PsfontsMap(object):
    """
    A psfonts.map formatted file, mapping TeX fonts to PS fonts.

    Usage::

     >>> map = PsfontsMap(find_tex_file('pdftex.map'))
     >>> entry = map[b'ptmbo8r']
     >>> entry.texname
     b'ptmbo8r'
     >>> entry.psname
     b'Times-Bold'
     >>> entry.encoding
     '/usr/local/texlive/2008/texmf-dist/fonts/enc/dvips/base/8r.enc'
     >>> entry.effects
     {'slant': 0.16700000000000001}
     >>> entry.filename

    Parameters
    ----------

    filename : string or bytestring

    Notes
    -----

    For historical reasons, TeX knows many Type-1 fonts by different
    names than the outside world. (For one thing, the names have to
    fit in eight characters.) Also, TeX's native fonts are not Type-1
    but Metafont, which is nontrivial to convert to PostScript except
    as a bitmap. While high-quality conversions to Type-1 format exist
    and are shipped with modern TeX distributions, we need to know
    which Type-1 fonts are the counterparts of which native fonts. For
    these reasons a mapping is needed from internal font names to font
    file names.

    A texmf tree typically includes mapping files called e.g.
    :file:`psfonts.map`, :file:`pdftex.map`, or :file:`dvipdfm.map`.
    The file :file:`psfonts.map` is used by :program:`dvips`,
    :file:`pdftex.map` by :program:`pdfTeX`, and :file:`dvipdfm.map`
    by :program:`dvipdfm`. :file:`psfonts.map` might avoid embedding
    the 35 PostScript fonts (i.e., have no filename for them, as in
    the Times-Bold example above), while the pdf-related files perhaps
    only avoid the "Base 14" pdf fonts. But the user may have
    configured these files differently.
    """
    __slots__ = ('_font', '_filename')

    def __init__(self, filename):
        self._font = {}
        self._filename = filename
        if six.PY3 and isinstance(filename, bytes):
            encoding = sys.getfilesystemencoding() or 'utf-8'
            self._filename = filename.decode(encoding, errors='replace')
        with open(filename, 'rb') as file:
            self._parse(file)

    def __getitem__(self, texname):
        assert isinstance(texname, bytes)
        try:
            result = self._font[texname]
        except KeyError:
            fmt = ('A PostScript file for the font whose TeX name is "{0}" '
                   'could not be found in the file "{1}". The dviread module '
                   'can only handle fonts that have an associated PostScript '
                   'font file. '
                   'This problem can often be solved by installing '
                   'a suitable PostScript font package in your (TeX) '
                   'package manager.')
            msg = fmt.format(texname.decode('ascii'), self._filename)
            msg = textwrap.fill(msg, break_on_hyphens=False,
                                break_long_words=False)
            _log.info(msg)
            raise
        fn, enc = result.filename, result.encoding
        if fn is not None and not fn.startswith(b'/'):
            fn = find_tex_file(fn)
        if enc is not None and not enc.startswith(b'/'):
            enc = find_tex_file(result.encoding)
        return result._replace(filename=fn, encoding=enc)

    def _parse(self, file):
        """
        Parse the font mapping file.

        The format is, AFAIK: texname fontname [effects and filenames]
        Effects are PostScript snippets like ".177 SlantFont",
        filenames begin with one or two less-than signs. A filename
        ending in enc is an encoding file, other filenames are font
        files. This can be overridden with a left bracket: <[foobar
        indicates an encoding file named foobar.

        There is some difference between <foo.pfb and <<bar.pfb in
        subsetting, but I have no example of << in my TeX installation.
        """
        # If the map file specifies multiple encodings for a font, we
        # follow pdfTeX in choosing the last one specified. Such
        # entries are probably mistakes but they have occurred.
        # http://tex.stackexchange.com/questions/10826/
        # http://article.gmane.org/gmane.comp.tex.pdftex/4914

        empty_re = re.compile(br'%|\s*$')
        word_re = re.compile(
            br'''(?x) (?:
                 "<\[ (?P<enc1>  [^"]+    )" | # quoted encoding marked by [
                 "<   (?P<enc2>  [^"]+.enc)" | # quoted encoding, ends in .enc
                 "<<? (?P<file1> [^"]+    )" | # quoted font file name
                 "    (?P<eff1>  [^"]+    )" | # quoted effects or font name
                 <\[  (?P<enc3>  \S+      )  | # encoding marked by [
                 <    (?P<enc4>  \S+  .enc)  | # encoding, ends in .enc
                 <<?  (?P<file2> \S+      )  | # font file name
                      (?P<eff2>  \S+      )    # effects or font name
            )''')
        effects_re = re.compile(
            br'''(?x) (?P<slant> -?[0-9]*(?:\.[0-9]+)) \s* SlantFont
                    | (?P<extend>-?[0-9]*(?:\.[0-9]+)) \s* ExtendFont''')

        lines = (line.strip()
                 for line in file
                 if not empty_re.match(line))
        for line in lines:
            effects, encoding, filename = b'', None, None
            words = word_re.finditer(line)

            # The named groups are mutually exclusive and are
            # referenced below at an estimated order of probability of
            # occurrence based on looking at my copy of pdftex.map.
            # The font names are probably unquoted:
            w = next(words)
            texname = w.group('eff2') or w.group('eff1')
            w = next(words)
            psname = w.group('eff2') or w.group('eff1')

            for w in words:
                # Any effects are almost always quoted:
                eff = w.group('eff1') or w.group('eff2')
                if eff:
                    effects = eff
                    continue
                # Encoding files usually have the .enc suffix
                # and almost never need quoting:
                enc = (w.group('enc4') or w.group('enc3') or
                       w.group('enc2') or w.group('enc1'))
                if enc:
                    if encoding is not None:
                        _log.debug('Multiple encodings for %s = %s',
                                   texname, psname)
                    encoding = enc
                    continue
                # File names are probably unquoted:
                filename = w.group('file2') or w.group('file1')

            effects_dict = {}
            for match in effects_re.finditer(effects):
                slant = match.group('slant')
                if slant:
                    effects_dict['slant'] = float(slant)
                else:
                    effects_dict['extend'] = float(match.group('extend'))

            self._font[texname] = PsFont(
                texname=texname, psname=psname, effects=effects_dict,
                encoding=encoding, filename=filename)


class Encoding(object):
    """
    Parses a \\*.enc file referenced from a psfonts.map style file.
    The format this class understands is a very limited subset of
    PostScript.

    Usage (subject to change)::

      for name in Encoding(filename):
          whatever(name)

    Parameters
    ----------
    filename : string or bytestring

    Attributes
    ----------
    encoding : list
        List of character names
    """
    __slots__ = ('encoding',)

    def __init__(self, filename):
        with open(filename, 'rb') as file:
            _log.debug('Parsing TeX encoding %s', filename)
            self.encoding = self._parse(file)
            _log.debug('Result: %s', self.encoding)

    def __iter__(self):
        for name in self.encoding:
            yield name

    def _parse(self, file):
        result = []

        lines = (line.split(b'%', 1)[0].strip() for line in file)
        data = b''.join(lines)
        beginning = data.find(b'[')
        if beginning < 0:
            raise ValueError("Cannot locate beginning of encoding in {}"
                             .format(file))
        data = data[beginning:]
        end = data.find(b']')
        if end < 0:
            raise ValueError("Cannot locate end of encoding in {}"
                             .format(file))
        data = data[:end]

        return re.findall(br'/([^][{}<>\s]+)', data)


def find_tex_file(filename, format=None):
    """
    Find a file in the texmf tree.

    Calls :program:`kpsewhich` which is an interface to the kpathsea
    library [1]_. Most existing TeX distributions on Unix-like systems use
    kpathsea. It is also available as part of MikTeX, a popular
    distribution on Windows.

    Parameters
    ----------
    filename : string or bytestring
    format : string or bytestring
        Used as the value of the `--format` option to :program:`kpsewhich`.
        Could be e.g. 'tfm' or 'vf' to limit the search to that type of files.

    References
    ----------

    .. [1] `Kpathsea documentation <http://www.tug.org/kpathsea/>`_
        The library that :program:`kpsewhich` is part of.
    """

    if six.PY3:
        # we expect these to always be ascii encoded, but use utf-8
        # out of caution
        if isinstance(filename, bytes):
            filename = filename.decode('utf-8', errors='replace')
        if isinstance(format, bytes):
            format = format.decode('utf-8', errors='replace')

    cmd = ['kpsewhich']
    if format is not None:
        cmd += ['--format=' + format]
    cmd += [filename]
    _log.debug('find_tex_file(%s): %s', filename, cmd)
    # stderr is unused, but reading it avoids a subprocess optimization
    # that breaks EINTR handling in some Python versions:
    # http://bugs.python.org/issue12493
    # https://github.com/matplotlib/matplotlib/issues/633
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    result = pipe.communicate()[0].rstrip()
    _log.debug('find_tex_file result: %s', result)
    return result.decode('ascii')


# With multiple text objects per figure (e.g., tick labels) we may end
# up reading the same tfm and vf files many times, so we implement a
# simple cache. TODO: is this worth making persistent?

@lru_cache()
def _fontfile(cls, suffix, texname):
    filename = find_tex_file(texname + suffix)
    return cls(filename) if filename else None


_tfmfile = partial(_fontfile, Tfm, ".tfm")
_vffile = partial(_fontfile, Vf, ".vf")


if __name__ == '__main__':
    import sys
    fname = sys.argv[1]
    try:
        dpi = float(sys.argv[2])
    except IndexError:
        dpi = None
    with Dvi(fname, dpi) as dvi:
        fontmap = PsfontsMap(find_tex_file('pdftex.map'))
        for page in dvi:
            print('=== new page ===')
            fPrev = None
            for x, y, f, c, w in page.text:
                if f != fPrev:
                    print('font', f.texname, 'scaled', f._scale/pow(2.0, 20))
                    fPrev = f
                print(x, y, c, 32 <= c < 128 and chr(c) or '.', w)
            for x, y, w, h in page.boxes:
                print(x, y, 'BOX', w, h)
