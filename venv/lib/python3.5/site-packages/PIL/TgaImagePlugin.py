#
# The Python Imaging Library.
# $Id$
#
# TGA file handling
#
# History:
# 95-09-01 fl   created (reads 24-bit files only)
# 97-01-04 fl   support more TGA versions, including compressed images
# 98-07-04 fl   fixed orientation and alpha layer bugs
# 98-09-11 fl   fixed orientation for runlength decoder
#
# Copyright (c) Secret Labs AB 1997-98.
# Copyright (c) Fredrik Lundh 1995-97.
#
# See the README file for information on usage and redistribution.
#


from . import Image, ImageFile, ImagePalette
from ._binary import i8, i16le as i16, o8, o16le as o16

__version__ = "0.3"


#
# --------------------------------------------------------------------
# Read RGA file


MODES = {
    # map imagetype/depth to rawmode
    (1, 8):  "P",
    (3, 1):  "1",
    (3, 8):  "L",
    (2, 16): "BGR;5",
    (2, 24): "BGR",
    (2, 32): "BGRA",
}


##
# Image plugin for Targa files.

class TgaImageFile(ImageFile.ImageFile):

    format = "TGA"
    format_description = "Targa"

    def _open(self):

        # process header
        s = self.fp.read(18)

        idlen = i8(s[0])

        colormaptype = i8(s[1])
        imagetype = i8(s[2])

        depth = i8(s[16])

        flags = i8(s[17])

        self.size = i16(s[12:]), i16(s[14:])

        # validate header fields
        if colormaptype not in (0, 1) or\
           self.size[0] <= 0 or self.size[1] <= 0 or\
           depth not in (1, 8, 16, 24, 32):
            raise SyntaxError("not a TGA file")

        # image mode
        if imagetype in (3, 11):
            self.mode = "L"
            if depth == 1:
                self.mode = "1"  # ???
        elif imagetype in (1, 9):
            self.mode = "P"
        elif imagetype in (2, 10):
            self.mode = "RGB"
            if depth == 32:
                self.mode = "RGBA"
        else:
            raise SyntaxError("unknown TGA mode")

        # orientation
        orientation = flags & 0x30
        if orientation == 0x20:
            orientation = 1
        elif not orientation:
            orientation = -1
        else:
            raise SyntaxError("unknown TGA orientation")

        self.info["orientation"] = orientation

        if imagetype & 8:
            self.info["compression"] = "tga_rle"

        if idlen:
            self.info["id_section"] = self.fp.read(idlen)

        if colormaptype:
            # read palette
            start, size, mapdepth = i16(s[3:]), i16(s[5:]), i16(s[7:])
            if mapdepth == 16:
                self.palette = ImagePalette.raw(
                    "BGR;16", b"\0"*2*start + self.fp.read(2*size))
            elif mapdepth == 24:
                self.palette = ImagePalette.raw(
                    "BGR", b"\0"*3*start + self.fp.read(3*size))
            elif mapdepth == 32:
                self.palette = ImagePalette.raw(
                    "BGRA", b"\0"*4*start + self.fp.read(4*size))

        # setup tile descriptor
        try:
            rawmode = MODES[(imagetype & 7, depth)]
            if imagetype & 8:
                # compressed
                self.tile = [("tga_rle", (0, 0)+self.size,
                              self.fp.tell(), (rawmode, orientation, depth))]
            else:
                self.tile = [("raw", (0, 0)+self.size,
                              self.fp.tell(), (rawmode, 0, orientation))]
        except KeyError:
            pass  # cannot decode

#
# --------------------------------------------------------------------
# Write TGA file


SAVE = {
    "1": ("1", 1, 0, 3),
    "L": ("L", 8, 0, 3),
    "P": ("P", 8, 1, 1),
    "RGB": ("BGR", 24, 0, 2),
    "RGBA": ("BGRA", 32, 0, 2),
}


def _save(im, fp, filename):

    try:
        rawmode, bits, colormaptype, imagetype = SAVE[im.mode]
    except KeyError:
        raise IOError("cannot write mode %s as TGA" % im.mode)

    if colormaptype:
        colormapfirst, colormaplength, colormapentry = 0, 256, 24
    else:
        colormapfirst, colormaplength, colormapentry = 0, 0, 0

    if im.mode == "RGBA":
        flags = 8
    else:
        flags = 0

    orientation = im.info.get("orientation", -1)
    if orientation > 0:
        flags = flags | 0x20

    fp.write(b"\000" +
             o8(colormaptype) +
             o8(imagetype) +
             o16(colormapfirst) +
             o16(colormaplength) +
             o8(colormapentry) +
             o16(0) +
             o16(0) +
             o16(im.size[0]) +
             o16(im.size[1]) +
             o8(bits) +
             o8(flags))

    if colormaptype:
        fp.write(im.im.getpalette("RGB", "BGR"))

    ImageFile._save(
        im, fp, [("raw", (0, 0) + im.size, 0, (rawmode, 0, orientation))])

    # write targa version 2 footer
    fp.write(b"\000" * 8 + b"TRUEVISION-XFILE." + b"\000")

#
# --------------------------------------------------------------------
# Registry


Image.register_open(TgaImageFile.format, TgaImageFile)
Image.register_save(TgaImageFile.format, _save)

Image.register_extension(TgaImageFile.format, ".tga")
