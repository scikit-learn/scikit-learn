#
# The Python Imaging Library.
# $Id$
#
# GD file handling
#
# History:
# 1996-04-12 fl   Created
#
# Copyright (c) 1997 by Secret Labs AB.
# Copyright (c) 1996 by Fredrik Lundh.
#
# See the README file for information on usage and redistribution.
#


# NOTE: This format cannot be automatically recognized, so the
# class is not registered for use with Image.open().  To open a
# gd file, use the GdImageFile.open() function instead.

# THE GD FORMAT IS NOT DESIGNED FOR DATA INTERCHANGE.  This
# implementation is provided for convenience and demonstrational
# purposes only.


from . import ImageFile, ImagePalette
from ._binary import i16be as i16
from ._util import isPath

__version__ = "0.1"

try:
    import builtins
except ImportError:
    import __builtin__
    builtins = __builtin__


##
# Image plugin for the GD uncompressed format.  Note that this format
# is not supported by the standard <b>Image.open</b> function.  To use
# this plugin, you have to import the <b>GdImageFile</b> module and
# use the <b>GdImageFile.open</b> function.

class GdImageFile(ImageFile.ImageFile):

    format = "GD"
    format_description = "GD uncompressed images"

    def _open(self):

        # Header
        s = self.fp.read(775)

        self.mode = "L"  # FIXME: "P"
        self.size = i16(s[0:2]), i16(s[2:4])

        # transparency index
        tindex = i16(s[5:7])
        if tindex < 256:
            self.info["transparent"] = tindex

        self.palette = ImagePalette.raw("RGB", s[7:])

        self.tile = [("raw", (0, 0)+self.size, 775, ("L", 0, -1))]


def open(fp, mode="r"):
    """
    Load texture from a GD image file.

    :param filename: GD file name, or an opened file handle.
    :param mode: Optional mode.  In this version, if the mode argument
        is given, it must be "r".
    :returns: An image instance.
    :raises IOError: If the image could not be read.
    """
    if mode != "r":
        raise ValueError("bad mode")

    if isPath(fp):
        filename = fp
        fp = builtins.open(fp, "rb")
    else:
        filename = ""

    try:
        return GdImageFile(fp, filename)
    except SyntaxError:
        raise IOError("cannot identify this image file")
