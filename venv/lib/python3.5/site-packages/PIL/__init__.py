#
# The Python Imaging Library.
# $Id$
#
# package placeholder
#
# Copyright (c) 1999 by Secret Labs AB.
#
# See the README file for information on usage and redistribution.
#

# ;-)

from . import version

VERSION = '1.1.7'  # PIL Version
PILLOW_VERSION = version.__version__

__version__ = PILLOW_VERSION

_plugins = ['BlpImagePlugin',
            'BmpImagePlugin',
            'BufrStubImagePlugin',
            'CurImagePlugin',
            'DcxImagePlugin',
            'DdsImagePlugin',
            'EpsImagePlugin',
            'FitsStubImagePlugin',
            'FliImagePlugin',
            'FpxImagePlugin',
            'FtexImagePlugin',
            'GbrImagePlugin',
            'GifImagePlugin',
            'GribStubImagePlugin',
            'Hdf5StubImagePlugin',
            'IcnsImagePlugin',
            'IcoImagePlugin',
            'ImImagePlugin',
            'ImtImagePlugin',
            'IptcImagePlugin',
            'JpegImagePlugin',
            'Jpeg2KImagePlugin',
            'McIdasImagePlugin',
            'MicImagePlugin',
            'MpegImagePlugin',
            'MpoImagePlugin',
            'MspImagePlugin',
            'PalmImagePlugin',
            'PcdImagePlugin',
            'PcxImagePlugin',
            'PdfImagePlugin',
            'PixarImagePlugin',
            'PngImagePlugin',
            'PpmImagePlugin',
            'PsdImagePlugin',
            'SgiImagePlugin',
            'SpiderImagePlugin',
            'SunImagePlugin',
            'TgaImagePlugin',
            'TiffImagePlugin',
            'WebPImagePlugin',
            'WmfImagePlugin',
            'XbmImagePlugin',
            'XpmImagePlugin',
            'XVThumbImagePlugin']
