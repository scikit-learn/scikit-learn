# -*- coding: utf-8 -*-
# Copyright (c) 2014-2020, imageio contributors
# imageio is distributed under the terms of the (new) BSD License.

# This docstring is used at the index of the documentation pages, and
# gets inserted into a slightly larger description (in setup.py) for
# the page on Pypi:
"""
Imageio is a Python library that provides an easy interface to read and
write a wide range of image data, including animated images, volumetric
data, and scientific formats. It is cross-platform, runs on Python 3.9+,
and is easy to install.

Main website: https://imageio.readthedocs.io/
"""

# flake8: noqa

__version__ = "2.36.1"

import warnings

# Load some bits from core
from .core import FormatManager, RETURN_BYTES

# Instantiate the old format manager
formats = FormatManager()
show_formats = formats.show

from . import v2
from . import v3
from . import plugins

# import config after core to avoid circular import
from . import config

# import all APIs into the top level (meta API)
from .v2 import (
    imread as imread_v2,
    mimread,
    volread,
    mvolread,
    imwrite,
    mimwrite,
    volwrite,
    mvolwrite,
    # aliases
    get_reader as read,
    get_writer as save,
    imwrite as imsave,
    mimwrite as mimsave,
    volwrite as volsave,
    mvolwrite as mvolsave,
    # misc
    help,
    get_reader,
    get_writer,
)
from .v3 import (
    imopen,
    # imread,  # Will take over once v3 is released
    # imwrite, # Will take over once v3 is released
    imiter,
)


def imread(uri, format=None, **kwargs):
    """imread(uri, format=None, **kwargs)

    Reads an image from the specified file. Returns a numpy array, which
    comes with a dict of meta data at its 'meta' attribute.

    Note that the image data is returned as-is, and may not always have
    a dtype of uint8 (and thus may differ from what e.g. PIL returns).

    Parameters
    ----------
    uri : {str, pathlib.Path, bytes, file}
        The resource to load the image from, e.g. a filename, pathlib.Path,
        http address or file object, see the docs for more info.
    format : str
        The format to use to read the file. By default imageio selects
        the appropriate for you based on the filename and its contents.
    kwargs : ...
        Further keyword arguments are passed to the reader. See :func:`.help`
        to see what arguments are available for a particular format.
    """

    warnings.warn(
        "Starting with ImageIO v3 the behavior of this function will switch to that of"
        " iio.v3.imread. To keep the current behavior (and make this warning disappear)"
        " use `import imageio.v2 as imageio` or call `imageio.v2.imread` directly.",
        DeprecationWarning,
        stacklevel=2,
    )

    return imread_v2(uri, format=format, **kwargs)


__all__ = [
    "v2",
    "v3",
    "config",
    "plugins",
    # v3 API
    "imopen",
    "imread",
    "imwrite",
    "imiter",
    # v2 API
    "mimread",
    "volread",
    "mvolread",
    "imwrite",
    "mimwrite",
    "volwrite",
    "mvolwrite",
    # v2 aliases
    "read",
    "save",
    "imsave",
    "mimsave",
    "volsave",
    "mvolsave",
    # functions to deprecate
    "help",
    "get_reader",
    "get_writer",
    "formats",
    "show_formats",
]
