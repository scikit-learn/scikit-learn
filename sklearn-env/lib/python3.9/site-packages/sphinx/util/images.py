"""
    sphinx.util.images
    ~~~~~~~~~~~~~~~~~~

    Image utility functions for Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import base64
import imghdr
from collections import OrderedDict
from os import path
from typing import IO, BinaryIO, NamedTuple, Optional, Tuple

import imagesize

try:
    from PIL import Image
except ImportError:
    Image = None

mime_suffixes = OrderedDict([
    ('.gif', 'image/gif'),
    ('.jpg', 'image/jpeg'),
    ('.png', 'image/png'),
    ('.pdf', 'application/pdf'),
    ('.svg', 'image/svg+xml'),
    ('.svgz', 'image/svg+xml'),
    ('.ai', 'application/illustrator'),
])


class DataURI(NamedTuple):
    mimetype: str
    charset: str
    data: bytes


def get_image_size(filename: str) -> Optional[Tuple[int, int]]:
    try:
        size = imagesize.get(filename)
        if size[0] == -1:
            size = None
        elif isinstance(size[0], float) or isinstance(size[1], float):
            size = (int(size[0]), int(size[1]))

        if size is None and Image:  # fallback to Pillow
            with Image.open(filename) as im:
                size = im.size

        return size
    except Exception:
        return None


def guess_mimetype_for_stream(stream: IO, default: Optional[str] = None) -> Optional[str]:
    imgtype = imghdr.what(stream)
    if imgtype:
        return 'image/' + imgtype
    else:
        return default


def guess_mimetype(filename: str = '', default: Optional[str] = None) -> Optional[str]:
    _, ext = path.splitext(filename.lower())
    if ext in mime_suffixes:
        return mime_suffixes[ext]
    elif path.exists(filename):
        with open(filename, 'rb') as f:
            return guess_mimetype_for_stream(f, default=default)

    return default


def get_image_extension(mimetype: str) -> Optional[str]:
    for ext, _mimetype in mime_suffixes.items():
        if mimetype == _mimetype:
            return ext

    return None


def parse_data_uri(uri: str) -> Optional[DataURI]:
    if not uri.startswith('data:'):
        return None

    # data:[<MIME-type>][;charset=<encoding>][;base64],<data>
    mimetype = 'text/plain'
    charset = 'US-ASCII'

    properties, data = uri[5:].split(',', 1)
    for prop in properties.split(';'):
        if prop == 'base64':
            pass  # skip
        elif prop.startswith('charset='):
            charset = prop[8:]
        elif prop:
            mimetype = prop

    image_data = base64.b64decode(data)
    return DataURI(mimetype, charset, image_data)


def test_svg(h: bytes, f: Optional[BinaryIO]) -> Optional[str]:
    """An additional imghdr library helper; test the header is SVG's or not."""
    try:
        if '<svg' in h.decode().lower():
            return 'svg+xml'
    except UnicodeDecodeError:
        pass

    return None


# install test_svg() to imghdr
# refs: https://docs.python.org/3.6/library/imghdr.html#imghdr.tests
imghdr.tests.append(test_svg)
