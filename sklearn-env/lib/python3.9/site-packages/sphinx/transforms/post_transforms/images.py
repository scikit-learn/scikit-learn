"""
    sphinx.transforms.post_transforms.images
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Docutils transforms used by Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import os
import re
from math import ceil
from typing import Any, Dict, List, Optional, Tuple

from docutils import nodes

from sphinx.application import Sphinx
from sphinx.locale import __
from sphinx.transforms import SphinxTransform
from sphinx.util import epoch_to_rfc1123, logging, requests, rfc1123_to_epoch, sha1
from sphinx.util.images import get_image_extension, guess_mimetype, parse_data_uri
from sphinx.util.osutil import ensuredir

logger = logging.getLogger(__name__)

MAX_FILENAME_LEN = 32
CRITICAL_PATH_CHAR_RE = re.compile('[:;<>|*" ]')


class BaseImageConverter(SphinxTransform):
    def apply(self, **kwargs: Any) -> None:
        for node in self.document.findall(nodes.image):
            if self.match(node):
                self.handle(node)

    def match(self, node: nodes.image) -> bool:
        return True

    def handle(self, node: nodes.image) -> None:
        pass

    @property
    def imagedir(self) -> str:
        return os.path.join(self.app.doctreedir, 'images')


class ImageDownloader(BaseImageConverter):
    default_priority = 100

    def match(self, node: nodes.image) -> bool:
        if self.app.builder.supported_image_types == []:
            return False
        elif self.app.builder.supported_remote_images:
            return False
        else:
            return '://' in node['uri']

    def handle(self, node: nodes.image) -> None:
        try:
            basename = os.path.basename(node['uri'])
            if '?' in basename:
                basename = basename.split('?')[0]
            if basename == '' or len(basename) > MAX_FILENAME_LEN:
                filename, ext = os.path.splitext(node['uri'])
                basename = sha1(filename.encode()).hexdigest() + ext
            basename = re.sub(CRITICAL_PATH_CHAR_RE, "_", basename)

            dirname = node['uri'].replace('://', '/').translate({ord("?"): "/",
                                                                 ord("&"): "/"})
            if len(dirname) > MAX_FILENAME_LEN:
                dirname = sha1(dirname.encode()).hexdigest()
            ensuredir(os.path.join(self.imagedir, dirname))
            path = os.path.join(self.imagedir, dirname, basename)

            headers = {}
            if os.path.exists(path):
                timestamp: float = ceil(os.stat(path).st_mtime)
                headers['If-Modified-Since'] = epoch_to_rfc1123(timestamp)

            r = requests.get(node['uri'], headers=headers)
            if r.status_code >= 400:
                logger.warning(__('Could not fetch remote image: %s [%d]') %
                               (node['uri'], r.status_code))
            else:
                self.app.env.original_image_uri[path] = node['uri']

                if r.status_code == 200:
                    with open(path, 'wb') as f:
                        f.write(r.content)

                last_modified = r.headers.get('last-modified')
                if last_modified:
                    timestamp = rfc1123_to_epoch(last_modified)
                    os.utime(path, (timestamp, timestamp))

                mimetype = guess_mimetype(path, default='*')
                if mimetype != '*' and os.path.splitext(basename)[1] == '':
                    # append a suffix if URI does not contain suffix
                    ext = get_image_extension(mimetype)
                    newpath = os.path.join(self.imagedir, dirname, basename + ext)
                    os.replace(path, newpath)
                    self.app.env.original_image_uri.pop(path)
                    self.app.env.original_image_uri[newpath] = node['uri']
                    path = newpath
                node['candidates'].pop('?')
                node['candidates'][mimetype] = path
                node['uri'] = path
                self.app.env.images.add_file(self.env.docname, path)
        except Exception as exc:
            logger.warning(__('Could not fetch remote image: %s [%s]') % (node['uri'], exc))


class DataURIExtractor(BaseImageConverter):
    default_priority = 150

    def match(self, node: nodes.image) -> bool:
        if self.app.builder.supported_remote_images == []:
            return False
        elif self.app.builder.supported_data_uri_images is True:
            return False
        else:
            return node['uri'].startswith('data:')

    def handle(self, node: nodes.image) -> None:
        image = parse_data_uri(node['uri'])
        ext = get_image_extension(image.mimetype)
        if ext is None:
            logger.warning(__('Unknown image format: %s...'), node['uri'][:32],
                           location=node)
            return

        ensuredir(os.path.join(self.imagedir, 'embeded'))
        digest = sha1(image.data).hexdigest()
        path = os.path.join(self.imagedir, 'embeded', digest + ext)
        self.app.env.original_image_uri[path] = node['uri']

        with open(path, 'wb') as f:
            f.write(image.data)

        node['candidates'].pop('?')
        node['candidates'][image.mimetype] = path
        node['uri'] = path
        self.app.env.images.add_file(self.env.docname, path)


def get_filename_for(filename: str, mimetype: str) -> str:
    basename = os.path.basename(filename)
    basename = re.sub(CRITICAL_PATH_CHAR_RE, "_", basename)
    return os.path.splitext(basename)[0] + get_image_extension(mimetype)


class ImageConverter(BaseImageConverter):
    """A base class for image converters.

    An image converter is kind of Docutils transform module.  It is used to
    convert image files which are not supported by a builder to the
    appropriate format for that builder.

    For example, :py:class:`LaTeX builder <.LaTeXBuilder>` supports PDF,
    PNG and JPEG as image formats.  However it does not support SVG images.
    For such case, using image converters allows to embed these
    unsupported images into the document.  One of the image converters;
    :ref:`sphinx.ext.imgconverter <sphinx.ext.imgconverter>` can convert
    a SVG image to PNG format using Imagemagick internally.

    There are three steps to make your custom image converter:

    1. Make a subclass of ``ImageConverter`` class
    2. Override ``conversion_rules``, ``is_available()`` and ``convert()``
    3. Register your image converter to Sphinx using
       :py:meth:`.Sphinx.add_post_transform`
    """
    default_priority = 200

    #: The converter is available or not.  Will be filled at the first call of
    #: the build.  The result is shared in the same process.
    #:
    #: .. todo:: This should be refactored not to store the state without class
    #:           variable.
    available: Optional[bool] = None

    #: A conversion rules the image converter supports.
    #: It is represented as a list of pair of source image format (mimetype) and
    #: destination one::
    #:
    #:     conversion_rules = [
    #:         ('image/svg+xml', 'image/png'),
    #:         ('image/gif', 'image/png'),
    #:         ('application/pdf', 'image/png'),
    #:     ]
    conversion_rules: List[Tuple[str, str]] = []

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    def match(self, node: nodes.image) -> bool:
        if not self.app.builder.supported_image_types:
            return False
        elif set(self.guess_mimetypes(node)) & set(self.app.builder.supported_image_types):
            # builder supports the image; no need to convert
            return False
        elif node['uri'].startswith('data:'):
            # all data URI MIME types are assumed to be supported
            return False
        elif self.available is None:
            # store the value to the class variable to share it during the build
            self.__class__.available = self.is_available()

        if not self.available:
            return False
        else:
            rule = self.get_conversion_rule(node)
            if rule:
                return True
            else:
                return False

    def get_conversion_rule(self, node: nodes.image) -> Tuple[str, str]:
        for candidate in self.guess_mimetypes(node):
            for supported in self.app.builder.supported_image_types:
                rule = (candidate, supported)
                if rule in self.conversion_rules:
                    return rule

        return None

    def is_available(self) -> bool:
        """Return the image converter is available or not."""
        raise NotImplementedError()

    def guess_mimetypes(self, node: nodes.image) -> List[str]:
        if '?' in node['candidates']:
            return []
        elif '*' in node['candidates']:
            return [guess_mimetype(node['uri'])]
        else:
            return node['candidates'].keys()

    def handle(self, node: nodes.image) -> None:
        _from, _to = self.get_conversion_rule(node)

        if _from in node['candidates']:
            srcpath = node['candidates'][_from]
        else:
            srcpath = node['candidates']['*']

        filename = get_filename_for(srcpath, _to)
        ensuredir(self.imagedir)
        destpath = os.path.join(self.imagedir, filename)

        abs_srcpath = os.path.join(self.app.srcdir, srcpath)
        if self.convert(abs_srcpath, destpath):
            if '*' in node['candidates']:
                node['candidates']['*'] = destpath
            else:
                node['candidates'][_to] = destpath
            node['uri'] = destpath

            self.env.original_image_uri[destpath] = srcpath
            self.env.images.add_file(self.env.docname, destpath)

    def convert(self, _from: str, _to: str) -> bool:
        """Convert an image file to the expected format.

        *_from* is a path of the source image file, and *_to* is a path
        of the destination file.
        """
        raise NotImplementedError()


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_post_transform(ImageDownloader)
    app.add_post_transform(DataURIExtractor)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
