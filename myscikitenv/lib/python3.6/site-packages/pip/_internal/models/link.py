import posixpath
import re

from pip._vendor.six.moves.urllib import parse as urllib_parse

from pip._internal.download import path_to_url
from pip._internal.utils.misc import (
    WHEEL_EXTENSION, redact_password_from_url, splitext,
)
from pip._internal.utils.models import KeyBasedCompareMixin
from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from typing import Optional, Tuple, Union
    from pip._internal.index import HTMLPage


class Link(KeyBasedCompareMixin):
    """Represents a parsed link from a Package Index's simple URL
    """

    def __init__(self, url, comes_from=None, requires_python=None):
        # type: (str, Optional[Union[str, HTMLPage]], Optional[str]) -> None
        """
        url:
            url of the resource pointed to (href of the link)
        comes_from:
            instance of HTMLPage where the link was found, or string.
        requires_python:
            String containing the `Requires-Python` metadata field, specified
            in PEP 345. This may be specified by a data-requires-python
            attribute in the HTML link tag, as described in PEP 503.
        """

        # url can be a UNC windows share
        if url.startswith('\\\\'):
            url = path_to_url(url)

        self.url = url
        self.comes_from = comes_from
        self.requires_python = requires_python if requires_python else None

        super(Link, self).__init__(
            key=(self.url),
            defining_class=Link
        )

    def __str__(self):
        if self.requires_python:
            rp = ' (requires-python:%s)' % self.requires_python
        else:
            rp = ''
        if self.comes_from:
            return '%s (from %s)%s' % (redact_password_from_url(self.url),
                                       self.comes_from, rp)
        else:
            return redact_password_from_url(str(self.url))

    def __repr__(self):
        return '<Link %s>' % self

    @property
    def filename(self):
        # type: () -> str
        _, netloc, path, _, _ = urllib_parse.urlsplit(self.url)
        name = posixpath.basename(path.rstrip('/')) or netloc
        name = urllib_parse.unquote(name)
        assert name, ('URL %r produced no filename' % self.url)
        return name

    @property
    def scheme(self):
        # type: () -> str
        return urllib_parse.urlsplit(self.url)[0]

    @property
    def netloc(self):
        # type: () -> str
        return urllib_parse.urlsplit(self.url)[1]

    @property
    def path(self):
        # type: () -> str
        return urllib_parse.unquote(urllib_parse.urlsplit(self.url)[2])

    def splitext(self):
        # type: () -> Tuple[str, str]
        return splitext(posixpath.basename(self.path.rstrip('/')))

    @property
    def ext(self):
        # type: () -> str
        return self.splitext()[1]

    @property
    def url_without_fragment(self):
        # type: () -> str
        scheme, netloc, path, query, fragment = urllib_parse.urlsplit(self.url)
        return urllib_parse.urlunsplit((scheme, netloc, path, query, None))

    _egg_fragment_re = re.compile(r'[#&]egg=([^&]*)')

    @property
    def egg_fragment(self):
        # type: () -> Optional[str]
        match = self._egg_fragment_re.search(self.url)
        if not match:
            return None
        return match.group(1)

    _subdirectory_fragment_re = re.compile(r'[#&]subdirectory=([^&]*)')

    @property
    def subdirectory_fragment(self):
        # type: () -> Optional[str]
        match = self._subdirectory_fragment_re.search(self.url)
        if not match:
            return None
        return match.group(1)

    _hash_re = re.compile(
        r'(sha1|sha224|sha384|sha256|sha512|md5)=([a-f0-9]+)'
    )

    @property
    def hash(self):
        # type: () -> Optional[str]
        match = self._hash_re.search(self.url)
        if match:
            return match.group(2)
        return None

    @property
    def hash_name(self):
        # type: () -> Optional[str]
        match = self._hash_re.search(self.url)
        if match:
            return match.group(1)
        return None

    @property
    def show_url(self):
        # type: () -> Optional[str]
        return posixpath.basename(self.url.split('#', 1)[0].split('?', 1)[0])

    @property
    def is_wheel(self):
        # type: () -> bool
        return self.ext == WHEEL_EXTENSION

    @property
    def is_artifact(self):
        # type: () -> bool
        """
        Determines if this points to an actual artifact (e.g. a tarball) or if
        it points to an "abstract" thing like a path or a VCS location.
        """
        from pip._internal.vcs import vcs

        if self.scheme in vcs.all_schemes:
            return False

        return True
