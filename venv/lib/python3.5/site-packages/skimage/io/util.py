try:
    from urllib.request import urlopen  # Python 3
except ImportError:
    from urllib2 import urlopen  # Python 2

import os
import re
import tempfile
from contextlib import contextmanager
import six


URL_REGEX = re.compile(r'http://|https://|ftp://|file://|file:\\')


def is_url(filename):
    """Return True if string is an http or ftp path."""
    return (isinstance(filename, six.string_types) and
            URL_REGEX.match(filename) is not None)


@contextmanager
def file_or_url_context(resource_name):
    """Yield name of file from the given resource (i.e. file or url)."""
    if is_url(resource_name):
        _, ext = os.path.splitext(resource_name)
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=ext) as f:
                u = urlopen(resource_name)
                f.write(u.read())
            # f must be closed before yielding
            yield f.name
        finally:
            os.remove(f.name)
    else:
        yield resource_name
