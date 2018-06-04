# Copyright 2017 Palantir Technologies, Inc.
"""A collection of URI utilities with logic built on the VSCode URI library.

https://github.com/Microsoft/vscode-uri/blob/e59cab84f5df6265aed18ae5f43552d3eef13bb9/lib/index.ts
"""
import re
from urllib import parse
from pyls import IS_WIN

RE_DRIVE_LETTER_PATH = re.compile(r'^\/[a-zA-Z]:')


def urlparse(uri):
    """Parse and decode the parts of a URI."""
    scheme, netloc, path, params, query, fragment = parse.urlparse(uri)
    return (
        parse.unquote(scheme),
        parse.unquote(netloc),
        parse.unquote(path),
        parse.unquote(params),
        parse.unquote(query),
        parse.unquote(fragment)
    )


def urlunparse(parts):
    """Unparse and encode parts of a URI."""
    scheme, netloc, path, params, query, fragment = parts

    # Avoid encoding the windows drive letter colon
    if RE_DRIVE_LETTER_PATH.match(path):
        quoted_path = path[:3] + parse.quote(path[3:])
    else:
        quoted_path = parse.quote(path)

    return parse.urlunparse((
        parse.quote(scheme),
        parse.quote(netloc),
        quoted_path,
        parse.quote(params),
        parse.quote(query),
        parse.quote(fragment)
    ))


def to_fs_path(uri):
    """Returns the filesystem path of the given URI.

    Will handle UNC paths and normalize windows drive letters to lower-case. Also
    uses the platform specific path separator. Will *not* validate the path for
    invalid characters and semantics. Will *not* look at the scheme of this URI.
    """
    # scheme://netloc/path;parameters?query#fragment
    scheme, netloc, path, _params, _query, _fragment = urlparse(uri)

    if netloc and path and scheme == 'file':
        # unc path: file://shares/c$/far/boo
        value = "//{}{}".format(netloc, path)

    elif RE_DRIVE_LETTER_PATH.match(path):
        # windows drive letter: file:///C:/far/boo
        value = path[1].lower() + path[2:]

    else:
        # Other path
        value = path

    if IS_WIN:
        value = value.replace('/', '\\')

    return value


def from_fs_path(path):
    """Returns a URI for the given filesystem path."""
    scheme = 'file'
    params, query, fragment = '', '', ''
    path, netloc = _normalize_win_path(path)
    return urlunparse((scheme, netloc, path, params, query, fragment))


def uri_with(uri, scheme=None, netloc=None, path=None, params=None, query=None, fragment=None):
    """Return a URI with the given part(s) replaced.

    Parts are decoded / encoded.
    """
    old_scheme, old_netloc, old_path, old_params, old_query, old_fragment = urlparse(uri)
    path, _netloc = _normalize_win_path(path)
    return urlunparse((
        scheme or old_scheme,
        netloc or old_netloc,
        path or old_path,
        params or old_params,
        query or old_query,
        fragment or old_fragment
    ))


def _normalize_win_path(path):
    netloc = ''

    # normalize to fwd-slashes on windows,
    # on other systems bwd-slaches are valid
    # filename character, eg /f\oo/ba\r.txt
    if IS_WIN:
        path = path.replace('\\', '/')

    # check for authority as used in UNC shares
    # or use the path as given
    if path[:2] == '//':
        idx = path.index('/', 2)
        if idx == -1:
            netloc = path[2:]
        else:
            netloc = path[2:idx]
            path = path[idx:]
    else:
        path = path

    # Ensure that path starts with a slash
    # or that it is at least a slash
    if not path.startswith('/'):
        path = '/' + path

    # Normalize drive paths to lower case
    if RE_DRIVE_LETTER_PATH.match(path):
        path = path[0] + path[1].lower() + path[2:]

    return path, netloc
