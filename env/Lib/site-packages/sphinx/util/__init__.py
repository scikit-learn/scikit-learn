"""Utility functions for Sphinx."""

from __future__ import annotations

import hashlib
import os
import posixpath
import re
from typing import Any

from sphinx.errors import ExtensionError as _ExtensionError
from sphinx.errors import FiletypeNotFoundError
from sphinx.util import _files, _importer, logging
from sphinx.util import index_entries as _index_entries
from sphinx.util._lines import parse_line_num_spec as parselinenos  # NoQA: F401
from sphinx.util._uri import encode_uri  # NoQA: F401
from sphinx.util._uri import is_url as isurl  # NoQA: F401
from sphinx.util.console import strip_colors  # NoQA: F401
from sphinx.util.matching import patfilter  # NoQA: F401
from sphinx.util.nodes import (  # NoQA: F401
    caption_ref_re,
    explicit_title_re,
    nested_parse_with_titles,
    split_explicit_title,
)

# import other utilities; partly for backwards compatibility, so don't
# prune unused ones indiscriminately
from sphinx.util.osutil import (  # NoQA: F401
    SEP,
    copyfile,
    ensuredir,
    make_filename,
    os_path,
    relative_uri,
)

logger = logging.getLogger(__name__)

# Generally useful regular expressions.
ws_re: re.Pattern[str] = re.compile(r'\s+')
url_re: re.Pattern[str] = re.compile(r'(?P<schema>.+)://.*')


# High-level utility functions.


def docname_join(basedocname: str, docname: str) -> str:
    return posixpath.normpath(posixpath.join('/' + basedocname, '..', docname))[1:]


def get_filetype(
    source_suffix: dict[str, str], filename: str | os.PathLike[str]
) -> str:
    for suffix, filetype in source_suffix.items():
        if os.fspath(filename).endswith(suffix):
            # If default filetype (None), considered as restructuredtext.
            return filetype or 'restructuredtext'
    raise FiletypeNotFoundError


def _md5(data: bytes = b'', **_kw: Any) -> hashlib._Hash:
    """Deprecated wrapper around hashlib.md5

    To be removed in Sphinx 9.0
    """
    return hashlib.md5(data, usedforsecurity=False)


def _sha1(data: bytes = b'', **_kw: Any) -> hashlib._Hash:
    """Deprecated wrapper around hashlib.sha1

    To be removed in Sphinx 9.0
    """
    return hashlib.sha1(data, usedforsecurity=False)


# deprecated name -> (object to return, canonical path or empty string)
_DEPRECATED_OBJECTS: dict[str, tuple[Any, str, tuple[int, int]]] = {
    'split_index_msg': (
        _index_entries.split_index_msg,
        'sphinx.util.index_entries.split_index_msg',
        (9, 0),
    ),
    'split_into': (
        _index_entries.split_index_msg,
        'sphinx.util.index_entries.split_into',
        (9, 0),
    ),
    'ExtensionError': (_ExtensionError, 'sphinx.errors.ExtensionError', (9, 0)),
    'md5': (_md5, '', (9, 0)),
    'sha1': (_sha1, '', (9, 0)),
    'import_object': (_importer.import_object, '', (10, 0)),
    'FilenameUniqDict': (_files.FilenameUniqDict, '', (10, 0)),
    'DownloadFiles': (_files.DownloadFiles, '', (10, 0)),
}


def __getattr__(name: str) -> Any:
    if name not in _DEPRECATED_OBJECTS:
        msg = f'module {__name__!r} has no attribute {name!r}'
        raise AttributeError(msg)

    from sphinx.deprecation import _deprecation_warning

    deprecated_object, canonical_name, remove = _DEPRECATED_OBJECTS[name]
    _deprecation_warning(__name__, name, canonical_name, remove=remove)
    return deprecated_object
