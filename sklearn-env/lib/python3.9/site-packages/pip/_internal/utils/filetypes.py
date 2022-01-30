"""Filetype information.
"""
from pip._internal.utils.misc import splitext
from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from typing import Tuple

WHEEL_EXTENSION = '.whl'
BZ2_EXTENSIONS = ('.tar.bz2', '.tbz')  # type: Tuple[str, ...]
XZ_EXTENSIONS = ('.tar.xz', '.txz', '.tlz',
                 '.tar.lz', '.tar.lzma')  # type: Tuple[str, ...]
ZIP_EXTENSIONS = ('.zip', WHEEL_EXTENSION)  # type: Tuple[str, ...]
TAR_EXTENSIONS = ('.tar.gz', '.tgz', '.tar')  # type: Tuple[str, ...]
ARCHIVE_EXTENSIONS = (
    ZIP_EXTENSIONS + BZ2_EXTENSIONS + TAR_EXTENSIONS + XZ_EXTENSIONS
)


def is_archive_file(name):
    # type: (str) -> bool
    """Return True if `name` is a considered as an archive file."""
    ext = splitext(name)[1].lower()
    if ext in ARCHIVE_EXTENSIONS:
        return True
    return False
