"""Stuff that differs in different Python versions and platform
distributions."""

# The following comment should be removed at some point in the future.
# mypy: disallow-untyped-defs=False

import codecs
import locale
import logging
import os
import sys

from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from typing import Optional, Union


__all__ = ["console_to_str", "get_path_uid", "stdlib_pkgs", "WINDOWS"]


logger = logging.getLogger(__name__)


def has_tls():
    # type: () -> bool
    try:
        import _ssl  # noqa: F401  # ignore unused
        return True
    except ImportError:
        pass

    from pip._vendor.urllib3.util import IS_PYOPENSSL
    return IS_PYOPENSSL


def str_to_display(data, desc=None):
    # type: (Union[bytes, str], Optional[str]) -> str
    """
    For display or logging purposes, convert a bytes object (or text) to
    text (e.g. unicode in Python 2) safe for output.

    :param desc: An optional phrase describing the input data, for use in
        the log message if a warning is logged. Defaults to "Bytes object".

    This function should never error out and so can take a best effort
    approach. It is okay to be lossy if needed since the return value is
    just for display.

    We assume the data is in the locale preferred encoding. If it won't
    decode properly, we warn the user but decode as best we can.

    We also ensure that the output can be safely written to standard output
    without encoding errors.
    """
    if isinstance(data, str):
        return data

    # Otherwise, data is a bytes object (str in Python 2).
    # First, get the encoding we assume. This is the preferred
    # encoding for the locale, unless that is not found, or
    # it is ASCII, in which case assume UTF-8
    encoding = locale.getpreferredencoding()
    if (not encoding) or codecs.lookup(encoding).name == "ascii":
        encoding = "utf-8"

    # Now try to decode the data - if we fail, warn the user and
    # decode with replacement.
    try:
        decoded_data = data.decode(encoding)
    except UnicodeDecodeError:
        logger.warning(
            '%s does not appear to be encoded as %s',
            desc or 'Bytes object',
            encoding,
        )
        decoded_data = data.decode(encoding, errors="backslashreplace")

    # Make sure we can print the output, by encoding it to the output
    # encoding with replacement of unencodable characters, and then
    # decoding again.
    # We use stderr's encoding because it's less likely to be
    # redirected and if we don't find an encoding we skip this
    # step (on the assumption that output is wrapped by something
    # that won't fail).
    # The double getattr is to deal with the possibility that we're
    # being called in a situation where sys.__stderr__ doesn't exist,
    # or doesn't have an encoding attribute. Neither of these cases
    # should occur in normal pip use, but there's no harm in checking
    # in case people use pip in (unsupported) unusual situations.
    output_encoding = getattr(getattr(sys, "__stderr__", None),
                              "encoding", None)

    if output_encoding:
        output_encoded = decoded_data.encode(
            output_encoding,
            errors="backslashreplace"
        )
        decoded_data = output_encoded.decode(output_encoding)

    return decoded_data


def console_to_str(data):
    # type: (bytes) -> str
    """Return a string, safe for output, of subprocess output.
    """
    return str_to_display(data, desc='Subprocess output')


def get_path_uid(path):
    # type: (str) -> int
    """
    Return path's uid.

    Does not follow symlinks:
        https://github.com/pypa/pip/pull/935#discussion_r5307003

    Placed this function in compat due to differences on AIX and
    Jython, that should eventually go away.

    :raises OSError: When path is a symlink or can't be read.
    """
    if hasattr(os, 'O_NOFOLLOW'):
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        file_uid = os.fstat(fd).st_uid
        os.close(fd)
    else:  # AIX and Jython
        # WARNING: time of check vulnerability, but best we can do w/o NOFOLLOW
        if not os.path.islink(path):
            # older versions of Jython don't have `os.fstat`
            file_uid = os.stat(path).st_uid
        else:
            # raise OSError for parity with os.O_NOFOLLOW above
            raise OSError(
                "{} is a symlink; Will not return uid for symlinks".format(
                    path)
            )
    return file_uid


# packages in the stdlib that may have installation metadata, but should not be
# considered 'installed'.  this theoretically could be determined based on
# dist.location (py27:`sysconfig.get_paths()['stdlib']`,
# py26:sysconfig.get_config_vars('LIBDEST')), but fear platform variation may
# make this ineffective, so hard-coding
stdlib_pkgs = {"python", "wsgiref", "argparse"}


# windows detection, covers cpython and ironpython
WINDOWS = (sys.platform.startswith("win") or
           (sys.platform == 'cli' and os.name == 'nt'))
