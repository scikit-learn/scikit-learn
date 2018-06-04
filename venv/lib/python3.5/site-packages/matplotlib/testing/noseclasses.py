"""
The module testing.noseclasses is deprecated as of 2.1
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
try:
    from ._nose.plugins.knownfailure import KnownFailure as _KnownFailure
    has_nose = True
except ImportError:
    has_nose = False
    _KnownFailure = object

from .. import cbook

cbook.warn_deprecated(
    since="2.1",
    message="The noseclass module has been deprecated in 2.1 and will "
            "be removed in matplotlib 2.3.")


@cbook.deprecated("2.1")
class KnownFailure(_KnownFailure):
    def __init__(self):
        if not has_nose:
            raise ImportError("Need nose for this plugin.")
