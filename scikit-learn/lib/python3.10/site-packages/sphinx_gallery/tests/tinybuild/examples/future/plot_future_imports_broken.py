"""
Test without ``__future__`` imports
-----------------------------------

Test that ``__future__`` imports inside sphinx_gallery modules does not affect the
parsing of this script.
"""

import sys

PY3_OR_LATER = sys.version_info[0] >= 3

# SyntaxError on Python 2
print(3 / 2, end="")

# Need to make this example fail on Python 3 as well (currently no way to say
# that an example is expected to fail only on Python 2)
if PY3_OR_LATER:
    raise RuntimeError("Forcing this example to fail on Python 3")
