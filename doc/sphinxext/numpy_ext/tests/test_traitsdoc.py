from __future__ import division, absolute_import, print_function

import sys
from nose import SkipTest

def test_import():
    if sys.version_info[0] >= 3:
        raise SkipTest("traitsdoc not ported to Python3")
    import numpydoc.traitsdoc

# No tests at the moment...
