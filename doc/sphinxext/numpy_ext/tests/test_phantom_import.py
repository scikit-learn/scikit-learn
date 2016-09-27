from __future__ import division, absolute_import, print_function

import sys
from nose import SkipTest

def test_import():
    if sys.version_info[0] >= 3:
        raise SkipTest("phantom_import not ported to Py3")

    import numpydoc.phantom_import

# No tests at the moment...
