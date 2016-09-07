from __future__ import division, absolute_import, print_function

import sys
from nose import SkipTest

def test_import():
    if sys.version_info[0] >= 3:
        raise SkipTest("plot_directive not ported to Python 3 (use the one from Matplotlib instead)")
    import numpydoc.plot_directive

# No tests at the moment...
