from __future__ import division, absolute_import, print_function

import os
import sys
import numpy as np
from . import util

from numpy.testing import run_module_suite, assert_array_equal, dec

def _path(*a):
    return os.path.join(*((os.path.dirname(__file__),) + a))

class TestCommonBlock(util.F2PyTest):
    sources = [_path('src', 'common', 'block.f')]

    @dec.knownfailureif(sys.platform=='win32', msg='Fails with MinGW64 Gfortran (Issue #9673)')
    def test_common_block(self):
        self.module.initcb()
        assert_array_equal(self.module.block.long_bn,
                           np.array(1.0, dtype=np.float64))
        assert_array_equal(self.module.block.string_bn,
                           np.array('2', dtype='|S1'))
        assert_array_equal(self.module.block.ok,
                           np.array(3, dtype=np.int32))

if __name__ == "__main__":
    run_module_suite()
