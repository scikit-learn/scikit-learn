from __future__ import division, print_function, absolute_import

import pytest
from numpy.testing import assert_array_equal
from scipy._lib._numpy_compat import suppress_warnings
import scipy.ndimage as ndi

import os

try:
    from PIL import Image
    pil_missing = False
except ImportError:
    pil_missing = True


@pytest.mark.skipif(pil_missing, reason="The Python Image Library could not be found.")
def test_imread():
    lp = os.path.join(os.path.dirname(__file__), 'dots.png')
    with suppress_warnings() as sup:
        # PIL causes a Py3k ResourceWarning
        sup.filter(message="unclosed file")
        sup.filter(DeprecationWarning)
        img = ndi.imread(lp, mode="RGB")
    assert_array_equal(img.shape, (300, 420, 3))

    with suppress_warnings() as sup:
        # PIL causes a Py3k ResourceWarning
        sup.filter(message="unclosed file")
        sup.filter(DeprecationWarning)
        img = ndi.imread(lp, flatten=True)
    assert_array_equal(img.shape, (300, 420))

    with open(lp, 'rb') as fobj:
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning)
            img = ndi.imread(fobj, mode="RGB")
        assert_array_equal(img.shape, (300, 420, 3))
