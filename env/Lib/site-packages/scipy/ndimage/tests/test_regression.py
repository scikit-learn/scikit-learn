from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal

import scipy.ndimage as ndimage


def test_byte_order_median():
    """Regression test for #413: median_filter does not handle bytes orders."""
    a = np.arange(9, dtype='<f4').reshape(3, 3)
    ref = ndimage.filters.median_filter(a,(3, 3))
    b = np.arange(9, dtype='>f4').reshape(3, 3)
    t = ndimage.filters.median_filter(b, (3, 3))
    assert_array_almost_equal(ref, t)


def test_zoom_output_shape():
    """Ticket #643"""
    x = np.arange(12).reshape((3,4))
    ndimage.zoom(x, 2, output=np.zeros((6,8)))


def test_ticket_742():
    def SE(img, thresh=.7, size=4):
        mask = img > thresh
        rank = len(mask.shape)
        la, co = ndimage.label(mask,
                               ndimage.generate_binary_structure(rank, rank))
        slices = ndimage.find_objects(la)

    if np.dtype(np.intp) != np.dtype('i'):
        shape = (3,1240,1240)
        a = np.random.rand(np.prod(shape)).reshape(shape)
        # shouldn't crash
        SE(a)


def test_gh_issue_3025():
    """Github issue #3025 - improper merging of labels"""
    d = np.zeros((60,320))
    d[:,:257] = 1
    d[:,260:] = 1
    d[36,257] = 1
    d[35,258] = 1
    d[35,259] = 1
    assert ndimage.label(d, np.ones((3,3)))[1] == 1
