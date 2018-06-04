from __future__ import print_function, division

import numpy as np
from skimage.restoration import inpaint

from skimage._shared import testing
from skimage._shared.testing import assert_allclose


def test_inpaint_biharmonic_2d():
    img = np.tile(np.square(np.linspace(0, 1, 5)), (5, 1))
    mask = np.zeros_like(img)
    mask[2, 2:] = 1
    mask[1, 3:] = 1
    mask[0, 4:] = 1
    img[np.where(mask)] = 0
    out = inpaint.inpaint_biharmonic(img, mask)
    ref = np.array(
        [[0., 0.0625, 0.25000000, 0.5625000, 0.73925058],
         [0., 0.0625, 0.25000000, 0.5478048, 0.76557821],
         [0., 0.0625, 0.25842878, 0.5623079, 0.85927796],
         [0., 0.0625, 0.25000000, 0.5625000, 1.00000000],
         [0., 0.0625, 0.25000000, 0.5625000, 1.00000000]]
    )
    assert_allclose(ref, out)


def test_inpaint_biharmonic_3d():
    img = np.tile(np.square(np.linspace(0, 1, 5)), (5, 1))
    img = np.dstack((img, img.T))
    mask = np.zeros_like(img)
    mask[2, 2:, :] = 1
    mask[1, 3:, :] = 1
    mask[0, 4:, :] = 1
    img[np.where(mask)] = 0
    out = inpaint.inpaint_biharmonic(img, mask)
    ref = np.dstack((
        np.array(
            [[0.0000, 0.0625, 0.25000000, 0.56250000, 0.53752796],
             [0.0000, 0.0625, 0.25000000, 0.44443780, 0.53762210],
             [0.0000, 0.0625, 0.23693666, 0.46621112, 0.68615592],
             [0.0000, 0.0625, 0.25000000, 0.56250000, 1.00000000],
             [0.0000, 0.0625, 0.25000000, 0.56250000, 1.00000000]]),
        np.array(
            [[0.0000, 0.0000, 0.00000000, 0.00000000, 0.19621902],
             [0.0625, 0.0625, 0.06250000, 0.17470756, 0.30140091],
             [0.2500, 0.2500, 0.27241289, 0.35155440, 0.43068654],
             [0.5625, 0.5625, 0.56250000, 0.56250000, 0.56250000],
             [1.0000, 1.0000, 1.00000000, 1.00000000, 1.00000000]])
    ))
    assert_allclose(ref, out)


def test_invalid_input():
    img, mask = np.zeros([]), np.zeros([])
    with testing.raises(ValueError):
        inpaint.inpaint_biharmonic(img, mask)

    img, mask = np.zeros((2, 2)), np.zeros((4, 1))
    with testing.raises(ValueError):
        inpaint.inpaint_biharmonic(img, mask)

    img = np.ma.array(np.zeros((2, 2)), mask=[[0, 0], [0, 0]])
    mask = np.zeros((2, 2))
    with testing.raises(TypeError):
        inpaint.inpaint_biharmonic(img, mask)
