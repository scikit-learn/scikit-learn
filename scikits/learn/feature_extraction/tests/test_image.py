# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
import scipy as sp
from scipy import ndimage

import nose

from ..image import img_to_graph, grid_to_graph
from ...utils.graph import cs_graph_components

def test_img_to_graph():
    x, y = np.mgrid[:4, :4] - 10
    grad_x = img_to_graph(x)
    grad_y = img_to_graph(y)
    nose.tools.assert_equal(grad_x.nnz, grad_y.nnz)
    # Negative elements are the diagonal: the elements of the original
    # image. Positive elements are the values of the gradient, they
    # shoudl all be equal on grad_x and grad_y
    np.testing.assert_array_equal(grad_x.data[grad_x.data > 0],
                                  grad_y.data[grad_y.data > 0])


def test_connect_regions():
    lena = sp.lena()
    for thr in (50, 150):
        mask = lena > thr
        graph = img_to_graph(lena, mask)
        nose.tools.assert_equal(ndimage.label(mask)[1],
                                cs_graph_components(graph)[0])


def test_connect_regions_with_grid():
    lena = sp.lena()
    mask = lena > 50
    graph = grid_to_graph(*lena.shape, mask=mask)
    nose.tools.assert_equal(ndimage.label(mask)[1],
                            cs_graph_components(graph)[0])

    mask = lena > 150
    graph = grid_to_graph(*lena.shape, mask=mask, dtype=None)
    nose.tools.assert_equal(ndimage.label(mask)[1],
                            cs_graph_components(graph)[0])

