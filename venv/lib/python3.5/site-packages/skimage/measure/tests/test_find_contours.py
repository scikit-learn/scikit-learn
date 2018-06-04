import numpy as np
from skimage.measure import find_contours

from skimage._shared import testing
from skimage._shared.testing import assert_array_equal


a = np.ones((8, 8), dtype=np.float32)
a[1:-1, 1] = 0
a[1, 1:-1] = 0

x, y = np.mgrid[-1:1:5j, -1:1:5j]
r = np.sqrt(x**2 + y**2)


def test_binary():
    ref = [[6. ,  1.5],
           [5. ,  1.5],
           [4. ,  1.5],
           [3. ,  1.5],
           [2. ,  1.5],
           [1.5,  2. ],
           [1.5,  3. ],
           [1.5,  4. ],
           [1.5,  5. ],
           [1.5,  6. ],
           [1. ,  6.5],
           [0.5,  6. ],
           [0.5,  5. ],
           [0.5,  4. ],
           [0.5,  3. ],
           [0.5,  2. ],
           [0.5,  1. ],
           [1. ,  0.5],
           [2. ,  0.5],
           [3. ,  0.5],
           [4. ,  0.5],
           [5. ,  0.5],
           [6. ,  0.5],
           [6.5,  1. ],
           [6. ,  1.5]]

    contours = find_contours(a, 0.5, positive_orientation='high')
    assert len(contours) == 1
    assert_array_equal(contours[0][::-1], ref)


def test_float():
    contours = find_contours(r, 0.5)
    assert len(contours) == 1
    assert_array_equal(contours[0],
                    [[ 2.,  3.],
                     [ 1.,  2.],
                     [ 2.,  1.],
                     [ 3.,  2.],
                     [ 2.,  3.]])


def test_memory_order():
    contours = find_contours(np.ascontiguousarray(r), 0.5)
    assert len(contours) == 1

    contours = find_contours(np.asfortranarray(r), 0.5)
    assert len(contours) == 1


def test_invalid_input():
    with testing.raises(ValueError):
        find_contours(r, 0.5, 'foo', 'bar')
    with testing.raises(ValueError):
        find_contours(r[..., None], 0.5)
