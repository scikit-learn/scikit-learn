from skimage import transform
import numpy as np
from numpy import testing


def test_seam_carving():
    img = np.array([[0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 1],
                    [1, 0, 0, 1, 0]], dtype=np.float)

    expected = np.array([[0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 1],
                        [0, 0, 1, 0]], dtype=np.float)

    energy = 1 - img

    out = transform.seam_carve(img, energy, 'vertical', 1, border=0)
    testing.assert_equal(out, expected)

    img = img.T
    energy = energy.T

    out = transform.seam_carve(img, energy, 'horizontal', 1, border=0)
    testing.assert_equal(out, expected.T)


if __name__ == '__main__':
    np.testing.run_module_suite()
