import numpy as np
import pytest

import skimage as ski
from skimage.transform import ThinPlateSplineTransform

SRC = np.array([[0, 0], [0, 5], [5, 5], [5, 0]])

DST = np.array([[5, 0], [0, 0], [0, 5], [5, 5]])


class TestThinPlateSplineTransform:
    def test_call_before_estimation(self):
        tps = ThinPlateSplineTransform()
        assert tps.src is None
        with pytest.raises(ValueError, match="Transformation is undefined"):
            tps(SRC)

    def test_call_invalid_coords_shape(self):
        tps = ThinPlateSplineTransform()
        tps.estimate(SRC, DST)
        coords = np.array([1, 2, 3])
        with pytest.raises(
            ValueError, match=r"Input `coords` must have shape \(N, 2\)"
        ):
            tps(coords)

    def test_call_on_SRC(self):
        tps = ThinPlateSplineTransform()
        tps.estimate(SRC, DST)
        result = tps(SRC)
        np.testing.assert_allclose(result, DST, atol=1e-15)

    def test_tps_transform_inverse(self):
        tps = ThinPlateSplineTransform()
        tps.estimate(SRC, DST)
        with pytest.raises(NotImplementedError):
            tps.inverse()

    def test_tps_estimation_faulty_input(self):
        src = np.array([[0, 0], [0, 5], [5, 5], [5, 0]])
        dst = np.array([[5, 0], [0, 0], [0, 5]])

        tps = ThinPlateSplineTransform()
        assert tps.src is None

        with pytest.raises(ValueError, match="Shape of `src` and `dst` didn't match"):
            tps.estimate(src, dst)

        less_than_3pts = np.array([[0, 0], [0, 5]])
        with pytest.raises(ValueError, match="Need at least 3 points"):
            tps.estimate(less_than_3pts, dst)
        with pytest.raises(ValueError, match="Need at least 3 points"):
            tps.estimate(src, less_than_3pts)
        with pytest.raises(ValueError, match="Need at least 3 points"):
            tps.estimate(less_than_3pts, less_than_3pts)

        not_2d = np.array([0, 1, 2, 3])
        with pytest.raises(ValueError, match=".*`src` must be a 2-dimensional array"):
            tps.estimate(not_2d, dst)
        with pytest.raises(ValueError, match=".*`dst` must be a 2-dimensional array"):
            tps.estimate(src, not_2d)

        # When the estimation fails, the instance attributes remain unchanged
        assert tps.src is None

    def test_rotate(self):
        image = ski.data.astronaut()
        desired = ski.transform.rotate(image, angle=90)

        src = np.array([[0, 0], [0, 511], [511, 511], [511, 0]])
        dst = np.array([[511, 0], [0, 0], [0, 511], [511, 511]])
        tps = ThinPlateSplineTransform()
        tps.estimate(src, dst)
        result = ski.transform.warp(image, tps)

        np.testing.assert_allclose(result, desired, atol=1e-13)
