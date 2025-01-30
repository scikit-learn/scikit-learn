import numpy as np
import pytest
from numpy.testing import assert_array_equal
import scipy.ndimage as ndi

from skimage import io, draw
from skimage._shared.testing import fetch
from skimage.data import binary_blobs
from skimage.morphology import medial_axis, skeletonize, thin
from skimage.morphology._skeletonize import G123_LUT, G123P_LUT, _generate_thin_luts


class TestSkeletonize:
    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_no_foreground(self, method):
        image = np.zeros((5, 5))
        result = skeletonize(image, method=method)
        assert_array_equal(result, np.zeros((5, 5)))

    @pytest.mark.parametrize(
        "ndim,method", [(1, "zhang"), (3, "zhang"), (1, "lee"), (4, "lee")]
    )
    def test_wrong_ndim(self, ndim, method):
        image = np.zeros((5,) * ndim, dtype=bool)
        with pytest.raises(ValueError):
            skeletonize(image, method=method)

    def test_wrong_method(self):
        image = np.ones((5, 5), dtype=bool)
        with pytest.raises(ValueError):
            skeletonize(image, method="foo")

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_skeletonize_all_foreground(self, method):
        image = np.ones((3, 4), dtype=bool)
        result = skeletonize(image, method=method)
        if method == "zhang":
            expected = np.array([[0, 0, 1, 0], [1, 1, 0, 0], [0, 0, 0, 0]], dtype=bool)
        else:  # "lee"
            expected = np.array([[0, 0, 0, 0], [1, 1, 1, 1], [0, 0, 0, 0]], dtype=bool)
        assert_array_equal(result, expected)

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_single_point(self, method):
        image = np.zeros((5, 5), dtype=bool)
        image[3, 3] = 1
        result = skeletonize(image, method=method)
        assert_array_equal(result, image)

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_vec_1d(self, method):
        # Corner case of a 2D image, which is a 1D vector
        image = np.ones((5, 1), dtype=bool)
        result = skeletonize(image, method=method)
        assert_array_equal(result, image)

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_already_thinned(self, method):
        image = np.array(
            [
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1],
                [0, 1, 1, 1, 0],
                [1, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        result = skeletonize(image, method=method)
        assert_array_equal(result, image)

    def test_output(self):
        image = io.imread(fetch("data/bw_text.png"), as_gray=True)

        # make black the foreground
        image = image == 0
        result = skeletonize(image)

        expected = np.load(fetch("data/bw_text_skeleton.npy"))
        assert_array_equal(result, expected)

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    @pytest.mark.parametrize("dtype", [bool, float, int])
    def test_num_neighbors(self, method, dtype):
        # an empty image
        image = np.zeros((300, 300), dtype=dtype)

        # foreground object 1
        image[10:-10, 10:100] = 1
        image[-100:-10, 10:-10] = 2
        image[10:-10, -100:-10] = 3

        # foreground object 2
        rs, cs = draw.line(250, 150, 10, 280)
        for i in range(10):
            image[rs + i, cs] = 4
        rs, cs = draw.line(10, 150, 250, 280)
        for i in range(20):
            image[rs + i, cs] = 5

        # foreground object 3
        ir, ic = np.indices(image.shape)
        circle1 = (ic - 135) ** 2 + (ir - 150) ** 2 < 30**2
        circle2 = (ic - 135) ** 2 + (ir - 150) ** 2 < 20**2
        image[circle1] = 1
        image[circle2] = 0
        result = skeletonize(image, method=method).astype(np.uint8)

        # there should never be a 2x2 block of foreground pixels in a skeleton
        mask = np.array([[1, 1], [1, 1]], np.uint8)
        blocks = ndi.correlate(result, mask, mode="constant")
        assert not np.any(blocks == 4)

    def test_lut_fix(self):
        image = np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 1, 1, 0, 0],
                [0, 0, 0, 1, 1, 0],
                [0, 0, 0, 0, 1, 1],
                [0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        result = skeletonize(image)
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        assert np.all(result == expected)

    @pytest.mark.parametrize("ndim,method", [(2, "zhang"), (2, "lee"), (3, "lee")])
    @pytest.mark.parametrize("dtype", [bool, np.uint8])
    def test_input_not_modified(self, method, ndim, dtype):
        # Skeletonize must not modify the input image
        image = np.ones((3,) * ndim, dtype=dtype)
        image = np.pad(image, 1)
        original = image.copy()
        _ = skeletonize(image, method=method)
        np.testing.assert_array_equal(image, original)

    @pytest.mark.parametrize("method", ["zhang", "lee"])
    def test_input_float_conv(self, method):
        # Check that the floats are correctly handled. Also check non-contiguous input
        image = np.random.random((16, 16))[::2, ::2]
        image[image < 0.5] = 0.0

        original = image.copy()
        result = skeletonize(image, method=method)

        assert result.dtype == bool
        assert_array_equal(image, original)

    def test_two_hole_image_vs_fiji(self):
        # Test a simple 2D image against FIJI
        image = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0],
                [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
                [0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
                [0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
                [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
                [0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        result = skeletonize(image, method="lee")
        assert_array_equal(result, expected)

    def test_3d_vs_fiji(self):
        # Generate an image with blobs and compare its skeleton
        # to the one generated by FIJI (Plugins>Skeleton->Skeletonize)
        image = binary_blobs(32, 0.05, n_dim=3, rng=1234)
        image = image[:-2, ...]

        result = skeletonize(image)
        expected = io.imread(fetch("data/_blobs_3d_fiji_skeleton.tif")).astype(bool)
        assert_array_equal(result, expected)


class TestThin:
    @property
    def input_image(self):
        # Image to test thinning with
        ii = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 1, 2, 3, 4, 5, 0],
                [0, 1, 0, 1, 1, 1, 0],
                [0, 1, 1, 1, 1, 1, 0],
                [0, 6, 1, 1, 1, 1, 0],
                [0, 1, 1, 1, 1, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=float,
        )
        return ii

    def test_all_zeros(self):
        image = np.zeros((10, 10), dtype=bool)
        assert np.all(thin(image) == False)

    @pytest.mark.parametrize("dtype", [bool, float, int])
    def test_thin_copies_input(self, dtype):
        """Ensure thinning does not modify the input image."""
        image = self.input_image.astype(dtype)
        original = image.copy()
        thin(image)
        np.testing.assert_array_equal(image, original)

    @pytest.mark.parametrize("dtype", [bool, float, int])
    def test_iter_1(self, dtype):
        image = self.input_image.astype(dtype)
        result = thin(image, 1).astype(bool)
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0],
                [0, 1, 0, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        assert_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", [bool, float, int])
    def test_noiter(self, dtype):
        image = self.input_image.astype(dtype)
        result = thin(image).astype(bool)
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0],
                [0, 1, 0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        assert_array_equal(result, expected)

    def test_baddim(self):
        for ii in [np.zeros(3, dtype=bool), np.zeros((3, 3, 3), dtype=bool)]:
            with pytest.raises(ValueError):
                thin(ii)

    def test_lut_generation(self):
        g123, g123p = _generate_thin_luts()

        assert_array_equal(g123, G123_LUT)
        assert_array_equal(g123p, G123P_LUT)


class TestMedialAxis:
    def test_all_zeros(self):
        result = medial_axis(np.zeros((10, 10), dtype=bool))
        assert np.all(result == False)

    def test_all_zeros_masked(self):
        result = medial_axis(
            np.zeros((10, 10), dtype=bool), np.zeros((10, 10), dtype=bool)
        )
        assert np.all(result == False)

    @pytest.mark.parametrize("dtype", [bool, float, int])
    def test_vertical_line(self, dtype):
        # Image is a thick vertical line (see gh-3861)
        image = np.zeros((9, 9), dtype=dtype)
        image[:, 2] = 1
        image[:, 3] = 2
        image[:, 4] = 3

        expected = np.full(image.shape, False)
        expected[:, 3] = True

        result = medial_axis(image)
        assert_array_equal(result, expected)

    def test_rectangle(self):
        image = np.zeros((9, 15), dtype=bool)
        image[1:-1, 1:-1] = True
        # Excepted are four diagonals from the corners, meeting in a horizontal line
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        result = medial_axis(image)
        assert np.all(result == expected)
        result, distance = medial_axis(image, return_distance=True)
        assert distance.max() == 4

    def test_rectange_with_hole(self):
        image = np.zeros((9, 15), dtype=bool)
        image[1:-1, 1:-1] = True
        image[4, 4:-4] = False
        expected = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=bool,
        )
        result = medial_axis(image)
        assert np.all(result == expected)

    def test_narrow_image(self):
        # Image is a 1-pixel thin strip
        image = np.zeros((1, 5), dtype=bool)
        image[:, 1:-1] = True
        result = medial_axis(image)
        assert np.all(result == image)
