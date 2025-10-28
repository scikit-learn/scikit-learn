import copy

import numpy as np
from packaging.version import Version

from .._shared.filters import gaussian
from .._shared.utils import check_nD
from .brief_cy import _brief_loop
from .util import (
    DescriptorExtractor,
    _mask_border_keypoints,
    _prepare_grayscale_input_2D,
)


np2 = Version(np.__version__) >= Version('2')


class BRIEF(DescriptorExtractor):
    """BRIEF binary descriptor extractor.

    BRIEF (Binary Robust Independent Elementary Features) is an efficient
    feature point descriptor. It is highly discriminative even when using
    relatively few bits and is computed using simple intensity difference
    tests.

    For each keypoint, intensity comparisons are carried out for a specifically
    distributed number N of pixel-pairs resulting in a binary descriptor of
    length N. For binary descriptors the Hamming distance can be used for
    feature matching, which leads to lower computational cost in comparison to
    the L2 norm.

    Parameters
    ----------
    descriptor_size : int, optional
        Size of BRIEF descriptor for each keypoint. Sizes 128, 256 and 512
        recommended by the authors. Default is 256.
    patch_size : int, optional
        Length of the two dimensional square patch sampling region around
        the keypoints. Default is 49.
    mode : {'normal', 'uniform'}, optional
        Probability distribution for sampling location of decision pixel-pairs
        around keypoints.
    rng : {`numpy.random.Generator`, int}, optional
        Pseudo-random number generator (RNG).
        By default, a PCG64 generator is used (see :func:`numpy.random.default_rng`).
        If `rng` is an int, it is used to seed the generator.

        The PRNG is used for the random sampling of the decision
        pixel-pairs. From a square window with length `patch_size`,
        pixel pairs are sampled using the `mode` parameter to build
        the descriptors using intensity comparison.

        For matching across images, the same `rng` should be used to construct
        descriptors. To facilitate this:

        (a) `rng` defaults to 1
        (b) Subsequent calls of the ``extract`` method will use the same rng/seed.
    sigma : float, optional
        Standard deviation of the Gaussian low-pass filter applied to the image
        to alleviate noise sensitivity, which is strongly recommended to obtain
        discriminative and good descriptors.

    Attributes
    ----------
    descriptors : (Q, `descriptor_size`) array of dtype bool
        2D ndarray of binary descriptors of size `descriptor_size` for Q
        keypoints after filtering out border keypoints with value at an
        index ``(i, j)`` either being ``True`` or ``False`` representing
        the outcome of the intensity comparison for i-th keypoint on j-th
        decision pixel-pair. It is ``Q == np.sum(mask)``.
    mask : (N,) array of dtype bool
        Mask indicating whether a keypoint has been filtered out
        (``False``) or is described in the `descriptors` array (``True``).

    Examples
    --------
    >>> from skimage.feature import (corner_harris, corner_peaks, BRIEF,
    ...                              match_descriptors)
    >>> import numpy as np
    >>> square1 = np.zeros((8, 8), dtype=np.int32)
    >>> square1[2:6, 2:6] = 1
    >>> square1
    array([[0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0]], dtype=int32)
    >>> square2 = np.zeros((9, 9), dtype=np.int32)
    >>> square2[2:7, 2:7] = 1
    >>> square2
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=int32)
    >>> keypoints1 = corner_peaks(corner_harris(square1), min_distance=1)
    >>> keypoints2 = corner_peaks(corner_harris(square2), min_distance=1)
    >>> extractor = BRIEF(patch_size=5)
    >>> extractor.extract(square1, keypoints1)
    >>> descriptors1 = extractor.descriptors
    >>> extractor.extract(square2, keypoints2)
    >>> descriptors2 = extractor.descriptors
    >>> matches = match_descriptors(descriptors1, descriptors2)
    >>> matches
    array([[0, 0],
           [1, 1],
           [2, 2],
           [3, 3]])
    >>> keypoints1[matches[:, 0]]
    array([[2, 2],
           [2, 5],
           [5, 2],
           [5, 5]])
    >>> keypoints2[matches[:, 1]]
    array([[2, 2],
           [2, 6],
           [6, 2],
           [6, 6]])

    """

    def __init__(
        self, descriptor_size=256, patch_size=49, mode='normal', sigma=1, rng=1
    ):
        mode = mode.lower()
        if mode not in ('normal', 'uniform'):
            raise ValueError("`mode` must be 'normal' or 'uniform'.")

        self.descriptor_size = descriptor_size
        self.patch_size = patch_size
        self.mode = mode
        self.sigma = sigma

        if isinstance(rng, np.random.Generator):
            # Spawn an independent RNG from parent RNG provided by the user.
            # This is necessary so that we can safely deepcopy the RNG.
            # See https://github.com/scikit-learn/scikit-learn/issues/16988#issuecomment-1518037853
            bg = rng._bit_generator
            ss = bg._seed_seq
            (child_ss,) = ss.spawn(1)
            self.rng = np.random.Generator(type(bg)(child_ss))
        elif rng is None:
            self.rng = np.random.default_rng(np.random.SeedSequence())
        else:
            self.rng = np.random.default_rng(rng)

        self.descriptors = None
        self.mask = None

    def extract(self, image, keypoints):
        """Extract BRIEF binary descriptors for given keypoints in image.

        Parameters
        ----------
        image : 2D array
            Input image.
        keypoints : (N, 2) array
            Keypoint coordinates as ``(row, col)``.

        """
        check_nD(image, 2)

        # Copy RNG so we can repeatedly call extract with the same random values
        rng = copy.deepcopy(self.rng)

        image = _prepare_grayscale_input_2D(image)

        # Gaussian low-pass filtering to alleviate noise sensitivity
        image = np.ascontiguousarray(gaussian(image, sigma=self.sigma, mode='reflect'))

        # Sampling pairs of decision pixels in patch_size x patch_size window
        desc_size = self.descriptor_size
        patch_size = self.patch_size
        if self.mode == 'normal':
            samples = (patch_size / 5.0) * rng.standard_normal(desc_size * 8)
            samples = np.array(samples, dtype=np.int32)
            samples = samples[
                (samples < (patch_size // 2)) & (samples > -(patch_size - 2) // 2)
            ]

            pos1 = samples[: desc_size * 2].reshape(desc_size, 2)
            pos2 = samples[desc_size * 2 : desc_size * 4].reshape(desc_size, 2)
        elif self.mode == 'uniform':
            samples = rng.integers(
                -(patch_size - 2) // 2, (patch_size // 2) + 1, (desc_size * 2, 2)
            )
            samples = np.array(samples, dtype=np.int32)
            pos1, pos2 = np.split(samples, 2)

        pos1 = np.ascontiguousarray(pos1)
        pos2 = np.ascontiguousarray(pos2)

        # Removing keypoints that are within (patch_size / 2) distance from the
        # image border
        self.mask = _mask_border_keypoints(image.shape, keypoints, patch_size // 2)

        keypoints = np.array(
            keypoints[self.mask, :],
            dtype=np.int64,
            order='C',
            copy=None if np2 else False,
        )

        self.descriptors = np.zeros(
            (keypoints.shape[0], desc_size), dtype=bool, order='C'
        )

        _brief_loop(image, self.descriptors.view(np.uint8), keypoints, pos1, pos2)
