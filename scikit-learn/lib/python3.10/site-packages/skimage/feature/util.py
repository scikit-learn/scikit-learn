import numpy as np

from ..util import img_as_float
from .._shared.utils import (
    _supported_float_type,
    check_nD,
)


class FeatureDetector:
    def __init__(self):
        self.keypoints_ = np.array([])

    def detect(self, image):
        """Detect keypoints in image.

        Parameters
        ----------
        image : 2D array
            Input image.

        """
        raise NotImplementedError()


class DescriptorExtractor:
    def __init__(self):
        self.descriptors_ = np.array([])

    def extract(self, image, keypoints):
        """Extract feature descriptors in image for given keypoints.

        Parameters
        ----------
        image : 2D array
            Input image.
        keypoints : (N, 2) array
            Keypoint locations as ``(row, col)``.

        """
        raise NotImplementedError()


def plot_matched_features(
    image0,
    image1,
    *,
    keypoints0,
    keypoints1,
    matches,
    ax,
    keypoints_color='k',
    matches_color=None,
    only_matches=False,
    alignment='horizontal',
):
    """Plot matched features between two images.

    .. versionadded:: 0.23

    Parameters
    ----------
    image0 : (N, M [, 3]) array
        First image.
    image1 : (N, M [, 3]) array
        Second image.
    keypoints0 : (K1, 2) array
        First keypoint coordinates as ``(row, col)``.
    keypoints1 : (K2, 2) array
        Second keypoint coordinates as ``(row, col)``.
    matches : (Q, 2) array
        Indices of corresponding matches in first and second sets of
        descriptors, where `matches[:, 0]` (resp. `matches[:, 1]`) contains
        the indices in the first (resp. second) set of descriptors.
    ax : matplotlib.axes.Axes
        The Axes object where the images and their matched features are drawn.
    keypoints_color : matplotlib color, optional
        Color for keypoint locations.
    matches_color : matplotlib color or sequence thereof, optional
        Single color or sequence of colors for each line defined by `matches`,
        which connect keypoint matches. See [1]_ for an overview of supported
        color formats. By default, colors are picked randomly.
    only_matches : bool, optional
        Set to True to plot matches only and not the keypoint locations.
    alignment : {'horizontal', 'vertical'}, optional
        Whether to show the two images side by side (`'horizontal'`), or one above
        the other (`'vertical'`).

    References
    ----------
    .. [1] https://matplotlib.org/stable/users/explain/colors/colors.html#specifying-colors

    Notes
    -----
    To make a sequence of colors passed to `matches_color` work for any number of
    `matches`, you can wrap that sequence in :func:`itertools.cycle`.
    """
    image0 = img_as_float(image0)
    image1 = img_as_float(image1)

    new_shape0 = list(image0.shape)
    new_shape1 = list(image1.shape)

    if image0.shape[0] < image1.shape[0]:
        new_shape0[0] = image1.shape[0]
    elif image0.shape[0] > image1.shape[0]:
        new_shape1[0] = image0.shape[0]

    if image0.shape[1] < image1.shape[1]:
        new_shape0[1] = image1.shape[1]
    elif image0.shape[1] > image1.shape[1]:
        new_shape1[1] = image0.shape[1]

    if new_shape0 != image0.shape:
        new_image0 = np.zeros(new_shape0, dtype=image0.dtype)
        new_image0[: image0.shape[0], : image0.shape[1]] = image0
        image0 = new_image0

    if new_shape1 != image1.shape:
        new_image1 = np.zeros(new_shape1, dtype=image1.dtype)
        new_image1[: image1.shape[0], : image1.shape[1]] = image1
        image1 = new_image1

    offset = np.array(image0.shape)
    if alignment == 'horizontal':
        image = np.concatenate([image0, image1], axis=1)
        offset[0] = 0
    elif alignment == 'vertical':
        image = np.concatenate([image0, image1], axis=0)
        offset[1] = 0
    else:
        mesg = (
            f"`plot_matched_features` accepts either 'horizontal' or 'vertical' for "
            f"alignment, but '{alignment}' was given. See "
            f"https://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.plot_matched_features "
            f"for details."
        )
        raise ValueError(mesg)

    if not only_matches:
        ax.scatter(
            keypoints0[:, 1],
            keypoints0[:, 0],
            facecolors='none',
            edgecolors=keypoints_color,
        )
        ax.scatter(
            keypoints1[:, 1] + offset[1],
            keypoints1[:, 0] + offset[0],
            facecolors='none',
            edgecolors=keypoints_color,
        )

    ax.imshow(image, cmap='gray')
    ax.axis((0, image0.shape[1] + offset[1], image0.shape[0] + offset[0], 0))

    number_of_matches = matches.shape[0]

    from matplotlib.colors import is_color_like

    if matches_color is None:
        rng = np.random.default_rng(seed=0)
        colors = [rng.random(3) for _ in range(number_of_matches)]
    elif is_color_like(matches_color):
        colors = [matches_color for _ in range(number_of_matches)]
    elif hasattr(matches_color, "__len__") and len(matches_color) == number_of_matches:
        # No need to check each color, matplotlib does so for us
        colors = matches_color
    else:
        error_message = (
            '`matches_color` needs to be a single color '
            'or a sequence of length equal to the number of matches.'
        )
        raise ValueError(error_message)

    for i, match in enumerate(matches):
        idx0, idx1 = match
        ax.plot(
            (keypoints0[idx0, 1], keypoints1[idx1, 1] + offset[1]),
            (keypoints0[idx0, 0], keypoints1[idx1, 0] + offset[0]),
            '-',
            color=colors[i],
        )


def _prepare_grayscale_input_2D(image):
    image = np.squeeze(image)
    check_nD(image, 2)
    image = img_as_float(image)
    float_dtype = _supported_float_type(image.dtype)
    return image.astype(float_dtype, copy=False)


def _prepare_grayscale_input_nD(image):
    image = np.squeeze(image)
    check_nD(image, range(2, 6))
    image = img_as_float(image)
    float_dtype = _supported_float_type(image.dtype)
    return image.astype(float_dtype, copy=False)


def _mask_border_keypoints(image_shape, keypoints, distance):
    """Mask coordinates that are within certain distance from the image border.

    Parameters
    ----------
    image_shape : (2,) array_like
        Shape of the image as ``(rows, cols)``.
    keypoints : (N, 2) array
        Keypoint coordinates as ``(rows, cols)``.
    distance : int
        Image border distance.

    Returns
    -------
    mask : (N,) bool array
        Mask indicating if pixels are within the image (``True``) or in the
        border region of the image (``False``).

    """

    rows = image_shape[0]
    cols = image_shape[1]

    mask = (
        ((distance - 1) < keypoints[:, 0])
        & (keypoints[:, 0] < (rows - distance + 1))
        & ((distance - 1) < keypoints[:, 1])
        & (keypoints[:, 1] < (cols - distance + 1))
    )

    return mask
