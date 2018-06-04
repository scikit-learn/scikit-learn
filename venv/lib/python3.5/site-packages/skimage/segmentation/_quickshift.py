import numpy as np

import scipy.ndimage as ndi

from ..util import img_as_float
from ..color import rgb2lab

from ._quickshift_cy import _quickshift_cython


def quickshift(image, ratio=1.0, kernel_size=5, max_dist=10,
               return_tree=False, sigma=0, convert2lab=True, random_seed=42):
    """Segments image using quickshift clustering in Color-(x,y) space.

    Produces an oversegmentation of the image using the quickshift mode-seeking
    algorithm.

    Parameters
    ----------
    image : (width, height, channels) ndarray
        Input image.
    ratio : float, optional, between 0 and 1
        Balances color-space proximity and image-space proximity.
        Higher values give more weight to color-space.
    kernel_size : float, optional
        Width of Gaussian kernel used in smoothing the
        sample density. Higher means fewer clusters.
    max_dist : float, optional
        Cut-off point for data distances.
        Higher means fewer clusters.
    return_tree : bool, optional
        Whether to return the full segmentation hierarchy tree and distances.
    sigma : float, optional
        Width for Gaussian smoothing as preprocessing. Zero means no smoothing.
    convert2lab : bool, optional
        Whether the input should be converted to Lab colorspace prior to
        segmentation. For this purpose, the input is assumed to be RGB.
    random_seed : int, optional
        Random seed used for breaking ties.

    Returns
    -------
    segment_mask : (width, height) ndarray
        Integer mask indicating segment labels.

    Notes
    -----
    The authors advocate to convert the image to Lab color space prior to
    segmentation, though this is not strictly necessary. For this to work, the
    image must be given in RGB format.

    References
    ----------
    .. [1] Quick shift and kernel methods for mode seeking,
           Vedaldi, A. and Soatto, S.
           European Conference on Computer Vision, 2008
    """

    image = img_as_float(np.atleast_3d(image))
    if convert2lab:
        if image.shape[2] != 3:
            ValueError("Only RGB images can be converted to Lab space.")
        image = rgb2lab(image)

    if kernel_size < 1:
        raise ValueError("`kernel_size` should be >= 1.")

    image = ndi.gaussian_filter(image, [sigma, sigma, 0])
    image = np.ascontiguousarray(image * ratio)

    segment_mask = _quickshift_cython(
        image, kernel_size=kernel_size, max_dist=max_dist,
        return_tree=return_tree, random_seed=random_seed)
    return segment_mask
