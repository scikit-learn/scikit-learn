"""
canny.py - Canny Edge detector

Reference: Canny, J., A Computational Approach To Edge Detection, IEEE Trans.
    Pattern Analysis and Machine Intelligence, 8:679-714, 1986
"""

import numpy as np
import scipy.ndimage as ndi

from ..util.dtype import dtype_limits
from .._shared.filters import gaussian
from .._shared.utils import _supported_float_type, check_nD
from ._canny_cy import _nonmaximum_suppression_bilinear


def _preprocess(image, mask, sigma, mode, cval):
    """Generate a smoothed image and an eroded mask.

    The image is smoothed using a gaussian filter ignoring masked
    pixels and the mask is eroded.

    Parameters
    ----------
    image : array
        Image to be smoothed.
    mask : array
        Mask with 1's for significant pixels, 0's for masked pixels.
    sigma : scalar or sequence of scalars
        Standard deviation for Gaussian kernel. The standard
        deviations of the Gaussian filter are given for each axis as a
        sequence, or as a single number, in which case it is equal for
        all axes.
    mode : str, {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}
        The ``mode`` parameter determines how the array borders are
        handled, where ``cval`` is the value when mode is equal to
        'constant'.
    cval : float, optional
        Value to fill past edges of input if `mode` is 'constant'.

    Returns
    -------
    smoothed_image : ndarray
        The smoothed array
    eroded_mask : ndarray
        The eroded mask.

    Notes
    -----
    This function calculates the fractional contribution of masked pixels
    by applying the function to the mask (which gets you the fraction of
    the pixel data that's due to significant points). We then mask the image
    and apply the function. The resulting values will be lower by the
    bleed-over fraction, so you can recalibrate by dividing by the function
    on the mask to recover the effect of smoothing from just the significant
    pixels.
    """
    gaussian_kwargs = dict(sigma=sigma, mode=mode, cval=cval, preserve_range=False)
    compute_bleedover = mode == 'constant' or mask is not None
    float_type = _supported_float_type(image.dtype)
    if mask is None:
        if compute_bleedover:
            mask = np.ones(image.shape, dtype=float_type)
        masked_image = image

        eroded_mask = np.ones(image.shape, dtype=bool)
        eroded_mask[:1, :] = 0
        eroded_mask[-1:, :] = 0
        eroded_mask[:, :1] = 0
        eroded_mask[:, -1:] = 0

    else:
        mask = mask.astype(bool, copy=False)
        masked_image = np.zeros_like(image)
        masked_image[mask] = image[mask]

        # Make the eroded mask. Setting the border value to zero will wipe
        # out the image edges for us.
        s = ndi.generate_binary_structure(2, 2)
        eroded_mask = ndi.binary_erosion(mask, s, border_value=0)

    if compute_bleedover:
        # Compute the fractional contribution of masked pixels by applying
        # the function to the mask (which gets you the fraction of the
        # pixel data that's due to significant points)
        bleed_over = (
            gaussian(mask.astype(float_type, copy=False), **gaussian_kwargs)
            + np.finfo(float_type).eps
        )

    # Smooth the masked image
    smoothed_image = gaussian(masked_image, **gaussian_kwargs)

    # Lower the result by the bleed-over fraction, so you can
    # recalibrate by dividing by the function on the mask to recover
    # the effect of smoothing from just the significant pixels.
    if compute_bleedover:
        smoothed_image /= bleed_over

    return smoothed_image, eroded_mask


def canny(
    image,
    sigma=1.0,
    low_threshold=None,
    high_threshold=None,
    mask=None,
    use_quantiles=False,
    *,
    mode='constant',
    cval=0.0,
):
    """Edge filter an image using the Canny algorithm.

    Parameters
    ----------
    image : 2D array
        Grayscale input image to detect edges on; can be of any dtype.
    sigma : float, optional
        Standard deviation of the Gaussian filter.
    low_threshold : float, optional
        Lower bound for hysteresis thresholding (linking edges).
        If None, low_threshold is set to 10% of dtype's max.
    high_threshold : float, optional
        Upper bound for hysteresis thresholding (linking edges).
        If None, high_threshold is set to 20% of dtype's max.
    mask : array, dtype=bool, optional
        Mask to limit the application of Canny to a certain area.
    use_quantiles : bool, optional
        If ``True`` then treat low_threshold and high_threshold as
        quantiles of the edge magnitude image, rather than absolute
        edge magnitude values. If ``True`` then the thresholds must be
        in the range [0, 1].
    mode : str, {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}
        The ``mode`` parameter determines how the array borders are
        handled during Gaussian filtering, where ``cval`` is the value when
        mode is equal to 'constant'.
    cval : float, optional
        Value to fill past edges of input if `mode` is 'constant'.

    Returns
    -------
    output : 2D array (image)
        The binary edge map.

    See also
    --------
    skimage.filters.sobel

    Notes
    -----
    The steps of the algorithm are as follows:

    * Smooth the image using a Gaussian with ``sigma`` width.

    * Apply the horizontal and vertical Sobel operators to get the gradients
      within the image. The edge strength is the norm of the gradient.

    * Thin potential edges to 1-pixel wide curves. First, find the normal
      to the edge at each point. This is done by looking at the
      signs and the relative magnitude of the X-Sobel and Y-Sobel
      to sort the points into 4 categories: horizontal, vertical,
      diagonal and antidiagonal. Then look in the normal and reverse
      directions to see if the values in either of those directions are
      greater than the point in question. Use interpolation to get a mix of
      points instead of picking the one that's the closest to the normal.

    * Perform a hysteresis thresholding: first label all points above the
      high threshold as edges. Then recursively label any point above the
      low threshold that is 8-connected to a labeled point as an edge.

    References
    ----------
    .. [1] Canny, J., A Computational Approach To Edge Detection, IEEE Trans.
           Pattern Analysis and Machine Intelligence, 8:679-714, 1986
           :DOI:`10.1109/TPAMI.1986.4767851`
    .. [2] William Green's Canny tutorial
           https://en.wikipedia.org/wiki/Canny_edge_detector

    Examples
    --------
    >>> from skimage import feature
    >>> rng = np.random.default_rng()
    >>> # Generate noisy image of a square
    >>> im = np.zeros((256, 256))
    >>> im[64:-64, 64:-64] = 1
    >>> im += 0.2 * rng.random(im.shape)
    >>> # First trial with the Canny filter, with the default smoothing
    >>> edges1 = feature.canny(im)
    >>> # Increase the smoothing for better results
    >>> edges2 = feature.canny(im, sigma=3)

    """

    # Regarding masks, any point touching a masked point will have a gradient
    # that is "infected" by the masked point, so it's enough to erode the
    # mask by one and then mask the output. We also mask out the border points
    # because who knows what lies beyond the edge of the image?

    if np.issubdtype(image.dtype, np.int64) or np.issubdtype(image.dtype, np.uint64):
        raise ValueError("64-bit integer images are not supported")

    check_nD(image, 2)
    dtype_max = dtype_limits(image, clip_negative=False)[1]

    if low_threshold is None:
        low_threshold = 0.1
    elif use_quantiles:
        if not (0.0 <= low_threshold <= 1.0):
            raise ValueError("Quantile thresholds must be between 0 and 1.")
    else:
        low_threshold /= dtype_max

    if high_threshold is None:
        high_threshold = 0.2
    elif use_quantiles:
        if not (0.0 <= high_threshold <= 1.0):
            raise ValueError("Quantile thresholds must be between 0 and 1.")
    else:
        high_threshold /= dtype_max

    if high_threshold < low_threshold:
        raise ValueError("low_threshold should be lower then high_threshold")

    # Image filtering
    smoothed, eroded_mask = _preprocess(image, mask, sigma, mode, cval)

    # Gradient magnitude estimation
    jsobel = ndi.sobel(smoothed, axis=1)
    isobel = ndi.sobel(smoothed, axis=0)
    magnitude = isobel * isobel
    magnitude += jsobel * jsobel
    np.sqrt(magnitude, out=magnitude)

    if use_quantiles:
        low_threshold, high_threshold = np.percentile(
            magnitude, [100.0 * low_threshold, 100.0 * high_threshold]
        )

    # Non-maximum suppression
    low_masked = _nonmaximum_suppression_bilinear(
        isobel, jsobel, magnitude, eroded_mask, low_threshold
    )

    # Double thresholding and edge tracking
    #
    # Segment the low-mask, then only keep low-segments that have
    # some high_mask component in them
    #
    low_mask = low_masked > 0
    strel = np.ones((3, 3), bool)
    labels, count = ndi.label(low_mask, strel)
    if count == 0:
        return low_mask

    high_mask = low_mask & (low_masked >= high_threshold)
    nonzero_sums = np.unique(labels[high_mask])
    good_label = np.zeros((count + 1,), bool)
    good_label[nonzero_sums] = True
    output_mask = good_label[labels]
    return output_mask
