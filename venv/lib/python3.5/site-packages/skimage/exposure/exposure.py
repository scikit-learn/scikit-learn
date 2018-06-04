from __future__ import division
import numpy as np

from ..color import rgb2gray
from ..util.dtype import dtype_range, dtype_limits
from .._shared.utils import warn


__all__ = ['histogram', 'cumulative_distribution', 'equalize_hist',
           'rescale_intensity', 'adjust_gamma', 'adjust_log', 'adjust_sigmoid']


DTYPE_RANGE = dtype_range.copy()
DTYPE_RANGE.update((d.__name__, limits) for d, limits in dtype_range.items())
DTYPE_RANGE.update({'uint10': (0, 2 ** 10 - 1),
                    'uint12': (0, 2 ** 12 - 1),
                    'uint14': (0, 2 ** 14 - 1),
                    'bool': dtype_range[np.bool_],
                    'float': dtype_range[np.float64]})


def histogram(image, nbins=256):
    """Return histogram of image.

    Unlike `numpy.histogram`, this function returns the centers of bins and
    does not rebin integer arrays. For integer arrays, each integer value has
    its own bin, which improves speed and intensity-resolution.

    The histogram is computed on the flattened image: for color images, the
    function should be used separately on each channel to obtain a histogram
    for each color channel.

    Parameters
    ----------
    image : array
        Input image.
    nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    hist : array
        The values of the histogram.
    bin_centers : array
        The values at the center of the bins.

    See Also
    --------
    cumulative_distribution

    Examples
    --------
    >>> from skimage import data, exposure, img_as_float
    >>> image = img_as_float(data.camera())
    >>> np.histogram(image, bins=2)
    (array([107432, 154712]), array([ 0. ,  0.5,  1. ]))
    >>> exposure.histogram(image, nbins=2)
    (array([107432, 154712]), array([ 0.25,  0.75]))
    """
    sh = image.shape
    if len(sh) == 3 and sh[-1] < 4:
        warn("This might be a color image. The histogram will be "
             "computed on the flattened image. You can instead "
             "apply this function to each color channel.")

    # For integer types, histogramming with bincount is more efficient.
    if np.issubdtype(image.dtype, np.integer):
        offset = 0
        image_min = np.min(image)
        if image_min < 0:
            offset = image_min
            image_range = np.max(image).astype(np.int64) - image_min
            # get smallest dtype that can hold both minimum and offset maximum
            offset_dtype = np.promote_types(np.min_scalar_type(image_range),
                                            np.min_scalar_type(image_min))
            if image.dtype != offset_dtype:
                # prevent overflow errors when offsetting
                image = image.astype(offset_dtype)
            image = image - offset
        hist = np.bincount(image.ravel())
        bin_centers = np.arange(len(hist)) + offset

        # clip histogram to start with a non-zero bin
        idx = np.nonzero(hist)[0][0]
        return hist[idx:], bin_centers[idx:]
    else:
        hist, bin_edges = np.histogram(image.flat, bins=nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
        return hist, bin_centers


def cumulative_distribution(image, nbins=256):
    """Return cumulative distribution function (cdf) for the given image.

    Parameters
    ----------
    image : array
        Image array.
    nbins : int
        Number of bins for image histogram.

    Returns
    -------
    img_cdf : array
        Values of cumulative distribution function.
    bin_centers : array
        Centers of bins.

    See Also
    --------
    histogram

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Cumulative_distribution_function

    Examples
    --------
    >>> from skimage import data, exposure, img_as_float
    >>> image = img_as_float(data.camera())
    >>> hi = exposure.histogram(image)
    >>> cdf = exposure.cumulative_distribution(image)
    >>> np.alltrue(cdf[0] == np.cumsum(hi[0])/float(image.size))
    True
    """
    hist, bin_centers = histogram(image, nbins)
    img_cdf = hist.cumsum()
    img_cdf = img_cdf / float(img_cdf[-1])
    return img_cdf, bin_centers


def equalize_hist(image, nbins=256, mask=None):
    """Return image after histogram equalization.

    Parameters
    ----------
    image : array
        Image array.
    nbins : int, optional
        Number of bins for image histogram. Note: this argument is
        ignored for integer images, for which each integer is its own
        bin.
    mask: ndarray of bools or 0s and 1s, optional
        Array of same shape as `image`. Only points at which mask == True
        are used for the equalization, which is applied to the whole image.

    Returns
    -------
    out : float array
        Image array after histogram equalization.

    Notes
    -----
    This function is adapted from [1]_ with the author's permission.

    References
    ----------
    .. [1] http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html
    .. [2] http://en.wikipedia.org/wiki/Histogram_equalization

    """
    if mask is not None:
        mask = np.array(mask, dtype=bool)
        cdf, bin_centers = cumulative_distribution(image[mask], nbins)
    else:
        cdf, bin_centers = cumulative_distribution(image, nbins)
    out = np.interp(image.flat, bin_centers, cdf)
    return out.reshape(image.shape)


def intensity_range(image, range_values='image', clip_negative=False):
    """Return image intensity range (min, max) based on desired value type.

    Parameters
    ----------
    image : array
        Input image.
    range_values : str or 2-tuple
        The image intensity range is configured by this parameter.
        The possible values for this parameter are enumerated below.

        'image'
            Return image min/max as the range.
        'dtype'
            Return min/max of the image's dtype as the range.
        dtype-name
            Return intensity range based on desired `dtype`. Must be valid key
            in `DTYPE_RANGE`. Note: `image` is ignored for this range type.
        2-tuple
            Return `range_values` as min/max intensities. Note that there's no
            reason to use this function if you just want to specify the
            intensity range explicitly. This option is included for functions
            that use `intensity_range` to support all desired range types.

    clip_negative : bool
        If True, clip the negative range (i.e. return 0 for min intensity)
        even if the image dtype allows negative values.
    """
    if range_values == 'dtype':
        range_values = image.dtype.type

    if range_values == 'image':
        i_min = np.min(image)
        i_max = np.max(image)
    elif range_values in DTYPE_RANGE:
        i_min, i_max = DTYPE_RANGE[range_values]
        if clip_negative:
            i_min = 0
    else:
        i_min, i_max = range_values
    return i_min, i_max


def rescale_intensity(image, in_range='image', out_range='dtype'):
    """Return image after stretching or shrinking its intensity levels.

    The desired intensity range of the input and output, `in_range` and
    `out_range` respectively, are used to stretch or shrink the intensity range
    of the input image. See examples below.

    Parameters
    ----------
    image : array
        Image array.
    in_range, out_range : str or 2-tuple
        Min and max intensity values of input and output image.
        The possible values for this parameter are enumerated below.

        'image'
            Use image min/max as the intensity range.
        'dtype'
            Use min/max of the image's dtype as the intensity range.
        dtype-name
            Use intensity range based on desired `dtype`. Must be valid key
            in `DTYPE_RANGE`.
        2-tuple
            Use `range_values` as explicit min/max intensities.

    Returns
    -------
    out : array
        Image array after rescaling its intensity. This image is the same dtype
        as the input image.

    See Also
    --------
    equalize_hist

    Examples
    --------
    By default, the min/max intensities of the input image are stretched to
    the limits allowed by the image's dtype, since `in_range` defaults to
    'image' and `out_range` defaults to 'dtype':

    >>> image = np.array([51, 102, 153], dtype=np.uint8)
    >>> rescale_intensity(image)
    array([  0, 127, 255], dtype=uint8)

    It's easy to accidentally convert an image dtype from uint8 to float:

    >>> 1.0 * image
    array([  51.,  102.,  153.])

    Use `rescale_intensity` to rescale to the proper range for float dtypes:

    >>> image_float = 1.0 * image
    >>> rescale_intensity(image_float)
    array([ 0. ,  0.5,  1. ])

    To maintain the low contrast of the original, use the `in_range` parameter:

    >>> rescale_intensity(image_float, in_range=(0, 255))
    array([ 0.2,  0.4,  0.6])

    If the min/max value of `in_range` is more/less than the min/max image
    intensity, then the intensity levels are clipped:

    >>> rescale_intensity(image_float, in_range=(0, 102))
    array([ 0.5,  1. ,  1. ])

    If you have an image with signed integers but want to rescale the image to
    just the positive range, use the `out_range` parameter:

    >>> image = np.array([-10, 0, 10], dtype=np.int8)
    >>> rescale_intensity(image, out_range=(0, 127))
    array([  0,  63, 127], dtype=int8)

    """
    dtype = image.dtype.type

    imin, imax = intensity_range(image, in_range)
    omin, omax = intensity_range(image, out_range, clip_negative=(imin >= 0))

    image = np.clip(image, imin, imax)

    image = (image - imin) / float(imax - imin)
    return dtype(image * (omax - omin) + omin)


def _assert_non_negative(image):

    if np.any(image < 0):
        raise ValueError('Image Correction methods work correctly only on '
                         'images with non-negative values. Use '
                         'skimage.exposure.rescale_intensity.')


def adjust_gamma(image, gamma=1, gain=1):
    """Performs Gamma Correction on the input image.

    Also known as Power Law Transform.
    This function transforms the input image pixelwise according to the
    equation ``O = I**gamma`` after scaling each pixel to the range 0 to 1.

    Parameters
    ----------
    image : ndarray
        Input image.
    gamma : float
        Non negative real number. Default value is 1.
    gain : float
        The constant multiplier. Default value is 1.

    Returns
    -------
    out : ndarray
        Gamma corrected output image.

    See Also
    --------
    adjust_log

    Notes
    -----
    For gamma greater than 1, the histogram will shift towards left and
    the output image will be darker than the input image.

    For gamma less than 1, the histogram will shift towards right and
    the output image will be brighter than the input image.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Gamma_correction

    Examples
    --------
    >>> from skimage import data, exposure, img_as_float
    >>> image = img_as_float(data.moon())
    >>> gamma_corrected = exposure.adjust_gamma(image, 2)
    >>> # Output is darker for gamma > 1
    >>> image.mean() > gamma_corrected.mean()
    True
    """
    _assert_non_negative(image)
    dtype = image.dtype.type

    if gamma < 0:
        raise ValueError("Gamma should be a non-negative real number.")

    scale = float(dtype_limits(image, True)[1] - dtype_limits(image, True)[0])

    out = ((image / scale) ** gamma) * scale * gain
    return dtype(out)


def adjust_log(image, gain=1, inv=False):
    """Performs Logarithmic correction on the input image.

    This function transforms the input image pixelwise according to the
    equation ``O = gain*log(1 + I)`` after scaling each pixel to the range 0 to 1.
    For inverse logarithmic correction, the equation is ``O = gain*(2**I - 1)``.

    Parameters
    ----------
    image : ndarray
        Input image.
    gain : float
        The constant multiplier. Default value is 1.
    inv : float
        If True, it performs inverse logarithmic correction,
        else correction will be logarithmic. Defaults to False.

    Returns
    -------
    out : ndarray
        Logarithm corrected output image.

    See Also
    --------
    adjust_gamma

    References
    ----------
    .. [1] http://www.ece.ucsb.edu/Faculty/Manjunath/courses/ece178W03/EnhancePart1.pdf

    """
    _assert_non_negative(image)
    dtype = image.dtype.type
    scale = float(dtype_limits(image, True)[1] - dtype_limits(image, True)[0])

    if inv:
        out = (2 ** (image / scale) - 1) * scale * gain
        return dtype(out)

    out = np.log2(1 + image / scale) * scale * gain
    return dtype(out)


def adjust_sigmoid(image, cutoff=0.5, gain=10, inv=False):
    """Performs Sigmoid Correction on the input image.

    Also known as Contrast Adjustment.
    This function transforms the input image pixelwise according to the
    equation ``O = 1/(1 + exp*(gain*(cutoff - I)))`` after scaling each pixel
    to the range 0 to 1.

    Parameters
    ----------
    image : ndarray
        Input image.
    cutoff : float
        Cutoff of the sigmoid function that shifts the characteristic curve
        in horizontal direction. Default value is 0.5.
    gain : float
        The constant multiplier in exponential's power of sigmoid function.
        Default value is 10.
    inv : bool
        If True, returns the negative sigmoid correction. Defaults to False.

    Returns
    -------
    out : ndarray
        Sigmoid corrected output image.

    See Also
    --------
    adjust_gamma

    References
    ----------
    .. [1] Gustav J. Braun, "Image Lightness Rescaling Using Sigmoidal Contrast
           Enhancement Functions",
           http://www.cis.rit.edu/fairchild/PDFs/PAP07.pdf

    """
    _assert_non_negative(image)
    dtype = image.dtype.type
    scale = float(dtype_limits(image, True)[1] - dtype_limits(image, True)[0])

    if inv:
        out = (1 - 1 / (1 + np.exp(gain * (cutoff - image / scale)))) * scale
        return dtype(out)

    out = (1 / (1 + np.exp(gain * (cutoff - image / scale)))) * scale
    return dtype(out)


def is_low_contrast(image, fraction_threshold=0.05, lower_percentile=1,
                    upper_percentile=99, method='linear'):
    """Detemine if an image is low contrast.

    Parameters
    ----------
    image : array-like
        The image under test.
    fraction_threshold : float, optional
        The low contrast fraction threshold. An image is considered low-
        contrast when its range of brightness spans less than this
        fraction of its data type's full range. [1]_
    lower_percentile : float, optional
        Disregard values below this percentile when computing image contrast.
    upper_percentile : float, optional
        Disregard values above this percentile when computing image contrast.
    method : str, optional
        The contrast determination method.  Right now the only available
        option is "linear".

    Returns
    -------
    out : bool
        True when the image is determined to be low contrast.

    References
    ----------
    .. [1] http://scikit-image.org/docs/dev/user_guide/data_types.html

    Examples
    --------
    >>> image = np.linspace(0, 0.04, 100)
    >>> is_low_contrast(image)
    True
    >>> image[-1] = 1
    >>> is_low_contrast(image)
    True
    >>> is_low_contrast(image, upper_percentile=100)
    False
    """
    image = np.asanyarray(image)
    if image.ndim == 3 and image.shape[2] in [3, 4]:
        image = rgb2gray(image)

    dlimits = dtype_limits(image, clip_negative=False)
    limits = np.percentile(image, [lower_percentile, upper_percentile])
    ratio = (limits[1] - limits[0]) / (dlimits[1] - dlimits[0])

    return ratio < fraction_threshold
