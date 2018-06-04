import itertools
import math
import numpy as np
from scipy import ndimage as ndi
from collections import OrderedDict
from ..exposure import histogram
from .._shared.utils import assert_nD, warn, deprecated
from ..transform import integral_image
from ..util import crop, dtype_limits


__all__ = ['try_all_threshold',
           'threshold_adaptive',
           'threshold_otsu',
           'threshold_yen',
           'threshold_isodata',
           'threshold_li',
           'threshold_local',
           'threshold_minimum',
           'threshold_mean',
           'threshold_niblack',
           'threshold_sauvola',
           'threshold_triangle',
           'apply_hysteresis_threshold']


def _try_all(image, methods=None, figsize=None, num_cols=2, verbose=True):
    """Returns a figure comparing the outputs of different methods.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.
    methods : dict, optional
        Names and associated functions.
        Functions must take and return an image.
    figsize : tuple, optional
        Figure size (in inches).
    num_cols : int, optional
        Number of columns.
    verbose : bool, optional
        Print function name for each method.

    Returns
    -------
    fig, ax : tuple
        Matplotlib figure and axes.
    """
    from matplotlib import pyplot as plt

    num_rows = math.ceil((len(methods) + 1.) / num_cols)
    num_rows = int(num_rows)  # Python 2.7 support
    fig, ax = plt.subplots(num_rows, num_cols, figsize=figsize,
                           sharex=True, sharey=True)
    ax = ax.ravel()

    ax[0].imshow(image, cmap=plt.cm.gray)
    ax[0].set_title('Original')

    i = 1
    for name, func in methods.items():
        ax[i].imshow(func(image), cmap=plt.cm.gray)
        ax[i].set_title(name)
        i += 1
        if verbose:
            print(func.__orifunc__)

    for a in ax:
        a.axis('off')

    fig.tight_layout()
    return fig, ax


def try_all_threshold(image, figsize=(8, 5), verbose=True):
    """Returns a figure comparing the outputs of different thresholding methods.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.
    figsize : tuple, optional
        Figure size (in inches).
    verbose : bool, optional
        Print function name for each method.

    Returns
    -------
    fig, ax : tuple
        Matplotlib figure and axes.

    Notes
    -----
    The following algorithms are used:

    * isodata
    * li
    * mean
    * minimum
    * otsu
    * triangle
    * yen

    Examples
    --------
    >>> from skimage.data import text
    >>> fig, ax = try_all_threshold(text(), figsize=(10, 6), verbose=False)
    """
    def thresh(func):
        """
        A wrapper function to return a thresholded image.
        """
        def wrapper(im):
            return im > func(im)
        try:
            wrapper.__orifunc__ = func.__orifunc__
        except AttributeError:
            wrapper.__orifunc__ = func.__module__ + '.' + func.__name__
        return wrapper

    # Global algorithms.
    methods = OrderedDict({'Isodata': thresh(threshold_isodata),
                           'Li': thresh(threshold_li),
                           'Mean': thresh(threshold_mean),
                           'Minimum': thresh(threshold_minimum),
                           'Otsu': thresh(threshold_otsu),
                           'Triangle': thresh(threshold_triangle),
                           'Yen': thresh(threshold_yen)})

    return _try_all(image, figsize=figsize,
                    methods=methods, verbose=verbose)


def threshold_local(image, block_size, method='gaussian', offset=0,
                    mode='reflect', param=None):
    """Compute a threshold mask image based on local pixel neighborhood.

    Also known as adaptive or dynamic thresholding. The threshold value is
    the weighted mean for the local neighborhood of a pixel subtracted by a
    constant. Alternatively the threshold can be determined dynamically by a
    given function, using the 'generic' method.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.
    block_size : int
        Odd size of pixel neighborhood which is used to calculate the
        threshold value (e.g. 3, 5, 7, ..., 21, ...).
    method : {'generic', 'gaussian', 'mean', 'median'}, optional
        Method used to determine adaptive threshold for local neighbourhood in
        weighted mean image.

        * 'generic': use custom function (see `param` parameter)
        * 'gaussian': apply gaussian filter (see `param` parameter for custom\
                      sigma value)
        * 'mean': apply arithmetic mean filter
        * 'median': apply median rank filter

        By default the 'gaussian' method is used.
    offset : float, optional
        Constant subtracted from weighted mean of neighborhood to calculate
        the local threshold value. Default offset is 0.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
        Default is 'reflect'.
    param : {int, function}, optional
        Either specify sigma for 'gaussian' method or function object for
        'generic' method. This functions takes the flat array of local
        neighbourhood as a single argument and returns the calculated
        threshold for the centre pixel.

    Returns
    -------
    threshold : (N, M) ndarray
        Threshold image. All pixels in the input image higher than the
        corresponding pixel in the threshold image are considered foreground.

    References
    ----------
    .. [1] http://docs.opencv.org/modules/imgproc/doc/miscellaneous_transformations.html?highlight=threshold#adaptivethreshold

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()[:50, :50]
    >>> binary_image1 = image > threshold_local(image, 15, 'mean')
    >>> func = lambda arr: arr.mean()
    >>> binary_image2 = image > threshold_local(image, 15, 'generic',
    ...                                         param=func)
    """
    if block_size % 2 == 0:
        raise ValueError("The kwarg ``block_size`` must be odd! Given "
                         "``block_size`` {0} is even.".format(block_size))
    assert_nD(image, 2)
    thresh_image = np.zeros(image.shape, 'double')
    if method == 'generic':
        ndi.generic_filter(image, param, block_size,
                           output=thresh_image, mode=mode)
    elif method == 'gaussian':
        if param is None:
            # automatically determine sigma which covers > 99% of distribution
            sigma = (block_size - 1) / 6.0
        else:
            sigma = param
        ndi.gaussian_filter(image, sigma, output=thresh_image, mode=mode)
    elif method == 'mean':
        mask = 1. / block_size * np.ones((block_size,))
        # separation of filters to speedup convolution
        ndi.convolve1d(image, mask, axis=0, output=thresh_image, mode=mode)
        ndi.convolve1d(thresh_image, mask, axis=1,
                       output=thresh_image, mode=mode)
    elif method == 'median':
        ndi.median_filter(image, block_size, output=thresh_image, mode=mode)

    return thresh_image - offset


@deprecated('threshold_local', removed_version='0.15')
def threshold_adaptive(image, block_size, method='gaussian', offset=0,
                       mode='reflect', param=None):
    warn('The return value of `threshold_local` is a threshold image, while '
         '`threshold_adaptive` returned the *thresholded* image.')
    return image > threshold_local(image, block_size=block_size,
                                   method=method, offset=offset, mode=mode,
                                   param=param)


def threshold_otsu(image, nbins=256):
    """Return threshold value based on Otsu's method.

    Parameters
    ----------
    image : (N, M) ndarray
        Grayscale input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    Raises
    ------
    ValueError
         If `image` only contains a single grayscale value.

    References
    ----------
    .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_otsu(image)
    >>> binary = image <= thresh

    Notes
    -----
    The input image must be grayscale.
    """
    if len(image.shape) > 2 and image.shape[-1] in (3, 4):
        msg = "threshold_otsu is expected to work correctly only for " \
              "grayscale images; image shape {0} looks like an RGB image"
        warn(msg.format(image.shape))

    # Check if the image is multi-colored or not
    if image.min() == image.max():
        raise ValueError("threshold_otsu is expected to work with images "
                         "having more than one color. The input image seems "
                         "to have just one color {0}.".format(image.min()))

    hist, bin_centers = histogram(image.ravel(), nbins)
    hist = hist.astype(float)

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]
    return threshold


def threshold_yen(image, nbins=256):
    """Return threshold value based on Yen's method.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    References
    ----------
    .. [1] Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion
           for Automatic Multilevel Thresholding" IEEE Trans. on Image
           Processing, 4(3): 370-378. DOI:10.1109/83.366472
    .. [2] Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding
           Techniques and Quantitative Performance Evaluation" Journal of
           Electronic Imaging, 13(1): 146-165, DOI:10.1117/1.1631315
           http://www.busim.ee.boun.edu.tr/~sankur/SankurFolder/Threshold_survey.pdf
    .. [3] ImageJ AutoThresholder code, http://fiji.sc/wiki/index.php/Auto_Threshold

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_yen(image)
    >>> binary = image <= thresh
    """
    hist, bin_centers = histogram(image.ravel(), nbins)
    # On blank images (e.g. filled with 0) with int dtype, `histogram()`
    # returns `bin_centers` containing only one value. Speed up with it.
    if bin_centers.size == 1:
        return bin_centers[0]

    # Calculate probability mass function
    pmf = hist.astype(np.float32) / hist.sum()
    P1 = np.cumsum(pmf)  # Cumulative normalized histogram
    P1_sq = np.cumsum(pmf ** 2)
    # Get cumsum calculated from end of squared array:
    P2_sq = np.cumsum(pmf[::-1] ** 2)[::-1]
    # P2_sq indexes is shifted +1. I assume, with P1[:-1] it's help avoid '-inf'
    # in crit. ImageJ Yen implementation replaces those values by zero.
    crit = np.log(((P1_sq[:-1] * P2_sq[1:]) ** -1) *
                  (P1[:-1] * (1.0 - P1[:-1])) ** 2)
    return bin_centers[crit.argmax()]


def threshold_isodata(image, nbins=256, return_all=False):
    """Return threshold value(s) based on ISODATA method.

    Histogram-based threshold, known as Ridler-Calvard method or inter-means.
    Threshold values returned satisfy the following equality::

        threshold = (image[image <= threshold].mean() +
                     image[image > threshold].mean()) / 2.0

    That is, returned thresholds are intensities that separate the image into
    two groups of pixels, where the threshold intensity is midway between the
    mean intensities of these groups.

    For integer images, the above equality holds to within one; for floating-
    point images, the equality holds to within the histogram bin-width.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.
    return_all: bool, optional
        If False (default), return only the lowest threshold that satisfies
        the above equality. If True, return all valid thresholds.

    Returns
    -------
    threshold : float or int or array
        Threshold value(s).

    References
    ----------
    .. [1] Ridler, TW & Calvard, S (1978), "Picture thresholding using an
           iterative selection method"
           IEEE Transactions on Systems, Man and Cybernetics 8: 630-632,
           DOI:10.1109/TSMC.1978.4310039
    .. [2] Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding
           Techniques and Quantitative Performance Evaluation" Journal of
           Electronic Imaging, 13(1): 146-165,
           http://www.busim.ee.boun.edu.tr/~sankur/SankurFolder/Threshold_survey.pdf
           DOI:10.1117/1.1631315
    .. [3] ImageJ AutoThresholder code,
           http://fiji.sc/wiki/index.php/Auto_Threshold

    Examples
    --------
    >>> from skimage.data import coins
    >>> image = coins()
    >>> thresh = threshold_isodata(image)
    >>> binary = image > thresh
    """
    hist, bin_centers = histogram(image.ravel(), nbins)

    # image only contains one unique value
    if len(bin_centers) == 1:
        if return_all:
            return bin_centers
        else:
            return bin_centers[0]

    hist = hist.astype(np.float32)

    # csuml and csumh contain the count of pixels in that bin or lower, and
    # in all bins strictly higher than that bin, respectively
    csuml = np.cumsum(hist)
    csumh = np.cumsum(hist[::-1])[::-1] - hist

    # intensity_sum contains the total pixel intensity from each bin
    intensity_sum = hist * bin_centers

    # l and h contain average value of all pixels in that bin or lower, and
    # in all bins strictly higher than that bin, respectively.
    # Note that since exp.histogram does not include empty bins at the low or
    # high end of the range, csuml and csumh are strictly > 0, except in the
    # last bin of csumh, which is zero by construction.
    # So no worries about division by zero in the following lines, except
    # for the last bin, but we can ignore that because no valid threshold
    # can be in the top bin. So we just patch up csumh[-1] to not cause 0/0
    # errors.
    csumh[-1] = 1
    l = np.cumsum(intensity_sum) / csuml
    h = (np.cumsum(intensity_sum[::-1])[::-1] - intensity_sum) / csumh

    # isodata finds threshold values that meet the criterion t = (l + m)/2
    # where l is the mean of all pixels <= t and h is the mean of all pixels
    # > t, as calculated above. So we are looking for places where
    # (l + m) / 2 equals the intensity value for which those l and m figures
    # were calculated -- which is, of course, the histogram bin centers.
    # We only require this equality to be within the precision of the bin
    # width, of course.
    all_mean = (l + h) / 2.0
    bin_width = bin_centers[1] - bin_centers[0]

    # Look only at thresholds that are below the actual all_mean value,
    # for consistency with the threshold being included in the lower pixel
    # group. Otherwise can get thresholds that are not actually fixed-points
    # of the isodata algorithm. For float images, this matters less, since
    # there really can't be any guarantees anymore anyway.
    distances = all_mean - bin_centers
    thresholds = bin_centers[(distances >= 0) & (distances < bin_width)]

    if return_all:
        return thresholds
    else:
        return thresholds[0]


def threshold_li(image):
    """Return threshold value based on adaptation of Li's Minimum Cross Entropy method.

    Parameters
    ----------
    image : (N, M) ndarray
        Input image.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    References
    ----------
    .. [1] Li C.H. and Lee C.K. (1993) "Minimum Cross Entropy Thresholding"
           Pattern Recognition, 26(4): 617-625
           DOI:10.1016/0031-3203(93)90115-D
    .. [2] Li C.H. and Tam P.K.S. (1998) "An Iterative Algorithm for Minimum
           Cross Entropy Thresholding" Pattern Recognition Letters, 18(8): 771-776
           DOI:10.1016/S0167-8655(98)00057-9
    .. [3] Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding
           Techniques and Quantitative Performance Evaluation" Journal of
           Electronic Imaging, 13(1): 146-165
           DOI:10.1117/1.1631315
    .. [4] ImageJ AutoThresholder code, http://fiji.sc/wiki/index.php/Auto_Threshold

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_li(image)
    >>> binary = image > thresh
    """
    # Make sure image has more than one value
    if np.all(image == image.flat[0]):
        raise ValueError("threshold_li is expected to work with images "
                         "having more than one value. The input image seems "
                         "to have just one value {0}.".format(image.flat[0]))

    # Copy to ensure input image is not modified
    image = image.copy()
    # Requires positive image (because of log(mean))
    immin = np.min(image)
    image -= immin
    imrange = np.max(image)
    tolerance = 0.5 * imrange / 256

    # Calculate the mean gray-level
    mean = np.mean(image)

    # Initial estimate
    new_thresh = mean
    old_thresh = new_thresh + 2 * tolerance

    # Stop the iterations when the difference between the
    # new and old threshold values is less than the tolerance
    while abs(new_thresh - old_thresh) > tolerance:
        old_thresh = new_thresh
        threshold = old_thresh + tolerance   # range
        # Calculate the means of background and object pixels
        mean_back = image[image <= threshold].mean()
        mean_obj = image[image > threshold].mean()

        temp = (mean_back - mean_obj) / (np.log(mean_back) - np.log(mean_obj))

        if temp < 0:
            new_thresh = temp - tolerance
        else:
            new_thresh = temp + tolerance

    return threshold + immin


def threshold_minimum(image, nbins=256, max_iter=10000):
    """Return threshold value based on minimum method.

    The histogram of the input `image` is computed and smoothed until there are
    only two maxima. Then the minimum in between is the threshold value.

    Parameters
    ----------
    image : (M, N) ndarray
        Input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.
    max_iter: int, optional
        Maximum number of iterations to smooth the histogram.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    Raises
    ------
    RuntimeError
        If unable to find two local maxima in the histogram or if the
        smoothing takes more than 1e4 iterations.

    References
    ----------
    .. [1] C. A. Glasbey, "An analysis of histogram-based thresholding
           algorithms," CVGIP: Graphical Models and Image Processing,
           vol. 55, pp. 532-537, 1993.
    .. [2] Prewitt, JMS & Mendelsohn, ML (1966), "The analysis of cell
           images", Annals of the New York Academy of Sciences 128: 1035-1053
           DOI:10.1111/j.1749-6632.1965.tb11715.x

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_minimum(image)
    >>> binary = image > thresh
    """

    def find_local_maxima_idx(hist):
        # We can't use scipy.signal.argrelmax
        # as it fails on plateaus
        maximum_idxs = list()
        direction = 1

        for i in range(hist.shape[0] - 1):
            if direction > 0:
                if hist[i + 1] < hist[i]:
                    direction = -1
                    maximum_idxs.append(i)
            else:
                if hist[i + 1] > hist[i]:
                    direction = 1

        return maximum_idxs

    hist, bin_centers = histogram(image.ravel(), nbins)

    smooth_hist = np.copy(hist).astype(np.float64)

    for counter in range(max_iter):
        smooth_hist = ndi.uniform_filter1d(smooth_hist, 3)
        maximum_idxs = find_local_maxima_idx(smooth_hist)
        if len(maximum_idxs) < 3:
            break

    if len(maximum_idxs) != 2:
        raise RuntimeError('Unable to find two maxima in histogram')
    elif counter == max_iter - 1:
        raise RuntimeError('Maximum iteration reached for histogram'
                           'smoothing')

    # Find lowest point between the maxima
    threshold_idx = np.argmin(smooth_hist[maximum_idxs[0]:maximum_idxs[1] + 1])

    return bin_centers[maximum_idxs[0] + threshold_idx]


def threshold_mean(image):
    """Return threshold value based on the mean of grayscale values.

    Parameters
    ----------
    image : (N, M[, ..., P]) ndarray
        Grayscale input image.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    References
    ----------
    .. [1] C. A. Glasbey, "An analysis of histogram-based thresholding
        algorithms," CVGIP: Graphical Models and Image Processing,
        vol. 55, pp. 532-537, 1993.
        DOI:10.1006/cgip.1993.1040

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_mean(image)
    >>> binary = image > thresh
    """
    return np.mean(image)


def threshold_triangle(image, nbins=256):
    """Return threshold value based on the triangle algorithm.

    Parameters
    ----------
    image : (N, M[, ..., P]) ndarray
        Grayscale input image.
    nbins : int, optional
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    References
    ----------
    .. [1] Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
       Automatic Measurement of Sister Chromatid Exchange Frequency,
       Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
       DOI:10.1177/25.7.70454
    .. [2] ImageJ AutoThresholder code,
       http://fiji.sc/wiki/index.php/Auto_Threshold

    Examples
    --------
    >>> from skimage.data import camera
    >>> image = camera()
    >>> thresh = threshold_triangle(image)
    >>> binary = image > thresh
    """
    # nbins is ignored for integer arrays
    # so, we recalculate the effective nbins.
    hist, bin_centers = histogram(image.ravel(), nbins)
    nbins = len(hist)

    # Find peak, lowest and highest gray levels.
    arg_peak_height = np.argmax(hist)
    peak_height = hist[arg_peak_height]
    arg_low_level, arg_high_level = np.where(hist>0)[0][[0, -1]]

    # Flip is True if left tail is shorter.
    flip = arg_peak_height - arg_low_level < arg_high_level - arg_peak_height
    if flip:
        hist = hist[::-1]
        arg_low_level = nbins - arg_high_level - 1
        arg_peak_height = nbins - arg_peak_height - 1

    # If flip == True, arg_high_level becomes incorrect
    # but we don't need it anymore.
    del(arg_high_level)

    # Set up the coordinate system.
    width = arg_peak_height - arg_low_level
    x1 = np.arange(width)
    y1 = hist[x1 + arg_low_level]

    # Normalize.
    norm = np.sqrt(peak_height**2 + width**2)
    peak_height /= norm
    width /= norm

    # Maximize the length.
    # The ImageJ implementation includes an additional constant when calculating
    # the length, but here we omit it as it does not affect the location of the
    # minimum.
    length = peak_height * x1 - width * y1
    arg_level = np.argmax(length) + arg_low_level

    if flip:
        arg_level = nbins - arg_level - 1

    return bin_centers[arg_level]


def _mean_std(image, w):
    """Return local mean and standard deviation of each pixel using a
    neighborhood defined by a rectangular window with size w times w.
    The algorithm uses integral images to speedup computation. This is
    used by threshold_niblack and threshold_sauvola.

    Parameters
    ----------
    image : ndarray
        Input image.
    w : int
        Odd window size (e.g. 3, 5, 7, ..., 21, ...).

    Returns
    -------
    m : 2-D array of same size of image with local mean values.
    s : 2-D array of same size of image with local standard
        deviation values.

    References
    ----------
    .. [1] F. Shafait, D. Keysers, and T. M. Breuel, "Efficient
           implementation of local adaptive thresholding techniques
           using integral images." in Document Recognition and
           Retrieval XV, (San Jose, USA), Jan. 2008.
           DOI:10.1117/12.767755
    """
    if w == 1 or w % 2 == 0:
        raise ValueError(
            "Window size w = %s must be odd and greater than 1." % w)

    left_pad = w // 2 + 1
    right_pad = w // 2
    padded = np.pad(image.astype('float'), (left_pad, right_pad),
                    mode='reflect')
    padded_sq = padded * padded

    integral = integral_image(padded)
    integral_sq = integral_image(padded_sq)

    kern = np.zeros((w + 1,) * image.ndim)
    for indices in itertools.product(*([[0, -1]] * image.ndim)):
        kern[indices] = (-1) ** (image.ndim % 2 != np.sum(indices) % 2)

    sum_full = ndi.correlate(integral, kern, mode='constant')
    m = crop(sum_full, (left_pad, right_pad)) / (w ** image.ndim)
    sum_sq_full = ndi.correlate(integral_sq, kern, mode='constant')
    g2 = crop(sum_sq_full, (left_pad, right_pad)) / (w ** image.ndim)
    # Note: we use np.clip because g2 is not guaranteed to be greater than
    # m*m when floating point error is considered
    s = np.sqrt(np.clip(g2 - m * m, 0, None))
    return m, s


def threshold_niblack(image, window_size=15, k=0.2):
    """Applies Niblack local threshold to an array.

    A threshold T is calculated for every pixel in the image using the
    following formula::

        T = m(x,y) - k * s(x,y)

    where m(x,y) and s(x,y) are the mean and standard deviation of
    pixel (x,y) neighborhood defined by a rectangular window with size w
    times w centered around the pixel. k is a configurable parameter
    that weights the effect of standard deviation.

    Parameters
    ----------
    image: (N, M) ndarray
        Grayscale input image.
    window_size : int, optional
        Odd size of pixel neighborhood window (e.g. 3, 5, 7...).
    k : float, optional
        Value of parameter k in threshold formula.

    Returns
    -------
    threshold : (N, M) ndarray
        Threshold mask. All pixels with an intensity higher than
        this value are assumed to be foreground.

    Notes
    -----
    This algorithm is originally designed for text recognition.

    References
    ----------
    .. [1] Niblack, W (1986), An introduction to Digital Image
           Processing, Prentice-Hall.

    Examples
    --------
    >>> from skimage import data
    >>> image = data.page()
    >>> binary_image = threshold_niblack(image, window_size=7, k=0.1)
    """
    m, s = _mean_std(image, window_size)
    return m - k * s


def threshold_sauvola(image, window_size=15, k=0.2, r=None):
    """Applies Sauvola local threshold to an array. Sauvola is a
    modification of Niblack technique.

    In the original method a threshold T is calculated for every pixel
    in the image using the following formula::

        T = m(x,y) * (1 + k * ((s(x,y) / R) - 1))

    where m(x,y) and s(x,y) are the mean and standard deviation of
    pixel (x,y) neighborhood defined by a rectangular window with size w
    times w centered around the pixel. k is a configurable parameter
    that weights the effect of standard deviation.
    R is the maximum standard deviation of a greyscale image.

    Parameters
    ----------
    image: (N, M) ndarray
        Grayscale input image.
    window_size : int, optional
        Odd size of pixel neighborhood window (e.g. 3, 5, 7...).
    k : float, optional
        Value of the positive parameter k.
    r : float, optional
        Value of R, the dynamic range of standard deviation.
        If None, set to the half of the image dtype range.

    Returns
    -------
    threshold : (N, M) ndarray
        Threshold mask. All pixels with an intensity higher than
        this value are assumed to be foreground.

    Notes
    -----
    This algorithm is originally designed for text recognition.

    References
    ----------
    .. [1] J. Sauvola and M. Pietikainen, "Adaptive document image
           binarization," Pattern Recognition 33(2),
           pp. 225-236, 2000.
           DOI:10.1016/S0031-3203(99)00055-2

    Examples
    --------
    >>> from skimage import data
    >>> image = data.page()
    >>> t_sauvola = threshold_sauvola(image, window_size=15, k=0.2)
    >>> binary_image = image > t_sauvola
    """
    if r is None:
        imin, imax = dtype_limits(image, clip_negative=False)
        r = 0.5 * (imax - imin)
    m, s = _mean_std(image, window_size)
    return m * (1 + k * ((s / r) - 1))


def apply_hysteresis_threshold(image, low, high):
    """Apply hysteresis thresholding to `image`.

    This algorithm finds regions where `image` is greater than `high`
    OR `image` is greater than `low` *and* that region is connected to
    a region greater than `high`.

    Parameters
    ----------
    image : array, shape (M,[ N, ..., P])
        Grayscale input image.
    low : float, or array of same shape as `image`
        Lower threshold.
    high : float, or array of same shape as `image`
        Higher threshold.

    Returns
    -------
    thresholded : array of bool, same shape as `image`
        Array in which `True` indicates the locations where `image`
        was above the hysteresis threshold.

    Examples
    --------
    >>> image = np.array([1, 2, 3, 2, 1, 2, 1, 3, 2])
    >>> apply_hysteresis_threshold(image, 1.5, 2.5).astype(int)
    array([0, 1, 1, 1, 0, 0, 0, 1, 1])

    References
    ----------
    .. [1] J. Canny. A computational approach to edge detection.
           IEEE Transactions on Pattern Analysis and Machine Intelligence.
           1986; vol. 8, pp.679-698.
           DOI: 10.1109/TPAMI.1986.4767851
    """
    low = np.clip(low, a_min=None, a_max=high)  # ensure low always below high
    mask_low = image > low
    mask_high = image > high
    # Connected components of mask_low
    labels_low, num_labels = ndi.label(mask_low)
    # Check which connected components contain pixels from mask_high
    sums = ndi.sum(mask_high, labels_low, np.arange(num_labels + 1))
    connected_to_high = sums > 0
    thresholded = connected_to_high[labels_low]
    return thresholded
