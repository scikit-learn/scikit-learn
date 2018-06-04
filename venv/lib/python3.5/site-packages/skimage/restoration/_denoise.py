# coding: utf-8
import scipy.stats
import numpy as np
from math import ceil
from .. import img_as_float
from ..restoration._denoise_cy import _denoise_bilateral, _denoise_tv_bregman
from .._shared.utils import skimage_deprecation, warn
import pywt
import skimage.color as color
import numbers


def denoise_bilateral(image, win_size=None, sigma_color=None, sigma_spatial=1,
                      bins=10000, mode='constant', cval=0, multichannel=None):
    """Denoise image using bilateral filter.

    This is an edge-preserving, denoising filter. It averages pixels based on
    their spatial closeness and radiometric similarity [1]_.

    Spatial closeness is measured by the Gaussian function of the Euclidean
    distance between two pixels and a certain standard deviation
    (`sigma_spatial`).

    Radiometric similarity is measured by the Gaussian function of the
    Euclidean distance between two color values and a certain standard
    deviation (`sigma_color`).

    Parameters
    ----------
    image : ndarray, shape (M, N[, 3])
        Input image, 2D grayscale or RGB.
    win_size : int
        Window size for filtering.
        If win_size is not specified, it is calculated as
        ``max(5, 2 * ceil(3 * sigma_spatial) + 1)``.
    sigma_color : float
        Standard deviation for grayvalue/color distance (radiometric
        similarity). A larger value results in averaging of pixels with larger
        radiometric differences. Note, that the image will be converted using
        the `img_as_float` function and thus the standard deviation is in
        respect to the range ``[0, 1]``. If the value is ``None`` the standard
        deviation of the ``image`` will be used.
    sigma_spatial : float
        Standard deviation for range distance. A larger value results in
        averaging of pixels with larger spatial differences.
    bins : int
        Number of discrete values for Gaussian weights of color filtering.
        A larger value results in improved accuracy.
    mode : {'constant', 'edge', 'symmetric', 'reflect', 'wrap'}
        How to handle values outside the image borders. See
        `numpy.pad` for detail.
    cval : string
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.
    multichannel : bool
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension.

    Returns
    -------
    denoised : ndarray
        Denoised image.

    References
    ----------
    .. [1] http://users.soe.ucsc.edu/~manduchi/Papers/ICCV98.pdf

    Examples
    --------
    >>> from skimage import data, img_as_float
    >>> astro = img_as_float(data.astronaut())
    >>> astro = astro[220:300, 220:320]
    >>> noisy = astro + 0.6 * astro.std() * np.random.random(astro.shape)
    >>> noisy = np.clip(noisy, 0, 1)
    >>> denoised = denoise_bilateral(noisy, sigma_color=0.05, sigma_spatial=15)
    """
    if multichannel is None:
        warn('denoise_bilateral will default to multichannel=False in v0.15')
        multichannel = True

    if multichannel:
        if image.ndim != 3:
            if image.ndim == 2:
                raise ValueError("Use ``multichannel=False`` for 2D grayscale "
                                 "images. The last axis of the input image "
                                 "must be multiple color channels not another "
                                 "spatial dimension.")
            else:
                raise ValueError("Bilateral filter is only implemented for "
                                 "2D grayscale images (image.ndim == 2) and "
                                 "2D multichannel (image.ndim == 3) images, "
                                 "but the input image has {0} dimensions. "
                                 "".format(image.ndim))
        elif image.shape[2] not in (3, 4):
            if image.shape[2] > 4:
                msg = ("The last axis of the input image is interpreted as "
                       "channels. Input image with shape {0} has {1} channels "
                       "in last axis. ``denoise_bilateral`` is implemented "
                       "for 2D grayscale and color images only")
                warn(msg.format(image.shape, image.shape[2]))
            else:
                msg = "Input image must be grayscale, RGB, or RGBA; " \
                      "but has shape {0}."
                warn(msg.format(image.shape))
    else:
        if image.ndim > 2:
            raise ValueError("Bilateral filter is not implemented for "
                             "grayscale images of 3 or more dimensions, "
                             "but input image has {0} dimension. Use "
                             "``multichannel=True`` for 2-D RGB "
                             "images.".format(image.shape))

    if win_size is None:
        win_size = max(5, 2 * int(ceil(3 * sigma_spatial)) + 1)

    return _denoise_bilateral(image, win_size, sigma_color, sigma_spatial,
                              bins, mode, cval)


def denoise_tv_bregman(image, weight, max_iter=100, eps=1e-3, isotropic=True):
    """Perform total-variation denoising using split-Bregman optimization.

    Total-variation denoising (also know as total-variation regularization)
    tries to find an image with less total-variation under the constraint
    of being similar to the input image, which is controlled by the
    regularization parameter ([1]_, [2]_, [3]_, [4]_).

    Parameters
    ----------
    image : ndarray
        Input data to be denoised (converted using img_as_float`).
    weight : float
        Denoising weight. The smaller the `weight`, the more denoising (at
        the expense of less similarity to the `input`). The regularization
        parameter `lambda` is chosen as `2 * weight`.
    eps : float, optional
        Relative difference of the value of the cost function that determines
        the stop criterion. The algorithm stops when::

            SUM((u(n) - u(n-1))**2) < eps

    max_iter : int, optional
        Maximal number of iterations used for the optimization.
    isotropic : boolean, optional
        Switch between isotropic and anisotropic TV denoising.

    Returns
    -------
    u : ndarray
        Denoised image.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Total_variation_denoising
    .. [2] Tom Goldstein and Stanley Osher, "The Split Bregman Method For L1
           Regularized Problems",
           ftp://ftp.math.ucla.edu/pub/camreport/cam08-29.pdf
    .. [3] Pascal Getreuer, "Rudin–Osher–Fatemi Total Variation Denoising
           using Split Bregman" in Image Processing On Line on 2012–05–19,
           http://www.ipol.im/pub/art/2012/g-tvd/article_lr.pdf
    .. [4] http://www.math.ucsb.edu/~cgarcia/UGProjects/BregmanAlgorithms_JacquelineBush.pdf

    """
    return _denoise_tv_bregman(image, weight, max_iter, eps, isotropic)


def _denoise_tv_chambolle_nd(image, weight=0.1, eps=2.e-4, n_iter_max=200):
    """Perform total-variation denoising on n-dimensional images.

    Parameters
    ----------
    image : ndarray
        n-D input data to be denoised.
    weight : float, optional
        Denoising weight. The greater `weight`, the more denoising (at
        the expense of fidelity to `input`).
    eps : float, optional
        Relative difference of the value of the cost function that determines
        the stop criterion. The algorithm stops when:

            (E_(n-1) - E_n) < eps * E_0

    n_iter_max : int, optional
        Maximal number of iterations used for the optimization.

    Returns
    -------
    out : ndarray
        Denoised array of floats.

    Notes
    -----
    Rudin, Osher and Fatemi algorithm.

    """

    ndim = image.ndim
    p = np.zeros((image.ndim, ) + image.shape, dtype=image.dtype)
    g = np.zeros_like(p)
    d = np.zeros_like(image)
    i = 0
    while i < n_iter_max:
        if i > 0:
            # d will be the (negative) divergence of p
            d = -p.sum(0)
            slices_d = [slice(None), ] * ndim
            slices_p = [slice(None), ] * (ndim + 1)
            for ax in range(ndim):
                slices_d[ax] = slice(1, None)
                slices_p[ax+1] = slice(0, -1)
                slices_p[0] = ax
                d[slices_d] += p[slices_p]
                slices_d[ax] = slice(None)
                slices_p[ax+1] = slice(None)
            out = image + d
        else:
            out = image
        E = (d ** 2).sum()

        # g stores the gradients of out along each axis
        # e.g. g[0] is the first order finite difference along axis 0
        slices_g = [slice(None), ] * (ndim + 1)
        for ax in range(ndim):
            slices_g[ax+1] = slice(0, -1)
            slices_g[0] = ax
            g[slices_g] = np.diff(out, axis=ax)
            slices_g[ax+1] = slice(None)

        norm = np.sqrt((g ** 2).sum(axis=0))[np.newaxis, ...]
        E += weight * norm.sum()
        tau = 1. / (2.*ndim)
        norm *= tau / weight
        norm += 1.
        p -= tau * g
        p /= norm
        E /= float(image.size)
        if i == 0:
            E_init = E
            E_previous = E
        else:
            if np.abs(E_previous - E) < eps * E_init:
                break
            else:
                E_previous = E
        i += 1
    return out


def denoise_tv_chambolle(image, weight=0.1, eps=2.e-4, n_iter_max=200,
                         multichannel=False):
    """Perform total-variation denoising on n-dimensional images.

    Parameters
    ----------
    image : ndarray of ints, uints or floats
        Input data to be denoised. `image` can be of any numeric type,
        but it is cast into an ndarray of floats for the computation
        of the denoised image.
    weight : float, optional
        Denoising weight. The greater `weight`, the more denoising (at
        the expense of fidelity to `input`).
    eps : float, optional
        Relative difference of the value of the cost function that
        determines the stop criterion. The algorithm stops when:

            (E_(n-1) - E_n) < eps * E_0

    n_iter_max : int, optional
        Maximal number of iterations used for the optimization.
    multichannel : bool, optional
        Apply total-variation denoising separately for each channel. This
        option should be true for color images, otherwise the denoising is
        also applied in the channels dimension.

    Returns
    -------
    out : ndarray
        Denoised image.

    Notes
    -----
    Make sure to set the multichannel parameter appropriately for color images.

    The principle of total variation denoising is explained in
    http://en.wikipedia.org/wiki/Total_variation_denoising

    The principle of total variation denoising is to minimize the
    total variation of the image, which can be roughly described as
    the integral of the norm of the image gradient. Total variation
    denoising tends to produce "cartoon-like" images, that is,
    piecewise-constant images.

    This code is an implementation of the algorithm of Rudin, Fatemi and Osher
    that was proposed by Chambolle in [1]_.

    References
    ----------
    .. [1] A. Chambolle, An algorithm for total variation minimization and
           applications, Journal of Mathematical Imaging and Vision,
           Springer, 2004, 20, 89-97.

    Examples
    --------
    2D example on astronaut image:

    >>> from skimage import color, data
    >>> img = color.rgb2gray(data.astronaut())[:50, :50]
    >>> img += 0.5 * img.std() * np.random.randn(*img.shape)
    >>> denoised_img = denoise_tv_chambolle(img, weight=60)

    3D example on synthetic data:

    >>> x, y, z = np.ogrid[0:20, 0:20, 0:20]
    >>> mask = (x - 22)**2 + (y - 20)**2 + (z - 17)**2 < 8**2
    >>> mask = mask.astype(np.float)
    >>> mask += 0.2*np.random.randn(*mask.shape)
    >>> res = denoise_tv_chambolle(mask, weight=100)

    """

    im_type = image.dtype
    if not im_type.kind == 'f':
        image = img_as_float(image)

    if multichannel:
        out = np.zeros_like(image)
        for c in range(image.shape[-1]):
            out[..., c] = _denoise_tv_chambolle_nd(image[..., c], weight, eps,
                                                   n_iter_max)
    else:
        out = _denoise_tv_chambolle_nd(image, weight, eps, n_iter_max)
    return out


def _bayes_thresh(details, var):
    """BayesShrink threshold for a zero-mean details coeff array."""
    # Equivalent to:  dvar = np.var(details) for 0-mean details array
    dvar = np.mean(details*details)
    eps = np.finfo(details.dtype).eps
    thresh = var / np.sqrt(max(dvar - var, eps))
    return thresh


def _universal_thresh(img, sigma):
    """ Universal threshold used by the VisuShrink method """
    return sigma*np.sqrt(2*np.log(img.size))


def _sigma_est_dwt(detail_coeffs, distribution='Gaussian'):
    """Calculate the robust median estimator of the noise standard deviation.

    Parameters
    ----------
    detail_coeffs : ndarray
        The detail coefficients corresponding to the discrete wavelet
        transform of an image.
    distribution : str
        The underlying noise distribution.

    Returns
    -------
    sigma : float
        The estimated noise standard deviation (see section 4.2 of [1]_).

    References
    ----------
    .. [1] D. L. Donoho and I. M. Johnstone. "Ideal spatial adaptation
       by wavelet shrinkage." Biometrika 81.3 (1994): 425-455.
       DOI:10.1093/biomet/81.3.425
    """
    # Consider regions with detail coefficients exactly zero to be masked out
    detail_coeffs = detail_coeffs[np.nonzero(detail_coeffs)]

    if distribution.lower() == 'gaussian':
        # 75th quantile of the underlying, symmetric noise distribution
        denom = scipy.stats.norm.ppf(0.75)
        sigma = np.median(np.abs(detail_coeffs)) / denom
    else:
        raise ValueError("Only Gaussian noise estimation is currently "
                         "supported")
    return sigma


def _wavelet_threshold(image, wavelet, method=None, threshold=None,
                       sigma=None, mode='soft', wavelet_levels=None):
    """Perform wavelet thresholding.

    Parameters
    ----------
    image : ndarray (2d or 3d) of ints, uints or floats
        Input data to be denoised. `image` can be of any numeric type,
        but it is cast into an ndarray of floats for the computation
        of the denoised image.
    wavelet : string
        The type of wavelet to perform. Can be any of the options
        pywt.wavelist outputs. For example, this may be any of ``{db1, db2,
        db3, db4, haar}``.
    method : {'BayesShrink', 'VisuShrink'}, optional
        Thresholding method to be used. The currently supported methods are
        "BayesShrink" [1]_ and "VisuShrink" [2]_. If it is set to None, a
        user-specified ``threshold`` must be supplied instead.
    threshold : float, optional
        The thresholding value to apply during wavelet coefficient
        thresholding. The default value (None) uses the selected ``method`` to
        estimate appropriate threshold(s) for noise removal.
    sigma : float, optional
        The standard deviation of the noise. The noise is estimated when sigma
        is None (the default) by the method in [2]_.
    mode : {'soft', 'hard'}, optional
        An optional argument to choose the type of denoising performed. It
        noted that choosing soft thresholding given additive noise finds the
        best approximation of the original image.
    wavelet_levels : int or None, optional
        The number of wavelet decomposition levels to use.  The default is
        three less than the maximum number of possible decomposition levels
        (see Notes below).

    Returns
    -------
    out : ndarray
        Denoised image.

    References
    ----------
    .. [1] Chang, S. Grace, Bin Yu, and Martin Vetterli. "Adaptive wavelet
           thresholding for image denoising and compression." Image Processing,
           IEEE Transactions on 9.9 (2000): 1532-1546.
           DOI: 10.1109/83.862633
    .. [2] D. L. Donoho and I. M. Johnstone. "Ideal spatial adaptation
           by wavelet shrinkage." Biometrika 81.3 (1994): 425-455.
           DOI: 10.1093/biomet/81.3.425

    """
    wavelet = pywt.Wavelet(wavelet)

    # original_extent is used to workaround PyWavelets issue #80
    # odd-sized input results in an image with 1 extra sample after waverecn
    original_extent = [slice(s) for s in image.shape]

    # Determine the number of wavelet decomposition levels
    if wavelet_levels is None:
        # Determine the maximum number of possible levels for image
        dlen = wavelet.dec_len
        wavelet_levels = np.min(
            [pywt.dwt_max_level(s, dlen) for s in image.shape])

        # Skip coarsest wavelet scales (see Notes in docstring).
        wavelet_levels = max(wavelet_levels - 3, 1)

    coeffs = pywt.wavedecn(image, wavelet=wavelet, level=wavelet_levels)
    # Detail coefficients at each decomposition level
    dcoeffs = coeffs[1:]

    if sigma is None:
        # Estimate the noise via the method in [2]_
        detail_coeffs = dcoeffs[-1]['d' * image.ndim]
        sigma = _sigma_est_dwt(detail_coeffs, distribution='Gaussian')

    if method is not None and threshold is not None:
        warn(("Thresholding method {} selected.  The user-specified threshold "
              "will be ignored.").format(method))

    if threshold is None:
        var = sigma**2
        if method is None:
            raise ValueError(
                "If method is None, a threshold must be provided.")
        elif method == "BayesShrink":
            # The BayesShrink thresholds from [1]_ in docstring
            threshold = [{key: _bayes_thresh(level[key], var) for key in level}
                         for level in dcoeffs]
        elif method == "VisuShrink":
            # The VisuShrink thresholds from [2]_ in docstring
            threshold = _universal_thresh(image, sigma)
        else:
            raise ValueError("Unrecognized method: {}".format(method))

    if np.isscalar(threshold):
        # A single threshold for all coefficient arrays
        denoised_detail = [{key: pywt.threshold(level[key],
                                                value=threshold,
                                                mode=mode) for key in level}
                           for level in dcoeffs]
    else:
        # Dict of unique threshold coefficients for each detail coeff. array
        denoised_detail = [{key: pywt.threshold(level[key],
                                                value=thresh[key],
                                                mode=mode) for key in level}
                           for thresh, level in zip(threshold, dcoeffs)]
    denoised_coeffs = [coeffs[0]] + denoised_detail
    return pywt.waverecn(denoised_coeffs, wavelet)[original_extent]


def denoise_wavelet(image, sigma=None, wavelet='db1', mode='soft',
                    wavelet_levels=None, multichannel=False,
                    convert2ycbcr=False, method='BayesShrink'):
    """Perform wavelet denoising on an image.

    Parameters
    ----------
    image : ndarray ([M[, N[, ...P]][, C]) of ints, uints or floats
        Input data to be denoised. `image` can be of any numeric type,
        but it is cast into an ndarray of floats for the computation
        of the denoised image.
    sigma : float or list, optional
        The noise standard deviation used when computing the wavelet detail
        coefficient threshold(s). When None (default), the noise standard
        deviation is estimated via the method in [2]_.
    wavelet : string, optional
        The type of wavelet to perform and can be any of the options
        ``pywt.wavelist`` outputs. The default is `'db1'`. For example,
        ``wavelet`` can be any of ``{'db2', 'haar', 'sym9'}`` and many more.
    mode : {'soft', 'hard'}, optional
        An optional argument to choose the type of denoising performed. It
        noted that choosing soft thresholding given additive noise finds the
        best approximation of the original image.
    wavelet_levels : int or None, optional
        The number of wavelet decomposition levels to use.  The default is
        three less than the maximum number of possible decomposition levels.
    multichannel : bool, optional
        Apply wavelet denoising separately for each channel (where channels
        correspond to the final axis of the array).
    convert2ycbcr : bool, optional
        If True and multichannel True, do the wavelet denoising in the YCbCr
        colorspace instead of the RGB color space. This typically results in
        better performance for RGB images.
    method : {'BayesShrink', 'VisuShrink'}, optional
        Thresholding method to be used. The currently supported methods are
        "BayesShrink" [1]_ and "VisuShrink" [2]_. Defaults to "BayesShrink".

    Returns
    -------
    out : ndarray
        Denoised image.

    Notes
    -----
    The wavelet domain is a sparse representation of the image, and can be
    thought of similarly to the frequency domain of the Fourier transform.
    Sparse representations have most values zero or near-zero and truly random
    noise is (usually) represented by many small values in the wavelet domain.
    Setting all values below some threshold to 0 reduces the noise in the
    image, but larger thresholds also decrease the detail present in the image.

    If the input is 3D, this function performs wavelet denoising on each color
    plane separately. The output image is clipped between either [-1, 1] and
    [0, 1] depending on the input image range.

    When YCbCr conversion is done, every color channel is scaled between 0
    and 1, and `sigma` values are applied to these scaled color channels.

    Many wavelet coefficient thresholding approaches have been proposed.  By
    default, ``denoise_wavelet`` applies BayesShrink, which is an adaptive
    thresholding method that computes separate thresholds for each wavelet
    sub-band as described in [1]_.

    If ``method == "VisuShrink"``, a single "universal threshold" is applied to
    all wavelet detail coefficients as described in [2]_.  This threshold
    is designed to remove all Gaussian noise at a given ``sigma`` with high
    probability, but tends to produce images that appear overly smooth.

    References
    ----------
    .. [1] Chang, S. Grace, Bin Yu, and Martin Vetterli. "Adaptive wavelet
           thresholding for image denoising and compression." Image Processing,
           IEEE Transactions on 9.9 (2000): 1532-1546.
           DOI: 10.1109/83.862633
    .. [2] D. L. Donoho and I. M. Johnstone. "Ideal spatial adaptation
           by wavelet shrinkage." Biometrika 81.3 (1994): 425-455.
           DOI: 10.1093/biomet/81.3.425

    Examples
    --------
    >>> from skimage import color, data
    >>> img = img_as_float(data.astronaut())
    >>> img = color.rgb2gray(img)
    >>> img += 0.1 * np.random.randn(*img.shape)
    >>> img = np.clip(img, 0, 1)
    >>> denoised_img = denoise_wavelet(img, sigma=0.1)

    """
    if method not in ["BayesShrink", "VisuShrink"]:
        raise ValueError(
            ('Invalid method: {}. The currently supported methods are '
             '"BayesShrink" and "VisuShrink"').format(method))

    image = img_as_float(image)

    if multichannel:
        if isinstance(sigma, numbers.Number) or sigma is None:
            sigma = [sigma] * image.shape[-1]

    if multichannel:
        if convert2ycbcr:
            out = color.rgb2ycbcr(image)
            for i in range(3):
                # renormalizing this color channel to live in [0, 1]
                min, max = out[..., i].min(), out[..., i].max()
                channel = out[..., i] - min
                channel /= max - min
                out[..., i] = denoise_wavelet(channel, wavelet=wavelet,
                                              method=method, sigma=sigma[i],
                                              mode=mode,
                                              wavelet_levels=wavelet_levels)

                out[..., i] = out[..., i] * (max - min)
                out[..., i] += min
            out = color.ycbcr2rgb(out)
        else:
            out = np.empty_like(image)
            for c in range(image.shape[-1]):
                out[..., c] = _wavelet_threshold(image[..., c],
                                                 wavelet=wavelet,
                                                 method=method,
                                                 sigma=sigma[c], mode=mode,
                                                 wavelet_levels=wavelet_levels)
    else:
        out = _wavelet_threshold(image, wavelet=wavelet, method=method,
                                 sigma=sigma, mode=mode,
                                 wavelet_levels=wavelet_levels)

    clip_range = (-1, 1) if image.min() < 0 else (0, 1)
    return np.clip(out, *clip_range)


def estimate_sigma(image, average_sigmas=False, multichannel=False):
    """
    Robust wavelet-based estimator of the (Gaussian) noise standard deviation.

    Parameters
    ----------
    image : ndarray
        Image for which to estimate the noise standard deviation.
    average_sigmas : bool, optional
        If true, average the channel estimates of `sigma`.  Otherwise return
        a list of sigmas corresponding to each channel.
    multichannel : bool
        Estimate sigma separately for each channel.

    Returns
    -------
    sigma : float or list
        Estimated noise standard deviation(s).  If `multichannel` is True and
        `average_sigmas` is False, a separate noise estimate for each channel
        is returned.  Otherwise, the average of the individual channel
        estimates is returned.

    Notes
    -----
    This function assumes the noise follows a Gaussian distribution. The
    estimation algorithm is based on the median absolute deviation of the
    wavelet detail coefficients as described in section 4.2 of [1]_.

    References
    ----------
    .. [1] D. L. Donoho and I. M. Johnstone. "Ideal spatial adaptation
       by wavelet shrinkage." Biometrika 81.3 (1994): 425-455.
       DOI:10.1093/biomet/81.3.425

    Examples
    --------
    >>> import skimage.data
    >>> from skimage import img_as_float
    >>> img = img_as_float(skimage.data.camera())
    >>> sigma = 0.1
    >>> img = img + sigma * np.random.standard_normal(img.shape)
    >>> sigma_hat = estimate_sigma(img, multichannel=False)
    """
    if multichannel:
        nchannels = image.shape[-1]
        sigmas = [estimate_sigma(
            image[..., c], multichannel=False) for c in range(nchannels)]
        if average_sigmas:
            sigmas = np.mean(sigmas)
        return sigmas
    elif image.shape[-1] <= 4:
        msg = ("image is size {0} on the last axis, but multichannel is "
               "False.  If this is a color image, please set multichannel "
               "to True for proper noise estimation.")
        warn(msg.format(image.shape[-1]))
    coeffs = pywt.dwtn(image, wavelet='db2')
    detail_coeffs = coeffs['d' * image.ndim]
    return _sigma_est_dwt(detail_coeffs, distribution='Gaussian')
