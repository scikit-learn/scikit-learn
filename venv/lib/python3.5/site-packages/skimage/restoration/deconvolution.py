# -*- coding: utf-8 -*-
# deconvolution.py --- Image deconvolution

"""Implementations restoration functions"""

from __future__ import division

import numpy as np
import numpy.random as npr
from scipy.signal import fftconvolve, convolve

from . import uft

__keywords__ = "restoration, image, deconvolution"


def wiener(image, psf, balance, reg=None, is_real=True, clip=True):
    """Wiener-Hunt deconvolution

    Return the deconvolution with a Wiener-Hunt approach (i.e. with
    Fourier diagonalisation).

    Parameters
    ----------
    image : (M, N) ndarray
       Input degraded image
    psf : ndarray
       Point Spread Function. This is assumed to be the impulse
       response (input image space) if the data-type is real, or the
       transfer function (Fourier space) if the data-type is
       complex. There is no constraints on the shape of the impulse
       response. The transfer function must be of shape `(M, N)` if
       `is_real is True`, `(M, N // 2 + 1)` otherwise (see
       `np.fft.rfftn`).
    balance : float
       The regularisation parameter value that tunes the balance
       between the data adequacy that improve frequency restoration
       and the prior adequacy that reduce frequency restoration (to
       avoid noise artifacts).
    reg : ndarray, optional
       The regularisation operator. The Laplacian by default. It can
       be an impulse response or a transfer function, as for the
       psf. Shape constraint is the same as for the `psf` parameter.
    is_real : boolean, optional
       True by default. Specify if ``psf`` and ``reg`` are provided
       with hermitian hypothesis, that is only half of the frequency
       plane is provided (due to the redundancy of Fourier transform
       of real signal). It's apply only if ``psf`` and/or ``reg`` are
       provided as transfer function.  For the hermitian property see
       ``uft`` module or ``np.fft.rfftn``.
    clip : boolean, optional
       True by default. If True, pixel values of the result above 1 or
       under -1 are thresholded for skimage pipeline compatibility.

    Returns
    -------
    im_deconv : (M, N) ndarray
       The deconvolved image.

    Examples
    --------
    >>> from skimage import color, data, restoration
    >>> img = color.rgb2gray(data.astronaut())
    >>> from scipy.signal import convolve2d
    >>> psf = np.ones((5, 5)) / 25
    >>> img = convolve2d(img, psf, 'same')
    >>> img += 0.1 * img.std() * np.random.standard_normal(img.shape)
    >>> deconvolved_img = restoration.wiener(img, psf, 1100)

    Notes
    -----
    This function applies the Wiener filter to a noisy and degraded
    image by an impulse response (or PSF). If the data model is

    .. math:: y = Hx + n

    where :math:`n` is noise, :math:`H` the PSF and :math:`x` the
    unknown original image, the Wiener filter is

    .. math::
       \hat x = F^\dagger (|\Lambda_H|^2 + \lambda |\Lambda_D|^2)
       \Lambda_H^\dagger F y

    where :math:`F` and :math:`F^\dagger` are the Fourier and inverse
    Fourier transfroms respectively, :math:`\Lambda_H` the transfer
    function (or the Fourier transfrom of the PSF, see [Hunt] below)
    and :math:`\Lambda_D` the filter to penalize the restored image
    frequencies (Laplacian by default, that is penalization of high
    frequency). The parameter :math:`\lambda` tunes the balance
    between the data (that tends to increase high frequency, even
    those coming from noise), and the regularization.

    These methods are then specific to a prior model. Consequently,
    the application or the true image nature must corresponds to the
    prior model. By default, the prior model (Laplacian) introduce
    image smoothness or pixel correlation. It can also be interpreted
    as high-frequency penalization to compensate the instability of
    the solution with respect to the data (sometimes called noise
    amplification or "explosive" solution).

    Finally, the use of Fourier space implies a circulant property of
    :math:`H`, see [Hunt].

    References
    ----------
    .. [1] François Orieux, Jean-François Giovannelli, and Thomas
           Rodet, "Bayesian estimation of regularization and point
           spread function parameters for Wiener-Hunt deconvolution",
           J. Opt. Soc. Am. A 27, 1593-1607 (2010)

           http://www.opticsinfobase.org/josaa/abstract.cfm?URI=josaa-27-7-1593

           http://research.orieux.fr/files/papers/OGR-JOSA10.pdf

    .. [2] B. R. Hunt "A matrix theory proof of the discrete
           convolution theorem", IEEE Trans. on Audio and
           Electroacoustics, vol. au-19, no. 4, pp. 285-288, dec. 1971
    """
    if reg is None:
        reg, _ = uft.laplacian(image.ndim, image.shape, is_real=is_real)
    if not np.iscomplexobj(reg):
        reg = uft.ir2tf(reg, image.shape, is_real=is_real)

    if psf.shape != reg.shape:
        trans_func = uft.ir2tf(psf, image.shape, is_real=is_real)
    else:
        trans_func = psf

    wiener_filter = np.conj(trans_func) / (np.abs(trans_func) ** 2 +
                                           balance * np.abs(reg) ** 2)
    if is_real:
        deconv = uft.uirfft2(wiener_filter * uft.urfft2(image),
                             shape=image.shape)
    else:
        deconv = uft.uifft2(wiener_filter * uft.ufft2(image))

    if clip:
        deconv[deconv > 1] = 1
        deconv[deconv < -1] = -1

    return deconv


def unsupervised_wiener(image, psf, reg=None, user_params=None, is_real=True,
                        clip=True):
    """Unsupervised Wiener-Hunt deconvolution.

    Return the deconvolution with a Wiener-Hunt approach, where the
    hyperparameters are automatically estimated. The algorithm is a
    stochastic iterative process (Gibbs sampler) described in the
    reference below. See also ``wiener`` function.

    Parameters
    ----------
    image : (M, N) ndarray
       The input degraded image.
    psf : ndarray
       The impulse response (input image's space) or the transfer
       function (Fourier space). Both are accepted. The transfer
       function is automatically recognized as being complex
       (``np.iscomplexobj(psf)``).
    reg : ndarray, optional
       The regularisation operator. The Laplacian by default. It can
       be an impulse response or a transfer function, as for the psf.
    user_params : dict
       Dictionary of parameters for the Gibbs sampler. See below.
    clip : boolean, optional
       True by default. If true, pixel values of the result above 1 or
       under -1 are thresholded for skimage pipeline compatibility.

    Returns
    -------
    x_postmean : (M, N) ndarray
       The deconvolved image (the posterior mean).
    chains : dict
       The keys ``noise`` and ``prior`` contain the chain list of
       noise and prior precision respectively.

    Other parameters
    ----------------
    The keys of ``user_params`` are:

    threshold : float
       The stopping criterion: the norm of the difference between to
       successive approximated solution (empirical mean of object
       samples, see Notes section). 1e-4 by default.
    burnin : int
       The number of sample to ignore to start computation of the
       mean. 15 by default.
    min_iter : int
       The minimum number of iterations. 30 by default.
    max_iter : int
       The maximum number of iterations if ``threshold`` is not
       satisfied. 200 by default.
    callback : callable (None by default)
       A user provided callable to which is passed, if the function
       exists, the current image sample for whatever purpose. The user
       can store the sample, or compute other moments than the
       mean. It has no influence on the algorithm execution and is
       only for inspection.

    Examples
    --------
    >>> from skimage import color, data, restoration
    >>> img = color.rgb2gray(data.astronaut())
    >>> from scipy.signal import convolve2d
    >>> psf = np.ones((5, 5)) / 25
    >>> img = convolve2d(img, psf, 'same')
    >>> img += 0.1 * img.std() * np.random.standard_normal(img.shape)
    >>> deconvolved_img = restoration.unsupervised_wiener(img, psf)

    Notes
    -----
    The estimated image is design as the posterior mean of a
    probability law (from a Bayesian analysis). The mean is defined as
    a sum over all the possible images weighted by their respective
    probability. Given the size of the problem, the exact sum is not
    tractable. This algorithm use of MCMC to draw image under the
    posterior law. The practical idea is to only draw highly probable
    images since they have the biggest contribution to the mean. At the
    opposite, the less probable images are drawn less often since
    their contribution is low. Finally the empirical mean of these
    samples give us an estimation of the mean, and an exact
    computation with an infinite sample set.

    References
    ----------
    .. [1] François Orieux, Jean-François Giovannelli, and Thomas
           Rodet, "Bayesian estimation of regularization and point
           spread function parameters for Wiener-Hunt deconvolution",
           J. Opt. Soc. Am. A 27, 1593-1607 (2010)

           http://www.opticsinfobase.org/josaa/abstract.cfm?URI=josaa-27-7-1593

           http://research.orieux.fr/files/papers/OGR-JOSA10.pdf
    """
    params = {'threshold': 1e-4, 'max_iter': 200,
              'min_iter': 30, 'burnin': 15, 'callback': None}
    params.update(user_params or {})

    if reg is None:
        reg, _ = uft.laplacian(image.ndim, image.shape, is_real=is_real)
    if not np.iscomplexobj(reg):
        reg = uft.ir2tf(reg, image.shape, is_real=is_real)

    if psf.shape != reg.shape:
        trans_fct = uft.ir2tf(psf, image.shape,  is_real=is_real)
    else:
        trans_fct = psf

    # The mean of the object
    x_postmean = np.zeros(trans_fct.shape)
    # The previous computed mean in the iterative loop
    prev_x_postmean = np.zeros(trans_fct.shape)

    # Difference between two successive mean
    delta = np.NAN

    # Initial state of the chain
    gn_chain, gx_chain = [1], [1]

    # The correlation of the object in Fourier space (if size is big,
    # this can reduce computation time in the loop)
    areg2 = np.abs(reg) ** 2
    atf2 = np.abs(trans_fct) ** 2

    # The Fourier transfrom may change the image.size attribut, so we
    # store it.
    if is_real:
        data_spectrum = uft.urfft2(image.astype(np.float))
    else:
        data_spectrum = uft.ufft2(image.astype(np.float))

    # Gibbs sampling
    for iteration in range(params['max_iter']):
        # Sample of Eq. 27 p(circX^k | gn^k-1, gx^k-1, y).

        # weighting (correlation in direct space)
        precision = gn_chain[-1] * atf2 + gx_chain[-1] * areg2  # Eq. 29
        excursion = np.sqrt(0.5) / np.sqrt(precision) * (
            np.random.standard_normal(data_spectrum.shape) +
            1j * np.random.standard_normal(data_spectrum.shape))

        # mean Eq. 30 (RLS for fixed gn, gamma0 and gamma1 ...)
        wiener_filter = gn_chain[-1] * np.conj(trans_fct) / precision

        # sample of X in Fourier space
        x_sample = wiener_filter * data_spectrum + excursion
        if params['callback']:
            params['callback'](x_sample)

        # sample of Eq. 31 p(gn | x^k, gx^k, y)
        gn_chain.append(npr.gamma(image.size / 2,
                                  2 / uft.image_quad_norm(data_spectrum -
                                                          x_sample *
                                                          trans_fct)))

        # sample of Eq. 31 p(gx | x^k, gn^k-1, y)
        gx_chain.append(npr.gamma((image.size - 1) / 2,
                                  2 / uft.image_quad_norm(x_sample * reg)))

        # current empirical average
        if iteration > params['burnin']:
            x_postmean = prev_x_postmean + x_sample

        if iteration > (params['burnin'] + 1):
            current = x_postmean / (iteration - params['burnin'])
            previous = prev_x_postmean / (iteration - params['burnin'] - 1)

            delta = np.sum(np.abs(current - previous)) / \
                np.sum(np.abs(x_postmean)) / (iteration - params['burnin'])

        prev_x_postmean = x_postmean

        # stop of the algorithm
        if (iteration > params['min_iter']) and (delta < params['threshold']):
            break

    # Empirical average \approx POSTMEAN Eq. 44
    x_postmean = x_postmean / (iteration - params['burnin'])
    if is_real:
        x_postmean = uft.uirfft2(x_postmean, shape=image.shape)
    else:
        x_postmean = uft.uifft2(x_postmean)

    if clip:
        x_postmean[x_postmean > 1] = 1
        x_postmean[x_postmean < -1] = -1

    return (x_postmean, {'noise': gn_chain, 'prior': gx_chain})


def richardson_lucy(image, psf, iterations=50, clip=True):
    """Richardson-Lucy deconvolution.

    Parameters
    ----------
    image : ndarray
       Input degraded image (can be N dimensional).
    psf : ndarray
       The point spread function.
    iterations : int
       Number of iterations. This parameter plays the role of
       regularisation.
    clip : boolean, optional
       True by default. If true, pixel value of the result above 1 or
       under -1 are thresholded for skimage pipeline compatibility.

    Returns
    -------
    im_deconv : ndarray
       The deconvolved image.

    Examples
    --------
    >>> from skimage import color, data, restoration
    >>> camera = color.rgb2gray(data.camera())
    >>> from scipy.signal import convolve2d
    >>> psf = np.ones((5, 5)) / 25
    >>> camera = convolve2d(camera, psf, 'same')
    >>> camera += 0.1 * camera.std() * np.random.standard_normal(camera.shape)
    >>> deconvolved = restoration.richardson_lucy(camera, psf, 5)

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution
    """
    # compute the times for direct convolution and the fft method. The fft is of
    # complexity O(N log(N)) for each dimension and the direct method does
    # straight arithmetic (and is O(n*k) to add n elements k times)
    direct_time = np.prod(image.shape + psf.shape)
    fft_time =  np.sum([n*np.log(n) for n in image.shape + psf.shape])

    # see whether the fourier transform convolution method or the direct
    # convolution method is faster (discussed in scikit-image PR #1792)
    time_ratio = 40.032 * fft_time / direct_time

    if time_ratio <= 1 or len(image.shape) > 2:
        convolve_method = fftconvolve
    else:
        convolve_method = convolve

    image = image.astype(np.float)
    psf = psf.astype(np.float)
    im_deconv = 0.5 * np.ones(image.shape)
    psf_mirror = psf[::-1, ::-1]

    for _ in range(iterations):
        relative_blur = image / convolve_method(im_deconv, psf, 'same')
        im_deconv *= convolve_method(relative_blur, psf_mirror, 'same')

    if clip:
        im_deconv[im_deconv > 1] = 1
        im_deconv[im_deconv < -1] = -1

    return im_deconv
