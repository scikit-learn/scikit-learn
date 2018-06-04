import numpy as np

__all__ = ['frangi', 'hessian']


def _frangi_hessian_common_filter(image, scale_range, scale_step,
                                  beta1, beta2):
    """This is an intermediate function for Frangi and Hessian filters.

    Shares the common code for Frangi and Hessian functions.

    Parameters
    ----------
    image : (N, M) ndarray
        Array with input image data.
    scale_range : 2-tuple of floats, optional
        The range of sigmas used.
    scale_step : float, optional
        Step size between sigmas.
    beta1 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to deviation from a blob-like structure.
    beta2 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to areas of high variance/texture/structure.

    Returns
    -------
    filtered_list : list
        List of pre-filtered images.

    """
    # Import has to be here due to circular import error
    from ..feature import hessian_matrix, hessian_matrix_eigvals

    sigmas = np.arange(scale_range[0], scale_range[1], scale_step)
    if np.any(np.asarray(sigmas) < 0.0):
        raise ValueError("Sigma values less than zero are not valid")

    beta1 = 2 * beta1 ** 2
    beta2 = 2 * beta2 ** 2

    filtered_array = np.zeros(sigmas.shape + image.shape)
    lambdas_array = np.zeros(sigmas.shape + image.shape)

    # Filtering for all sigmas
    for i, sigma in enumerate(sigmas):
        # Make 2D hessian
        D = hessian_matrix(image, sigma, order='rc')

        # Correct for scale
        D = np.array(D) * (sigma ** 2)

        # Calculate (abs sorted) eigenvalues and vectors
        lambda1, lambda2 = hessian_matrix_eigvals(D)

        # Compute some similarity measures
        lambda1[lambda1 == 0] = 1e-10
        rb = (lambda2 / lambda1) ** 2
        s2 = lambda1 ** 2 + lambda2 ** 2

        # Compute the output image
        filtered = np.exp(-rb / beta1) * (np.ones(np.shape(image)) -
                                          np.exp(-s2 / beta2))

        # Store the results in 3D matrices
        filtered_array[i] = filtered
        lambdas_array[i] = lambda1
    return filtered_array, lambdas_array


def frangi(image, scale_range=(1, 10), scale_step=2, beta1=0.5, beta2=15,
           black_ridges=True):
    """Filter an image with the Frangi filter.

    This filter can be used to detect continuous edges, e.g. vessels,
    wrinkles, rivers. It can be used to calculate the fraction of the
    whole image containing such objects.

    Calculates the eigenvectors of the Hessian to compute the similarity of
    an image region to vessels, according to the method described in [1]_.

    Parameters
    ----------
    image : (N, M) ndarray
        Array with input image data.
    scale_range : 2-tuple of floats, optional
        The range of sigmas used.
    scale_step : float, optional
        Step size between sigmas.
    beta1 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to deviation from a blob-like structure.
    beta2 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to areas of high variance/texture/structure.
    black_ridges : boolean, optional
        When True (the default), the filter detects black ridges; when
        False, it detects white ridges.

    Returns
    -------
    out : (N, M) ndarray
        Filtered image (maximum of pixels across all scales).

    Notes
    -----
    Written by Marc Schrijver, 2/11/2001
    Re-Written by D. J. Kroon University of Twente (May 2009)

    References
    ----------
    .. [1] A. Frangi, W. Niessen, K. Vincken, and M. Viergever. "Multiscale
           vessel enhancement filtering," In LNCS, vol. 1496, pages 130-137,
           Germany, 1998. Springer-Verlag.
    .. [2] Kroon, D.J.: Hessian based Frangi vesselness filter.
    .. [3] http://mplab.ucsd.edu/tutorials/gabor.pdf.
    """
    filtered, lambdas = _frangi_hessian_common_filter(image,
                                                      scale_range, scale_step,
                                                      beta1, beta2)
    if black_ridges:
        filtered[lambdas < 0] = 0
    else:
        filtered[lambdas > 0] = 0

    # Return for every pixel the value of the scale(sigma) with the maximum
    # output pixel value
    return np.max(filtered, axis=0)


def hessian(image, scale_range=(1, 10), scale_step=2, beta1=0.5, beta2=15):
    """Filter an image with the Hessian filter.

    This filter can be used to detect continuous edges, e.g. vessels,
    wrinkles, rivers. It can be used to calculate the fraction of the whole
    image containing such objects.

    Almost equal to Frangi filter, but uses alternative method of smoothing.
    Refer to [1]_ to find the differences between Frangi and Hessian filters.

    Parameters
    ----------
    image : (N, M) ndarray
        Array with input image data.
    scale_range : 2-tuple of floats, optional
        The range of sigmas used.
    scale_step : float, optional
        Step size between sigmas.
    beta1 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to deviation from a blob-like structure.
    beta2 : float, optional
        Frangi correction constant that adjusts the filter's
        sensitivity to areas of high variance/texture/structure.

    Returns
    -------
    out : (N, M) ndarray
        Filtered image (maximum of pixels across all scales).

    Notes
    -----
    Written by Marc Schrijver, 2/11/2001
    Re-Written by D. J. Kroon University of Twente (May 2009)

    References
    ----------
    .. [1] Choon-Ching Ng, Moi Hoon Yap, Nicholas Costen and Baihua Li,
           "Automatic Wrinkle Detection using Hybrid Hessian Filter".
    """

    filtered, lambdas = _frangi_hessian_common_filter(image,
                                                      scale_range, scale_step,
                                                      beta1, beta2)
    filtered[lambdas < 0] = 0

    # Return for every pixel the value of the scale(sigma) with the maximum
    # output pixel value
    out = np.max(filtered, axis=0)
    out[out <= 0] = 1
    return out

