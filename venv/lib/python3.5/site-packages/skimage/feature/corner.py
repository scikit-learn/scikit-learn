from itertools import combinations_with_replacement

import numpy as np
from scipy import ndimage as ndi
from scipy import stats

from ..util import img_as_float
from ..feature import peak_local_max
from ..feature.util import _prepare_grayscale_input_2D
from ..feature.corner_cy import _corner_fast
from ._hessian_det_appx import _hessian_matrix_det
from ..transform import integral_image
from .._shared.utils import safe_as_int
from .corner_cy import _corner_moravec, _corner_orientations
from warnings import warn


def _compute_derivatives(image, mode='constant', cval=0):
    """Compute derivatives in x and y direction using the Sobel operator.

    Parameters
    ----------
    image : ndarray
        Input image.
    mode : {'constant', 'reflect', 'wrap', 'nearest', 'mirror'}, optional
        How to handle values outside the image borders.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.

    Returns
    -------
    imx : ndarray
        Derivative in x-direction.
    imy : ndarray
        Derivative in y-direction.

    """

    imy = ndi.sobel(image, axis=0, mode=mode, cval=cval)
    imx = ndi.sobel(image, axis=1, mode=mode, cval=cval)

    return imx, imy


def structure_tensor(image, sigma=1, mode='constant', cval=0):
    """Compute structure tensor using sum of squared differences.

    The structure tensor A is defined as::

        A = [Axx Axy]
            [Axy Ayy]

    which is approximated by the weighted sum of squared differences in a local
    window around each pixel in the image.

    Parameters
    ----------
    image : ndarray
        Input image.
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, which is used as a
        weighting function for the local summation of squared differences.
    mode : {'constant', 'reflect', 'wrap', 'nearest', 'mirror'}, optional
        How to handle values outside the image borders.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.

    Returns
    -------
    Axx : ndarray
        Element of the structure tensor for each pixel in the input image.
    Axy : ndarray
        Element of the structure tensor for each pixel in the input image.
    Ayy : ndarray
        Element of the structure tensor for each pixel in the input image.

    Examples
    --------
    >>> from skimage.feature import structure_tensor
    >>> square = np.zeros((5, 5))
    >>> square[2, 2] = 1
    >>> Axx, Axy, Ayy = structure_tensor(square, sigma=0.1)
    >>> Axx
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  1.,  0.],
           [ 0.,  4.,  0.,  4.,  0.],
           [ 0.,  1.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])

    """

    image = _prepare_grayscale_input_2D(image)

    imx, imy = _compute_derivatives(image, mode=mode, cval=cval)

    # structure tensore
    Axx = ndi.gaussian_filter(imx * imx, sigma, mode=mode, cval=cval)
    Axy = ndi.gaussian_filter(imx * imy, sigma, mode=mode, cval=cval)
    Ayy = ndi.gaussian_filter(imy * imy, sigma, mode=mode, cval=cval)

    return Axx, Axy, Ayy


def hessian_matrix(image, sigma=1, mode='constant', cval=0, order=None):
    """Compute Hessian matrix.

    The Hessian matrix is defined as::

        H = [Hrr Hrc]
            [Hrc Hcc]

    which is computed by convolving the image with the second derivatives
    of the Gaussian kernel in the respective x- and y-directions.

    Parameters
    ----------
    image : ndarray
        Input image.
    sigma : float
        Standard deviation used for the Gaussian kernel, which is used as
        weighting function for the auto-correlation matrix.
    mode : {'constant', 'reflect', 'wrap', 'nearest', 'mirror'}, optional
        How to handle values outside the image borders.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.
    order : {'xy', 'rc'}, optional
        This parameter allows for the use of reverse or forward order of
        the image axes in gradient computation. 'xy' indicates the usage
        of the last axis initially (Hxx, Hxy, Hyy), whilst 'rc' indicates
        the use of the first axis initially (Hrr, Hrc, Hcc).

    Returns
    -------
    Hrr : ndarray
        Element of the Hessian matrix for each pixel in the input image.
    Hrc : ndarray
        Element of the Hessian matrix for each pixel in the input image.
    Hcc : ndarray
        Element of the Hessian matrix for each pixel in the input image.

    Examples
    --------
    >>> from skimage.feature import hessian_matrix
    >>> square = np.zeros((5, 5))
    >>> square[2, 2] = 4
    >>> Hrr, Hrc, Hcc = hessian_matrix(square, sigma=0.1, order = 'rc')
    >>> Hrc
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0., -1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0., -1.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    """

    image = img_as_float(image)

    gaussian_filtered = ndi.gaussian_filter(image, sigma=sigma,
                                            mode=mode, cval=cval)

    if order is None:
        if image.ndim == 2:
            # The legacy 2D code followed (x, y) convention, so we swap the axis
            # order to maintain compatibility with old code
            warn('deprecation warning: the default order of the hessian matrix values '
                 'will be "row-column" instead of "xy" starting in skimage version 0.15. '
                 'Use order="rc" or order="xy" to set this explicitly')
            order = 'xy'
        else:
            order = 'rc'


    gradients = np.gradient(gaussian_filtered)
    axes = range(image.ndim)

    if order == 'rc':
        axes = reversed(axes)

    H_elems = [np.gradient(gradients[ax0], axis=ax1)
               for ax0, ax1 in combinations_with_replacement(axes, 2)]

    return H_elems


def _hessian_matrix_image(H_elems):
    """Convert the upper-diagonal elements of the Hessian matrix to a matrix.

    Parameters
    ----------
    H_elems : list of array
        The upper-diagonal elements of the Hessian matrix, as returned by
        `hessian_matrix`.

    Returns
    -------
    hessian_image : array
        An array of shape ``(M, N[, ...], image.ndim, image.ndim)``,
        containing the Hessian matrix corresponding to each coordinate.
    """
    image = H_elems[0]
    hessian_image = np.zeros(image.shape + (image.ndim, image.ndim))
    for idx, (row, col) in \
            enumerate(combinations_with_replacement(range(image.ndim), 2)):
        hessian_image[..., row, col] = H_elems[idx]
        hessian_image[..., col, row] = H_elems[idx]
    return hessian_image


def hessian_matrix_det(image, sigma=1, approximate=True):
    """Compute the approximate Hessian Determinant over an image.

    The 2D approximate method uses box filters over integral images to
    compute the approximate Hessian Determinant, as described in [1]_.

    Parameters
    ----------
    image : array
        The image over which to compute Hessian Determinant.
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, used for the Hessian
        matrix.
    approximate : bool, optional
        If ``True`` and the image is 2D, use a much faster approximate
        computation. This argument has no effect on 3D and higher images.

    Returns
    -------
    out : array
        The array of the Determinant of Hessians.

    References
    ----------
    .. [1] Herbert Bay, Andreas Ess, Tinne Tuytelaars, Luc Van Gool,
           "SURF: Speeded Up Robust Features"
           ftp://ftp.vision.ee.ethz.ch/publications/articles/eth_biwi_00517.pdf

    Notes
    -----
    For 2D images when ``approximate=True``, the running time of this method
    only depends on size of the image. It is independent of `sigma` as one
    would expect. The downside is that the result for `sigma` less than `3`
    is not accurate, i.e., not similar to the result obtained if someone
    computed the Hessian and took its determinant.
    """
    image = img_as_float(image)
    if image.ndim == 2 and approximate:
        integral = integral_image(image)
        return np.array(_hessian_matrix_det(integral, sigma))
    else:  # slower brute-force implementation for nD images
        hessian_mat_array = _hessian_matrix_image(hessian_matrix(image, sigma))
        return np.linalg.det(hessian_mat_array)


def _image_orthogonal_matrix22_eigvals(M00, M01, M11):
    l1 = (M00 + M11) / 2 + np.sqrt(4 * M01 ** 2 + (M00 - M11) ** 2) / 2
    l2 = (M00 + M11) / 2 - np.sqrt(4 * M01 ** 2 + (M00 - M11) ** 2) / 2
    return l1, l2


def structure_tensor_eigvals(Axx, Axy, Ayy):
    """Compute Eigen values of structure tensor.

    Parameters
    ----------
    Axx : ndarray
        Element of the structure tensor for each pixel in the input image.
    Axy : ndarray
        Element of the structure tensor for each pixel in the input image.
    Ayy : ndarray
        Element of the structure tensor for each pixel in the input image.

    Returns
    -------
    l1 : ndarray
        Larger eigen value for each input matrix.
    l2 : ndarray
        Smaller eigen value for each input matrix.

    Examples
    --------
    >>> from skimage.feature import structure_tensor, structure_tensor_eigvals
    >>> square = np.zeros((5, 5))
    >>> square[2, 2] = 1
    >>> Axx, Axy, Ayy = structure_tensor(square, sigma=0.1)
    >>> structure_tensor_eigvals(Axx, Axy, Ayy)[0]
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  2.,  4.,  2.,  0.],
           [ 0.,  4.,  0.,  4.,  0.],
           [ 0.,  2.,  4.,  2.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])

    """

    return _image_orthogonal_matrix22_eigvals(Axx, Axy, Ayy)


def hessian_matrix_eigvals(H_elems, Hxy=None, Hyy=None, Hxx=None):
    """Compute Eigenvalues of Hessian matrix.

    Parameters
    ----------
    H_elems : list of ndarray
        The upper-diagonal elements of the Hessian matrix, as returned
        by `hessian_matrix`.
    Hxy : ndarray, deprecated
        Element of the Hessian matrix for each pixel in the input image.
    Hyy : ndarray, deprecated
        Element of the Hessian matrix for each pixel in the input image.
    Hxx : ndarray, deprecated
        Element of the Hessian matrix for each pixel in the input image.

    Returns
    -------
    eigs : ndarray
        The eigenvalues of the Hessian matrix, in decreasing order. The
        eigenvalues are the leading dimension. That is, ``eigs[i, j, k]``
        contains the ith-largest eigenvalue at position (j, k).

    Examples
    --------
    >>> from skimage.feature import hessian_matrix, hessian_matrix_eigvals
    >>> square = np.zeros((5, 5))
    >>> square[2, 2] = 4
    >>> H_elems = hessian_matrix(square, sigma=0.1, order='rc')
    >>> hessian_matrix_eigvals(H_elems)[0]
    array([[ 0.,  0.,  2.,  0.,  0.],
           [ 0.,  1.,  0.,  1.,  0.],
           [ 2.,  0., -2.,  0.,  2.],
           [ 0.,  1.,  0.,  1.,  0.],
           [ 0.,  0.,  2.,  0.,  0.]])
    """
    if Hxy is not None:
        if Hxx is None:
            Hxx = H_elems
        H_elems = [Hxx, Hxy, Hyy]
        warn('The API of `hessian_matrix_eigvals` has changed. Use a list of '
             'elements instead of separate arguments. The old version of the '
             'API will be removed in version 0.16.')
    if len(H_elems) == 3:  # Use fast Cython code for 2D
        eigvals = np.array(_image_orthogonal_matrix22_eigvals(*H_elems))
    else:
        matrices = _hessian_matrix_image(H_elems)
        # eigvalsh returns eigenvalues in increasing order. We want decreasing
        eigvals = np.linalg.eigvalsh(matrices)[..., ::-1]
        leading_axes = tuple(range(eigvals.ndim - 1))
        eigvals = np.transpose(eigvals, (eigvals.ndim - 1,) + leading_axes)
    return eigvals


def shape_index(image, sigma=1, mode='constant', cval=0):
    """Compute the shape index.

    The shape index, as defined by Koenderink & van Doorn [1]_, is a
    single valued measure of local curvature, assuming the image as a 3D plane
    with intensities representing heights.

    It is derived from the eigen values of the Hessian, and its
    value ranges from -1 to 1 (and is undefined (=NaN) in *flat* regions),
    with following ranges representing following shapes:

    .. table:: Ranges of the shape index and corresponding shapes.

      ===================  =============
      Interval (s in ...)  Shape
      ===================  =============
      [  -1, -7/8)         Spherical cup
      [-7/8, -5/8)         Through
      [-5/8, -3/8)         Rut
      [-3/8, -1/8)         Saddle rut
      [-1/8, +1/8)         Saddle
      [+1/8, +3/8)         Saddle ridge
      [+3/8, +5/8)         Ridge
      [+5/8, +7/8)         Dome
      [+7/8,   +1]         Spherical cap
      ===================  =============

    Parameters
    ----------
    image : ndarray
        Input image.
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, which is used for
        smoothing the input data before Hessian eigen value calculation.
    mode : {'constant', 'reflect', 'wrap', 'nearest', 'mirror'}, optional
        How to handle values outside the image borders
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.

    Returns
    -------
    s : ndarray
        Shape index

    References
    ----------
    .. [1] Koenderink, J. J. & van Doorn, A. J.,
           "Surface shape and curvature scales",
           Image and Vision Computing, 1992, 10, 557-564.
           DOI:10.1016/0262-8856(92)90076-F

    Examples
    --------
    >>> from skimage.feature import shape_index
    >>> square = np.zeros((5, 5))
    >>> square[2, 2] = 4
    >>> s = shape_index(square, sigma=0.1)
    >>> s
    array([[ nan,  nan, -0.5,  nan,  nan],
           [ nan, -0. ,  nan, -0. ,  nan],
           [-0.5,  nan, -1. ,  nan, -0.5],
           [ nan, -0. ,  nan, -0. ,  nan],
           [ nan,  nan, -0.5,  nan,  nan]])
    """

    H = hessian_matrix(image, sigma=sigma, mode=mode, cval=cval, order='rc')
    l1, l2 = hessian_matrix_eigvals(H)

    return (2.0 / np.pi) * np.arctan((l2 + l1) / (l2 - l1))


def corner_kitchen_rosenfeld(image, mode='constant', cval=0):
    """Compute Kitchen and Rosenfeld corner measure response image.

    The corner measure is calculated as follows::

        (imxx * imy**2 + imyy * imx**2 - 2 * imxy * imx * imy)
            / (imx**2 + imy**2)

    Where imx and imy are the first and imxx, imxy, imyy the second
    derivatives.

    Parameters
    ----------
    image : ndarray
        Input image.
    mode : {'constant', 'reflect', 'wrap', 'nearest', 'mirror'}, optional
        How to handle values outside the image borders.
    cval : float, optional
        Used in conjunction with mode 'constant', the value outside
        the image boundaries.

    Returns
    -------
    response : ndarray
        Kitchen and Rosenfeld response image.

    """

    imx, imy = _compute_derivatives(image, mode=mode, cval=cval)
    imxx, imxy = _compute_derivatives(imx, mode=mode, cval=cval)
    imyx, imyy = _compute_derivatives(imy, mode=mode, cval=cval)

    numerator = (imxx * imy ** 2 + imyy * imx ** 2 - 2 * imxy * imx * imy)
    denominator = (imx ** 2 + imy ** 2)

    response = np.zeros_like(image, dtype=np.double)

    mask = denominator != 0
    response[mask] = numerator[mask] / denominator[mask]

    return response


def corner_harris(image, method='k', k=0.05, eps=1e-6, sigma=1):
    """Compute Harris corner measure response image.

    This corner detector uses information from the auto-correlation matrix A::

        A = [(imx**2)   (imx*imy)] = [Axx Axy]
            [(imx*imy)   (imy**2)]   [Axy Ayy]

    Where imx and imy are first derivatives, averaged with a gaussian filter.
    The corner measure is then defined as::

        det(A) - k * trace(A)**2

    or::

        2 * det(A) / (trace(A) + eps)

    Parameters
    ----------
    image : ndarray
        Input image.
    method : {'k', 'eps'}, optional
        Method to compute the response image from the auto-correlation matrix.
    k : float, optional
        Sensitivity factor to separate corners from edges, typically in range
        `[0, 0.2]`. Small values of k result in detection of sharp corners.
    eps : float, optional
        Normalisation factor (Noble's corner measure).
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, which is used as
        weighting function for the auto-correlation matrix.

    Returns
    -------
    response : ndarray
        Harris response image.

    References
    ----------
    .. [1] http://kiwi.cs.dal.ca/~dparks/CornerDetection/harris.htm
    .. [2] http://en.wikipedia.org/wiki/Corner_detection

    Examples
    --------
    >>> from skimage.feature import corner_harris, corner_peaks
    >>> square = np.zeros([10, 10])
    >>> square[2:8, 2:8] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> corner_peaks(corner_harris(square), min_distance=1)
    array([[2, 2],
           [2, 7],
           [7, 2],
           [7, 7]])

    """

    Axx, Axy, Ayy = structure_tensor(image, sigma)

    # determinant
    detA = Axx * Ayy - Axy ** 2
    # trace
    traceA = Axx + Ayy

    if method == 'k':
        response = detA - k * traceA ** 2
    else:
        response = 2 * detA / (traceA + eps)

    return response


def corner_shi_tomasi(image, sigma=1):
    """Compute Shi-Tomasi (Kanade-Tomasi) corner measure response image.

    This corner detector uses information from the auto-correlation matrix A::

        A = [(imx**2)   (imx*imy)] = [Axx Axy]
            [(imx*imy)   (imy**2)]   [Axy Ayy]

    Where imx and imy are first derivatives, averaged with a gaussian filter.
    The corner measure is then defined as the smaller eigenvalue of A::

        ((Axx + Ayy) - sqrt((Axx - Ayy)**2 + 4 * Axy**2)) / 2

    Parameters
    ----------
    image : ndarray
        Input image.
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, which is used as
        weighting function for the auto-correlation matrix.

    Returns
    -------
    response : ndarray
        Shi-Tomasi response image.

    References
    ----------
    .. [1] http://kiwi.cs.dal.ca/~dparks/CornerDetection/harris.htm
    .. [2] http://en.wikipedia.org/wiki/Corner_detection

    Examples
    --------
    >>> from skimage.feature import corner_shi_tomasi, corner_peaks
    >>> square = np.zeros([10, 10])
    >>> square[2:8, 2:8] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> corner_peaks(corner_shi_tomasi(square), min_distance=1)
    array([[2, 2],
           [2, 7],
           [7, 2],
           [7, 7]])

    """

    Axx, Axy, Ayy = structure_tensor(image, sigma)

    # minimum eigenvalue of A
    response = ((Axx + Ayy) - np.sqrt((Axx - Ayy) ** 2 + 4 * Axy ** 2)) / 2

    return response


def corner_foerstner(image, sigma=1):
    """Compute Foerstner corner measure response image.

    This corner detector uses information from the auto-correlation matrix A::

        A = [(imx**2)   (imx*imy)] = [Axx Axy]
            [(imx*imy)   (imy**2)]   [Axy Ayy]

    Where imx and imy are first derivatives, averaged with a gaussian filter.
    The corner measure is then defined as::

        w = det(A) / trace(A)           (size of error ellipse)
        q = 4 * det(A) / trace(A)**2    (roundness of error ellipse)

    Parameters
    ----------
    image : ndarray
        Input image.
    sigma : float, optional
        Standard deviation used for the Gaussian kernel, which is used as
        weighting function for the auto-correlation matrix.

    Returns
    -------
    w : ndarray
        Error ellipse sizes.
    q : ndarray
        Roundness of error ellipse.

    References
    ----------
    .. [1] http://www.ipb.uni-bonn.de/uploads/tx_ikgpublication/foerstner87.fast.pdf
    .. [2] http://en.wikipedia.org/wiki/Corner_detection

    Examples
    --------
    >>> from skimage.feature import corner_foerstner, corner_peaks
    >>> square = np.zeros([10, 10])
    >>> square[2:8, 2:8] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> w, q = corner_foerstner(square)
    >>> accuracy_thresh = 0.5
    >>> roundness_thresh = 0.3
    >>> foerstner = (q > roundness_thresh) * (w > accuracy_thresh) * w
    >>> corner_peaks(foerstner, min_distance=1)
    array([[2, 2],
           [2, 7],
           [7, 2],
           [7, 7]])

    """

    Axx, Axy, Ayy = structure_tensor(image, sigma)

    # determinant
    detA = Axx * Ayy - Axy ** 2
    # trace
    traceA = Axx + Ayy

    w = np.zeros_like(image, dtype=np.double)
    q = np.zeros_like(image, dtype=np.double)

    mask = traceA != 0

    w[mask] = detA[mask] / traceA[mask]
    q[mask] = 4 * detA[mask] / traceA[mask] ** 2

    return w, q


def corner_fast(image, n=12, threshold=0.15):
    """Extract FAST corners for a given image.

    Parameters
    ----------
    image : 2D ndarray
        Input image.
    n : int
        Minimum number of consecutive pixels out of 16 pixels on the circle
        that should all be either brighter or darker w.r.t testpixel.
        A point c on the circle is darker w.r.t test pixel p if
        `Ic < Ip - threshold` and brighter if `Ic > Ip + threshold`. Also
        stands for the n in `FAST-n` corner detector.
    threshold : float
        Threshold used in deciding whether the pixels on the circle are
        brighter, darker or similar w.r.t. the test pixel. Decrease the
        threshold when more corners are desired and vice-versa.

    Returns
    -------
    response : ndarray
        FAST corner response image.

    References
    ----------
    .. [1] Edward Rosten and Tom Drummond
           "Machine Learning for high-speed corner detection",
           http://www.edwardrosten.com/work/rosten_2006_machine.pdf
    .. [2] Wikipedia, "Features from accelerated segment test",
           https://en.wikipedia.org/wiki/Features_from_accelerated_segment_test

    Examples
    --------
    >>> from skimage.feature import corner_fast, corner_peaks
    >>> square = np.zeros((12, 12))
    >>> square[3:9, 3:9] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> corner_peaks(corner_fast(square, 9), min_distance=1)
    array([[3, 3],
           [3, 8],
           [8, 3],
           [8, 8]])

    """
    image = _prepare_grayscale_input_2D(image)

    image = np.ascontiguousarray(image)
    response = _corner_fast(image, n, threshold)
    return response


def corner_subpix(image, corners, window_size=11, alpha=0.99):
    """Determine subpixel position of corners.

    A statistical test decides whether the corner is defined as the
    intersection of two edges or a single peak. Depending on the classification
    result, the subpixel corner location is determined based on the local
    covariance of the grey-values. If the significance level for either
    statistical test is not sufficient, the corner cannot be classified, and
    the output subpixel position is set to NaN.

    Parameters
    ----------
    image : ndarray
        Input image.
    corners : (N, 2) ndarray
        Corner coordinates `(row, col)`.
    window_size : int, optional
        Search window size for subpixel estimation.
    alpha : float, optional
        Significance level for corner classification.

    Returns
    -------
    positions : (N, 2) ndarray
        Subpixel corner positions. NaN for "not classified" corners.

    References
    ----------
    .. [1] http://www.ipb.uni-bonn.de/uploads/tx_ikgpublication/\
           foerstner87.fast.pdf
    .. [2] http://en.wikipedia.org/wiki/Corner_detection

    Examples
    --------
    >>> from skimage.feature import corner_harris, corner_peaks, corner_subpix
    >>> img = np.zeros((10, 10))
    >>> img[:5, :5] = 1
    >>> img[5:, 5:] = 1
    >>> img.astype(int)
    array([[1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]])
    >>> coords = corner_peaks(corner_harris(img), min_distance=2)
    >>> coords_subpix = corner_subpix(img, coords, window_size=7)
    >>> coords_subpix
    array([[ 4.5,  4.5]])

    """

    # window extent in one direction
    wext = (window_size - 1) // 2

    image = np.pad(image, pad_width=wext, mode='constant', constant_values=0)

    # add pad width, make sure to not modify the input values in-place
    corners = safe_as_int(corners + wext)

    # normal equation arrays
    N_dot = np.zeros((2, 2), dtype=np.double)
    N_edge = np.zeros((2, 2), dtype=np.double)
    b_dot = np.zeros((2, ), dtype=np.double)
    b_edge = np.zeros((2, ), dtype=np.double)

    # critical statistical test values
    redundancy = window_size ** 2 - 2
    t_crit_dot = stats.f.isf(1 - alpha, redundancy, redundancy)
    t_crit_edge = stats.f.isf(alpha, redundancy, redundancy)

    # coordinates of pixels within window
    y, x = np.mgrid[- wext:wext + 1, - wext:wext + 1]

    corners_subpix = np.zeros_like(corners, dtype=np.double)

    for i, (y0, x0) in enumerate(corners):

        # crop window around corner + border for sobel operator
        miny = y0 - wext - 1
        maxy = y0 + wext + 2
        minx = x0 - wext - 1
        maxx = x0 + wext + 2
        window = image[miny:maxy, minx:maxx]

        winx, winy = _compute_derivatives(window, mode='constant', cval=0)

        # compute gradient suares and remove border
        winx_winx = (winx * winx)[1:-1, 1:-1]
        winx_winy = (winx * winy)[1:-1, 1:-1]
        winy_winy = (winy * winy)[1:-1, 1:-1]

        # sum of squared differences (mean instead of gaussian filter)
        Axx = np.sum(winx_winx)
        Axy = np.sum(winx_winy)
        Ayy = np.sum(winy_winy)

        # sum of squared differences weighted with coordinates
        # (mean instead of gaussian filter)
        bxx_x = np.sum(winx_winx * x)
        bxx_y = np.sum(winx_winx * y)
        bxy_x = np.sum(winx_winy * x)
        bxy_y = np.sum(winx_winy * y)
        byy_x = np.sum(winy_winy * x)
        byy_y = np.sum(winy_winy * y)

        # normal equations for subpixel position
        N_dot[0, 0] = Axx
        N_dot[0, 1] = N_dot[1, 0] = - Axy
        N_dot[1, 1] = Ayy

        N_edge[0, 0] = Ayy
        N_edge[0, 1] = N_edge[1, 0] = Axy
        N_edge[1, 1] = Axx

        b_dot[:] = bxx_y - bxy_x, byy_x - bxy_y
        b_edge[:] = byy_y + bxy_x, bxx_x + bxy_y

        # estimated positions
        try:
            est_dot = np.linalg.solve(N_dot, b_dot)
            est_edge = np.linalg.solve(N_edge, b_edge)
        except np.linalg.LinAlgError:
            # if image is constant the system is singular
            corners_subpix[i, :] = np.nan, np.nan
            continue

        # residuals
        ry_dot = y - est_dot[0]
        rx_dot = x - est_dot[1]
        ry_edge = y - est_edge[0]
        rx_edge = x - est_edge[1]
        # squared residuals
        rxx_dot = rx_dot * rx_dot
        rxy_dot = rx_dot * ry_dot
        ryy_dot = ry_dot * ry_dot
        rxx_edge = rx_edge * rx_edge
        rxy_edge = rx_edge * ry_edge
        ryy_edge = ry_edge * ry_edge

        # determine corner class (dot or edge)
        # variance for different models
        var_dot = np.sum(winx_winx * ryy_dot - 2 * winx_winy * rxy_dot
                         + winy_winy * rxx_dot)
        var_edge = np.sum(winy_winy * ryy_edge + 2 * winx_winy * rxy_edge
                          + winx_winx * rxx_edge)

        # test value (F-distributed)
        if var_dot < np.spacing(1) and var_edge < np.spacing(1):
            t = np.nan
        elif var_dot == 0:
            t = np.inf
        else:
            t = var_edge / var_dot

        # 1 for edge, -1 for dot, 0 for "not classified"
        corner_class = int(t < t_crit_edge) - int(t > t_crit_dot)

        if corner_class == -1:
            corners_subpix[i, :] = y0 + est_dot[0], x0 + est_dot[1]
        elif corner_class == 0:
            corners_subpix[i, :] = np.nan, np.nan
        elif corner_class == 1:
            corners_subpix[i, :] = y0 + est_edge[0], x0 + est_edge[1]

    # subtract pad width
    corners_subpix -= wext

    return corners_subpix


def corner_peaks(image, min_distance=1, threshold_abs=None, threshold_rel=0.1,
                 exclude_border=True, indices=True, num_peaks=np.inf,
                 footprint=None, labels=None):
    """Find corners in corner measure response image.

    This differs from `skimage.feature.peak_local_max` in that it suppresses
    multiple connected peaks with the same accumulator value.

    Parameters
    ----------
    * : *
        See :py:meth:`skimage.feature.peak_local_max`.

    Examples
    --------
    >>> from skimage.feature import peak_local_max
    >>> response = np.zeros((5, 5))
    >>> response[2:4, 2:4] = 1
    >>> response
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    >>> peak_local_max(response)
    array([[3, 3],
           [3, 2],
           [2, 3],
           [2, 2]])
    >>> corner_peaks(response)
    array([[2, 2]])

    """

    peaks = peak_local_max(image, min_distance=min_distance,
                           threshold_abs=threshold_abs,
                           threshold_rel=threshold_rel,
                           exclude_border=exclude_border,
                           indices=False, num_peaks=num_peaks,
                           footprint=footprint, labels=labels)
    if min_distance > 0:
        coords = np.transpose(peaks.nonzero())
        for r, c in coords:
            if peaks[r, c]:
                peaks[r - min_distance:r + min_distance + 1,
                      c - min_distance:c + min_distance + 1] = False
                peaks[r, c] = True

    if indices is True:
        return np.transpose(peaks.nonzero())
    else:
        return peaks


def corner_moravec(image, window_size=1):
    """Compute Moravec corner measure response image.

    This is one of the simplest corner detectors and is comparatively fast but
    has several limitations (e.g. not rotation invariant).

    Parameters
    ----------
    image : ndarray
        Input image.
    window_size : int, optional
        Window size.

    Returns
    -------
    response : ndarray
        Moravec response image.

    References
    ----------
    .. [1] http://kiwi.cs.dal.ca/~dparks/CornerDetection/moravec.htm
    .. [2] http://en.wikipedia.org/wiki/Corner_detection

    Examples
    --------
    >>> from skimage.feature import corner_moravec
    >>> square = np.zeros([7, 7])
    >>> square[3, 3] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0]])
    >>> corner_moravec(square).astype(int)
    array([[0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 2, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0]])
    """
    return _corner_moravec(image, window_size)


def corner_orientations(image, corners, mask):
    """Compute the orientation of corners.

    The orientation of corners is computed using the first order central moment
    i.e. the center of mass approach. The corner orientation is the angle of
    the vector from the corner coordinate to the intensity centroid in the
    local neighborhood around the corner calculated using first order central
    moment.

    Parameters
    ----------
    image : 2D array
        Input grayscale image.
    corners : (N, 2) array
        Corner coordinates as ``(row, col)``.
    mask : 2D array
        Mask defining the local neighborhood of the corner used for the
        calculation of the central moment.

    Returns
    -------
    orientations : (N, 1) array
        Orientations of corners in the range [-pi, pi].

    References
    ----------
    .. [1] Ethan Rublee, Vincent Rabaud, Kurt Konolige and Gary Bradski
          "ORB : An efficient alternative to SIFT and SURF"
          http://www.vision.cs.chubu.ac.jp/CV-R/pdf/Rublee_iccv2011.pdf
    .. [2] Paul L. Rosin, "Measuring Corner Properties"
          http://users.cs.cf.ac.uk/Paul.Rosin/corner2.pdf

    Examples
    --------
    >>> from skimage.morphology import octagon
    >>> from skimage.feature import (corner_fast, corner_peaks,
    ...                              corner_orientations)
    >>> square = np.zeros((12, 12))
    >>> square[3:9, 3:9] = 1
    >>> square.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> corners = corner_peaks(corner_fast(square, 9), min_distance=1)
    >>> corners
    array([[3, 3],
           [3, 8],
           [8, 3],
           [8, 8]])
    >>> orientations = corner_orientations(square, corners, octagon(3, 2))
    >>> np.rad2deg(orientations)
    array([  45.,  135.,  -45., -135.])

    """
    return _corner_orientations(image, corners, mask)
