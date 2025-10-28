import numpy as np

from .._shared.utils import _supported_float_type
from ._rolling_ball_cy import apply_kernel, apply_kernel_nan


def rolling_ball(image, *, radius=100, kernel=None, nansafe=False, num_threads=None):
    """Estimate background intensity using the rolling-ball algorithm.

    This function is a generalization of the rolling-ball algorithm [1]_ to
    estimate the background intensity of an n-dimensional image. This is
    typically useful for background subtraction in case of uneven exposure.
    Think of the image as a landscape (where altitude is determined by
    intensity), under which a ball of given radius is rolled. At each
    position, the ball's apex gives the resulting background intensity.

    Parameters
    ----------
    image : ndarray
        The image to be filtered.
    radius : int, optional
        Radius of the ball-shaped kernel to be rolled under the
        image landscape. Used only if `kernel` is ``None``.
    kernel : ndarray, optional
        An alternative way to specify the rolling ball, as an arbitrary
        kernel. It must have the same number of axes as `image`.
    nansafe: bool, optional
        If ``False`` (default), the function assumes that none of the values
        in `image` are ``np.nan``, and uses a faster implementation.
    num_threads: int, optional
        The maximum number of threads to use. If ``None``, the function uses
        the OpenMP default value; typically, it is equal to the maximum number
        of virtual cores.
        Note: This is an upper limit to the number of threads. The exact number
        is determined by the system's OpenMP library.

    Returns
    -------
    background : ndarray
        The estimated background of the image.

    Notes
    -----
    This implementation assumes that dark pixels correspond to the background. If
    you have a bright background, invert the image before passing it to this
    function, e.g., using :func:`skimage.util.invert`.

    For this method to give meaningful results, the radius of the ball (or
    typical size of the kernel, in the general case) should be larger than the
    typical size of the image features of interest.

    This algorithm is sensitive to noise (in particular salt-and-pepper
    noise). If this is a problem in your image, you can apply mild
    Gaussian smoothing before passing the image to this function.

    This algorithm's complexity is polynomial in the radius, with degree equal
    to the image dimensionality (a 2D image is N^2, a 3D image is N^3, etc.),
    so it can take a long time as the radius grows beyond 30 or so ([2]_, [3]_).
    It is an exact N-dimensional calculation; if all you need is an
    approximation, faster options to consider are top-hat filtering [4]_ or
    downscaling-then-upscaling to reduce the size of the input processed.

    References
    ----------
    .. [1] Sternberg, Stanley R. "Biomedical image processing." Computer 1
           (1983): 22-34. :DOI:`10.1109/MC.1983.1654163`
    .. [2] https://github.com/scikit-image/scikit-image/issues/5193
    .. [3] https://github.com/scikit-image/scikit-image/issues/7423
    .. [4] https://forum.image.sc/t/59267/7

    Examples
    --------
    >>> import numpy as np
    >>> import skimage as ski
    >>> image = ski.data.coins()
    >>> background = ski.restoration.rolling_ball(image)
    >>> filtered_image = image - background

    >>> import numpy as np
    >>> import skimage as ski
    >>> image = ski.data.coins()
    >>> kernel = ski.restoration.ellipsoid_kernel((101, 101), 75)
    >>> background = ski.restoration.rolling_ball(image, kernel=kernel)
    >>> filtered_image = image - background
    """

    image = np.asarray(image)
    float_type = _supported_float_type(image.dtype)
    img = image.astype(float_type, copy=False)

    if num_threads is None:
        num_threads = 0

    if kernel is None:
        kernel = ball_kernel(radius, image.ndim)

    kernel = kernel.astype(float_type)
    kernel_shape = np.asarray(kernel.shape)
    kernel_center = kernel_shape // 2
    center_intensity = kernel[tuple(kernel_center)]

    intensity_difference = center_intensity - kernel
    intensity_difference[kernel == np.inf] = np.inf
    intensity_difference = intensity_difference.astype(img.dtype)
    intensity_difference = intensity_difference.reshape(-1)

    img = np.pad(
        img, kernel_center[:, np.newaxis], constant_values=np.inf, mode="constant"
    )

    func = apply_kernel_nan if nansafe else apply_kernel
    background = func(
        img.reshape(-1),
        intensity_difference,
        np.zeros_like(image, dtype=img.dtype).reshape(-1),
        np.array(image.shape, dtype=np.intp),
        np.array(img.shape, dtype=np.intp),
        kernel_shape.astype(np.intp),
        num_threads,
    )

    background = background.astype(image.dtype, copy=False)

    return background


def ball_kernel(radius, ndim):
    """Create a ball kernel for restoration.rolling_ball.

    Parameters
    ----------
    radius : int
        Radius of the ball.
    ndim : int
        Number of dimensions of the ball. ``ndim`` should match the
        dimensionality of the image the kernel will be applied to.

    Returns
    -------
    kernel : ndarray
        The kernel containing the surface intensity of the top half
        of the ellipsoid.

    See Also
    --------
    rolling_ball
    """

    kernel_coords = np.stack(
        np.meshgrid(
            *[np.arange(-x, x + 1) for x in [np.ceil(radius)] * ndim], indexing='ij'
        ),
        axis=-1,
    )

    sum_of_squares = np.sum(kernel_coords**2, axis=-1)
    distance_from_center = np.sqrt(sum_of_squares)
    kernel = np.sqrt(np.clip(radius**2 - sum_of_squares, 0, None))
    kernel[distance_from_center > radius] = np.inf

    return kernel


def ellipsoid_kernel(shape, intensity):
    """Create an ellipoid kernel for restoration.rolling_ball.

    Parameters
    ----------
    shape : array-like
        Length of the principal axis of the ellipsoid (excluding
        the intensity axis). The kernel needs to have the same
        dimensionality as the image it will be applied to.
    intensity : int
        Length of the intensity axis of the ellipsoid.

    Returns
    -------
    kernel : ndarray
        The kernel containing the surface intensity of the top half
        of the ellipsoid.

    See Also
    --------
    rolling_ball
    """

    shape = np.asarray(shape)
    semi_axis = np.clip(shape // 2, 1, None)

    kernel_coords = np.stack(
        np.meshgrid(*[np.arange(-x, x + 1) for x in semi_axis], indexing='ij'), axis=-1
    )

    intensity_scaling = 1 - np.sum((kernel_coords / semi_axis) ** 2, axis=-1)
    kernel = intensity * np.sqrt(np.clip(intensity_scaling, 0, None))
    kernel[intensity_scaling < 0] = np.inf

    return kernel
