import numpy as np
from .._shared.utils import warn
from ._nl_means_denoising import (
    _nl_means_denoising_2d,
    _nl_means_denoising_3d,
    _fast_nl_means_denoising_2d,
    _fast_nl_means_denoising_3d)


def denoise_nl_means(image, patch_size=7, patch_distance=11, h=0.1,
                     multichannel=None, fast_mode=True, sigma=0.):
    """
    Perform non-local means denoising on 2-D or 3-D grayscale images, and
    2-D RGB images.

    Parameters
    ----------
    image : 2D or 3D ndarray
        Input image to be denoised, which can be 2D or 3D, and grayscale
        or RGB (for 2D images only, see ``multichannel`` parameter).
    patch_size : int, optional
        Size of patches used for denoising.
    patch_distance : int, optional
        Maximal distance in pixels where to search patches used for denoising.
    h : float, optional
        Cut-off distance (in gray levels). The higher h, the more permissive
        one is in accepting patches. A higher h results in a smoother image,
        at the expense of blurring features. For a Gaussian noise of standard
        deviation sigma, a rule of thumb is to choose the value of h to be
        sigma of slightly less.
    multichannel : bool, optional
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension. Set to ``False`` for 3-D images.
    fast_mode : bool, optional
        If True (default value), a fast version of the non-local means
        algorithm is used. If False, the original version of non-local means is
        used. See the Notes section for more details about the algorithms.
    sigma : float, optional
        The standard deviation of the (Gaussian) noise.  If provided, a more
        robust computation of patch weights is computed that takes the expected
        noise variance into account (see Notes below).

    Returns
    -------

    result : ndarray
        Denoised image, of same shape as `image`.

    Notes
    -----

    The non-local means algorithm is well suited for denoising images with
    specific textures. The principle of the algorithm is to average the value
    of a given pixel with values of other pixels in a limited neighbourhood,
    provided that the *patches* centered on the other pixels are similar enough
    to the patch centered on the pixel of interest.

    In the original version of the algorithm [1]_, corresponding to
    ``fast=False``, the computational complexity is::

        image.size * patch_size ** image.ndim * patch_distance ** image.ndim

    Hence, changing the size of patches or their maximal distance has a
    strong effect on computing times, especially for 3-D images.

    However, the default behavior corresponds to ``fast_mode=True``, for which
    another version of non-local means [2]_ is used, corresponding to a
    complexity of::

        image.size * patch_distance ** image.ndim

    The computing time depends only weakly on the patch size, thanks to
    the computation of the integral of patches distances for a given
    shift, that reduces the number of operations [1]_. Therefore, this
    algorithm executes faster than the classic algorith
    (``fast_mode=False``), at the expense of using twice as much memory.
    This implementation has been proven to be more efficient compared to
    other alternatives, see e.g. [3]_.

    Compared to the classic algorithm, all pixels of a patch contribute
    to the distance to another patch with the same weight, no matter
    their distance to the center of the patch. This coarser computation
    of the distance can result in a slightly poorer denoising
    performance. Moreover, for small images (images with a linear size
    that is only a few times the patch size), the classic algorithm can
    be faster due to boundary effects.

    The image is padded using the `reflect` mode of `skimage.util.pad`
    before denoising.

    If the noise standard deviation, `sigma`, is provided a more robust
    computation of patch weights is used.  Subtracting the known noise variance
    from the computed patch distances improves the estimates of patch
    similarity, giving a moderate improvement to denoising performance [4]_.
    It was also mentioned as an option for the fast variant of the algorithm in
    [3]_.

    When `sigma` is provided, a smaller `h` should typically be used to
    avoid oversmoothing.  The optimal value for `h` depends on the image
    content and noise level, but a reasonable starting point is
    ``h = 0.8 * sigma`` when `fast_mode` is `True`, or ``h = 0.6 * sigma`` when
    `fast_mode` is `False`.

    References
    ----------
    .. [1] A. Buades, B. Coll, & J-M. Morel. A non-local algorithm for image
           denoising. In CVPR 2005, Vol. 2, pp. 60-65, IEEE.
           DOI: 10.1109/CVPR.2005.38

    .. [2] J. Darbon, A. Cunha, T.F. Chan, S. Osher, and G.J. Jensen, Fast
           nonlocal filtering applied to electron cryomicroscopy, in 5th IEEE
           International Symposium on Biomedical Imaging: From Nano to Macro,
           2008, pp. 1331-1334.
           DOI: 10.1109/ISBI.2008.4541250

    .. [3] Jacques Froment. Parameter-Free Fast Pixelwise Non-Local Means
           Denoising. Image Processing On Line, 2014, vol. 4, pp. 300-326.
           DOI: 10.5201/ipol.2014.120

    .. [4] A. Buades, B. Coll, & J-M. Morel. Non-Local Means Denoising.
           Image Processing On Line, 2011, vol. 1, pp. 208-212.
           DOI: 10.5201/ipol.2011.bcm_nlm

    Examples
    --------
    >>> a = np.zeros((40, 40))
    >>> a[10:-10, 10:-10] = 1.
    >>> a += 0.3 * np.random.randn(*a.shape)
    >>> denoised_a = denoise_nl_means(a, 7, 5, 0.1)
    """
    if multichannel is None:
        warn('denoise_nl_means will default to multichannel=False in v0.15')
        multichannel = True

    if image.ndim == 2:
        image = image[..., np.newaxis]
        multichannel = True
    if image.ndim != 3:
        raise NotImplementedError("Non-local means denoising is only \
        implemented for 2D grayscale and RGB images or 3-D grayscale images.")
    nlm_kwargs = dict(s=patch_size, d=patch_distance, h=h, var=sigma * sigma)
    if multichannel:  # 2-D images
        if fast_mode:
            return np.squeeze(
                np.asarray(_fast_nl_means_denoising_2d(image, **nlm_kwargs)))
        else:
            return np.squeeze(
                np.asarray(_nl_means_denoising_2d(image, **nlm_kwargs)))
    else:  # 3-D grayscale
        if fast_mode:
            return np.asarray(_fast_nl_means_denoising_3d(image, **nlm_kwargs))
        else:
            return np.asarray(_nl_means_denoising_3d(image, **nlm_kwargs))
