import numpy as np

from ._felzenszwalb_cy import _felzenszwalb_cython


def felzenszwalb(image, scale=1, sigma=0.8, min_size=20, multichannel=True):
    """Computes Felsenszwalb's efficient graph based image segmentation.

    Produces an oversegmentation of a multichannel (i.e. RGB) image
    using a fast, minimum spanning tree based clustering on the image grid.
    The parameter ``scale`` sets an observation level. Higher scale means
    less and larger segments. ``sigma`` is the diameter of a Gaussian kernel,
    used for smoothing the image prior to segmentation.

    The number of produced segments as well as their size can only be
    controlled indirectly through ``scale``. Segment size within an image can
    vary greatly depending on local contrast.

    For RGB images, the algorithm uses the euclidean distance between pixels in
    color space.

    Parameters
    ----------
    image : (width, height, 3) or (width, height) ndarray
        Input image.
    scale : float
        Free parameter. Higher means larger clusters.
    sigma : float
        Width of Gaussian kernel used in preprocessing.
    min_size : int
        Minimum component size. Enforced using postprocessing.
    multichannel : bool, optional (default: True)
        Whether the last axis of the image is to be interpreted as multiple
        channels. A value of False, for a 3D image, is not currently supported.

    Returns
    -------
    segment_mask : (width, height) ndarray
        Integer mask indicating segment labels.

    References
    ----------
    .. [1] Efficient graph-based image segmentation, Felzenszwalb, P.F. and
           Huttenlocher, D.P.  International Journal of Computer Vision, 2004

    Examples
    --------
    >>> from skimage.segmentation import felzenszwalb
    >>> from skimage.data import coffee
    >>> img = coffee()
    >>> segments = felzenszwalb(img, scale=3.0, sigma=0.95, min_size=5)
    """

    if not multichannel and image.ndim > 2:
        raise ValueError("This algorithm works only on single or "
                         "multi-channel 2d images. ")

    image = np.atleast_3d(image)
    return _felzenszwalb_cython(image, scale=scale, sigma=sigma,
                                min_size=min_size)
