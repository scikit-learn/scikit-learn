import math
import numpy as np
from scipy import ndimage as ndi
from ..transform import resize
from ..util import img_as_float
from ._warps import _multichannel_default


def _smooth(image, sigma, mode, cval, multichannel=None):
    """Return image with each channel smoothed by the Gaussian filter."""
    multichannel = _multichannel_default(multichannel, image.ndim)
    smoothed = np.empty(image.shape, dtype=np.double)

    # apply Gaussian filter to all channels independently
    if multichannel:
        sigma = (sigma, )*(image.ndim - 1) + (0, )
    ndi.gaussian_filter(image, sigma, output=smoothed,
                        mode=mode, cval=cval)
    return smoothed


def _check_factor(factor):
    if factor <= 1:
        raise ValueError('scale factor must be greater than 1')


def pyramid_reduce(image, downscale=2, sigma=None, order=1,
                   mode='reflect', cval=0, multichannel=None):
    """Smooth and then downsample image.

    Parameters
    ----------
    image : ndarray
        Input image.
    downscale : float, optional
        Downscale factor.
    sigma : float, optional
        Sigma for Gaussian filter. Default is `2 * downscale / 6.0` which
        corresponds to a filter mask twice the size of the scale factor that
        covers more than 99% of the Gaussian distribution.
    order : int, optional
        Order of splines used in interpolation of downsampling. See
        `skimage.transform.warp` for detail.
    mode : {'reflect', 'constant', 'edge', 'symmetric', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
    cval : float, optional
        Value to fill past edges of input if mode is 'constant'.
    multichannel : bool, optional
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension. By default, is set to True for
        3D (2D+color) inputs, and False for others. Starting in release 0.16,
        this will always default to False.

    Returns
    -------
    out : array
        Smoothed and downsampled float image.

    References
    ----------
    .. [1] http://web.mit.edu/persci/people/adelson/pub_pdfs/pyramid83.pdf

    """
    multichannel = _multichannel_default(multichannel, image.ndim)
    _check_factor(downscale)

    image = img_as_float(image)

    out_shape = tuple([math.ceil(d / float(downscale)) for d in image.shape])
    if multichannel:
        out_shape = out_shape[:-1]

    if sigma is None:
        # automatically determine sigma which covers > 99% of distribution
        sigma = 2 * downscale / 6.0

    smoothed = _smooth(image, sigma, mode, cval, multichannel)
    out = resize(smoothed, out_shape, order=order, mode=mode, cval=cval,
                 anti_aliasing=False)

    return out


def pyramid_expand(image, upscale=2, sigma=None, order=1,
                   mode='reflect', cval=0, multichannel=None):
    """Upsample and then smooth image.

    Parameters
    ----------
    image : ndarray
        Input image.
    upscale : float, optional
        Upscale factor.
    sigma : float, optional
        Sigma for Gaussian filter. Default is `2 * upscale / 6.0` which
        corresponds to a filter mask twice the size of the scale factor that
        covers more than 99% of the Gaussian distribution.
    order : int, optional
        Order of splines used in interpolation of upsampling. See
        `skimage.transform.warp` for detail.
    mode : {'reflect', 'constant', 'edge', 'symmetric', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
    cval : float, optional
        Value to fill past edges of input if mode is 'constant'.
    multichannel : bool, optional
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension. By default, is set to True for
        3D (2D+color) inputs, and False for others. Starting in release 0.16,
        this will always default to False.


    Returns
    -------
    out : array
        Upsampled and smoothed float image.

    References
    ----------
    .. [1] http://web.mit.edu/persci/people/adelson/pub_pdfs/pyramid83.pdf

    """
    multichannel = _multichannel_default(multichannel, image.ndim)
    _check_factor(upscale)

    image = img_as_float(image)

    out_shape = tuple([math.ceil(upscale * d) for d in image.shape])
    if multichannel:
        out_shape = out_shape[:-1]

    if sigma is None:
        # automatically determine sigma which covers > 99% of distribution
        sigma = 2 * upscale / 6.0

    resized = resize(image, out_shape, order=order,
                     mode=mode, cval=cval, anti_aliasing=False)
    out = _smooth(resized, sigma, mode, cval, multichannel)

    return out


def pyramid_gaussian(image, max_layer=-1, downscale=2, sigma=None, order=1,
                     mode='reflect', cval=0, multichannel=None):
    """Yield images of the Gaussian pyramid formed by the input image.

    Recursively applies the `pyramid_reduce` function to the image, and yields
    the downscaled images.

    Note that the first image of the pyramid will be the original, unscaled
    image. The total number of images is `max_layer + 1`. In case all layers
    are computed, the last image is either a one-pixel image or the image where
    the reduction does not change its shape.

    Parameters
    ----------
    image : ndarray
        Input image.
    max_layer : int
        Number of layers for the pyramid. 0th layer is the original image.
        Default is -1 which builds all possible layers.
    downscale : float, optional
        Downscale factor.
    sigma : float, optional
        Sigma for Gaussian filter. Default is `2 * downscale / 6.0` which
        corresponds to a filter mask twice the size of the scale factor that
        covers more than 99% of the Gaussian distribution.
    order : int, optional
        Order of splines used in interpolation of downsampling. See
        `skimage.transform.warp` for detail.
    mode : {'reflect', 'constant', 'edge', 'symmetric', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
    cval : float, optional
        Value to fill past edges of input if mode is 'constant'.
    multichannel : bool, optional
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension. By default, is set to True for
        3D (2D+color) inputs, and False for others. Starting in release 0.16,
        this will always default to False.


    Returns
    -------
    pyramid : generator
        Generator yielding pyramid layers as float images.

    References
    ----------
    .. [1] http://web.mit.edu/persci/people/adelson/pub_pdfs/pyramid83.pdf

    """
    _check_factor(downscale)

    # cast to float for consistent data type in pyramid
    image = img_as_float(image)

    layer = 0
    current_shape = image.shape

    prev_layer_image = image
    yield image

    # build downsampled images until max_layer is reached or downscale process
    # does not change image size
    while layer != max_layer:
        layer += 1

        layer_image = pyramid_reduce(prev_layer_image, downscale, sigma, order,
                                     mode, cval, multichannel=multichannel)

        prev_shape = np.asarray(current_shape)
        prev_layer_image = layer_image
        current_shape = np.asarray(layer_image.shape)

        # no change to previous pyramid layer
        if np.all(current_shape == prev_shape):
            break

        yield layer_image


def pyramid_laplacian(image, max_layer=-1, downscale=2, sigma=None, order=1,
                      mode='reflect', cval=0, multichannel=None):
    """Yield images of the laplacian pyramid formed by the input image.

    Each layer contains the difference between the downsampled and the
    downsampled, smoothed image::

        layer = resize(prev_layer) - smooth(resize(prev_layer))

    Note that the first image of the pyramid will be the difference between the
    original, unscaled image and its smoothed version. The total number of
    images is `max_layer + 1`. In case all layers are computed, the last image
    is either a one-pixel image or the image where the reduction does not
    change its shape.

    Parameters
    ----------
    image : ndarray
        Input image.
    max_layer : int
        Number of layers for the pyramid. 0th layer is the original image.
        Default is -1 which builds all possible layers.
    downscale : float, optional
        Downscale factor.
    sigma : float, optional
        Sigma for Gaussian filter. Default is `2 * downscale / 6.0` which
        corresponds to a filter mask twice the size of the scale factor that
        covers more than 99% of the Gaussian distribution.
    order : int, optional
        Order of splines used in interpolation of downsampling. See
        `skimage.transform.warp` for detail.
    mode : {'reflect', 'constant', 'edge', 'symmetric', 'wrap'}, optional
        The mode parameter determines how the array borders are handled, where
        cval is the value when mode is equal to 'constant'.
    cval : float, optional
        Value to fill past edges of input if mode is 'constant'.
    multichannel : bool, optional
        Whether the last axis of the image is to be interpreted as multiple
        channels or another spatial dimension. By default, is set to True for
        3D (2D+color) inputs, and False for others. Starting in release 0.16,
        this will always default to False.


    Returns
    -------
    pyramid : generator
        Generator yielding pyramid layers as float images.

    References
    ----------
    .. [1] http://web.mit.edu/persci/people/adelson/pub_pdfs/pyramid83.pdf
    .. [2] http://sepwww.stanford.edu/data/media/public/sep/morgan/texturematch/paper_html/node3.html

    """
    multichannel = _multichannel_default(multichannel, image.ndim)
    _check_factor(downscale)

    # cast to float for consistent data type in pyramid
    image = img_as_float(image)

    if sigma is None:
        # automatically determine sigma which covers > 99% of distribution
        sigma = 2 * downscale / 6.0

    current_shape = image.shape

    smoothed_image = _smooth(image, sigma, mode, cval, multichannel)
    yield image - smoothed_image

    # build downsampled images until max_layer is reached or downscale process
    # does not change image size
    if max_layer == -1:
        max_layer = int(np.ceil(math.log(np.max(current_shape), downscale)))

    for layer in range(max_layer):

        out_shape = tuple(
            [math.ceil(d / float(downscale)) for d in current_shape])

        if multichannel:
            out_shape = out_shape[:-1]

        resized_image = resize(smoothed_image, out_shape, order=order,
                               mode=mode, cval=cval, anti_aliasing=False)
        smoothed_image = _smooth(resized_image, sigma, mode, cval,
                                 multichannel)
        current_shape = np.asarray(resized_image.shape)

        yield resized_image - smoothed_image
