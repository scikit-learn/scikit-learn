from ._seam_carving import _seam_carve_v
from .. import util
from .._shared import utils
import numpy as np


def seam_carve(image, energy_map, mode, num, border=1, force_copy=True):
    """ Carve vertical or horizontal seams off an image.

    Carves out vertical/horizontal seams from an image while using the given
    energy map to decide the importance of each pixel.

    Parameters
    ----------
    image : (M, N) or (M, N, 3) ndarray
        Input image whose seams are to be removed.
    energy_map : (M, N) ndarray
        The array to decide the importance of each pixel. The higher
        the value corresponding to a pixel, the more the algorithm will try
        to keep it in the image.
    mode : str {'horizontal', 'vertical'}
        Indicates whether seams are to be removed vertically or horizontally.
        Removing seams horizontally will decrease the height whereas removing
        vertically will decrease the width.
    num : int
        Number of seams are to be removed.
    border : int, optional
        The number of pixels in the right, left and bottom end of the image
        to be excluded from being considered for a seam. This is important as
        certain filters just ignore image boundaries and set them to `0`.
        By default border is set to `1`.
    force_copy : bool, optional
        If set, the `image` and `energy_map` are copied before being used by
        the method which modifies it in place. Set this to `False` if the
        original image and the energy map are no longer needed after
        this operation.

    Returns
    -------
    out : ndarray
        The cropped image with the seams removed.

    References
    ----------
    .. [1] Shai Avidan and Ariel Shamir
           "Seam Carving for Content-Aware Image Resizing"
           http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Avidan07.pdf
    """

    utils.assert_nD(image, (2, 3))
    image = util.img_as_float(image, force_copy)
    energy_map = util.img_as_float(energy_map, force_copy)

    if image.ndim == 2:
        image = image[..., np.newaxis]

    if mode == 'horizontal':
        image = np.swapaxes(image, 0, 1)
        energy_map = np.swapaxes(energy_map, 0, 1)

    image = np.ascontiguousarray(image)
    out = _seam_carve_v(image, energy_map, num, border)

    if mode == 'horizontal':
        out = np.swapaxes(out, 0, 1)

    return np.squeeze(out)
