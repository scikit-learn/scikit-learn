import numpy as np
from ..measure import label


def clear_border(labels, buffer_size=0, bgval=0, in_place=False):
    """Clear objects connected to the label image border.

    Parameters
    ----------
    labels : (M[, N[, ..., P]]) array of int or bool
        Imaging data labels.
    buffer_size : int, optional
        The width of the border examined.  By default, only objects
        that touch the outside of the image are removed.
    bgval : float or int, optional
        Cleared objects are set to this value.
    in_place : bool, optional
        Whether or not to manipulate the labels array in-place.

    Returns
    -------
    out : (M[, N[, ..., P]]) array
        Imaging data labels with cleared borders

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.segmentation import clear_border
    >>> labels = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0],
    ...                    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    ...                    [1, 0, 0, 1, 0, 1, 0, 0, 0],
    ...                    [0, 0, 1, 1, 1, 1, 1, 0, 0],
    ...                    [0, 1, 1, 1, 1, 1, 1, 1, 0],
    ...                    [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> clear_border(labels)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0],
           [0, 1, 1, 1, 1, 1, 1, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0]])

    """
    image = labels

    if any( ( buffer_size >= s for s in image.shape)):
        raise ValueError("buffer size may not be greater than image size")

    # create borders with buffer_size
    borders = np.zeros_like(image, dtype=np.bool_)
    ext = buffer_size + 1
    slstart = slice(ext)
    slend   = slice(-ext, None)
    slices  = [slice(s) for s in image.shape]
    for d in range(image.ndim):
        slicedim = list(slices)
        slicedim[d] = slstart
        borders[slicedim] = True
        slicedim[d] = slend
        borders[slicedim] = True

    # Re-label, in case we are dealing with a binary image
    # and to get consistent labeling
    labels = label(image, background=0)
    number = np.max(labels) + 1

    # determine all objects that are connected to borders
    borders_indices = np.unique(labels[borders])
    indices = np.arange(number + 1)
    # mask all label indices that are connected to borders
    label_mask = np.in1d(indices, borders_indices)
    # create mask for pixels to clear
    mask = label_mask[labels.ravel()].reshape(labels.shape)

    if not in_place:
        image = image.copy()

    # clear border pixels
    image[mask] = bgval

    return image
