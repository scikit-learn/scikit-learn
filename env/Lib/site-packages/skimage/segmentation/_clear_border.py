import numpy as np

from ..measure import label


def clear_border(labels, buffer_size=0, bgval=0, mask=None, *, out=None):
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
    mask : ndarray of bool, same shape as `image`, optional.
        Image data mask. Objects in labels image overlapping with
        False pixels of mask will be removed. If defined, the
        argument buffer_size will be ignored.
    out : ndarray
        Array of the same shape as `labels`, into which the
        output is placed. By default, a new array is created.

    Returns
    -------
    out : (M[, N[, ..., P]]) array
        Imaging data labels with cleared borders

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.segmentation import clear_border
    >>> labels = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0],
    ...                    [1, 1, 0, 0, 1, 0, 0, 1, 0],
    ...                    [1, 1, 0, 1, 0, 1, 0, 0, 0],
    ...                    [0, 0, 0, 1, 1, 1, 1, 0, 0],
    ...                    [0, 1, 1, 1, 1, 1, 1, 1, 0],
    ...                    [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> clear_border(labels)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 0, 0],
           [0, 1, 1, 1, 1, 1, 1, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> mask = np.array([[0, 0, 1, 1, 1, 1, 1, 1, 1],
    ...                  [0, 0, 1, 1, 1, 1, 1, 1, 1],
    ...                  [1, 1, 1, 1, 1, 1, 1, 1, 1],
    ...                  [1, 1, 1, 1, 1, 1, 1, 1, 1],
    ...                  [1, 1, 1, 1, 1, 1, 1, 1, 1],
    ...                  [1, 1, 1, 1, 1, 1, 1, 1, 1]]).astype(bool)
    >>> clear_border(labels, mask=mask)
    array([[0, 0, 0, 0, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 1, 0, 0, 1, 0],
           [0, 0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 0, 0],
           [0, 1, 1, 1, 1, 1, 1, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0]])

    """
    if any(buffer_size >= s for s in labels.shape) and mask is None:
        # ignore buffer_size if mask
        raise ValueError("buffer size may not be greater than labels size")

    if out is None:
        out = labels.copy()

    if mask is not None:
        err_msg = (
            f'labels and mask should have the same shape but '
            f'are {out.shape} and {mask.shape}'
        )
        if out.shape != mask.shape:
            raise (ValueError, err_msg)
        if mask.dtype != bool:
            raise TypeError("mask should be of type bool.")
        borders = ~mask
    else:
        # create borders with buffer_size
        borders = np.zeros_like(out, dtype=bool)
        ext = buffer_size + 1
        slstart = slice(ext)
        slend = slice(-ext, None)
        slices = [slice(None) for _ in out.shape]
        for d in range(out.ndim):
            slices[d] = slstart
            borders[tuple(slices)] = True
            slices[d] = slend
            borders[tuple(slices)] = True
            slices[d] = slice(None)

    # Re-label, in case we are dealing with a binary out
    # and to get consistent labeling
    labels, number = label(out, background=0, return_num=True)

    # determine all objects that are connected to borders
    borders_indices = np.unique(labels[borders])
    indices = np.arange(number + 1)
    # mask all label indices that are connected to borders
    label_mask = np.isin(indices, borders_indices)
    # create mask for pixels to clear
    mask = label_mask[labels.reshape(-1)].reshape(labels.shape)

    # clear border pixels
    out[mask] = bgval

    return out
