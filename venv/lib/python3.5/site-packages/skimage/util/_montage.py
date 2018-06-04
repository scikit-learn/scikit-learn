import numpy as np
from .. import exposure
from .._shared.utils import deprecated


__all__ = ['montage', 'montage2d']


def montage(arr_in, fill='mean', rescale_intensity=False, grid_shape=None,
            padding_width=0, multichannel=False):
    """Create a montage of several single- or multichannel images.

    Create a rectangular montage from an input array representing an ensemble
    of equally shaped single- (gray) or multichannel (color) images.

    For example, ``montage(arr_in)`` called with the following `arr_in`

    +---+---+---+
    | 1 | 2 | 3 |
    +---+---+---+

    will return

    +---+---+
    | 1 | 2 |
    +---+---+
    | 3 | * |
    +---+---+

    where the '*' patch will be determined by the `fill` parameter.

    Parameters
    ----------
    arr_in : (K, M, N[, C]) ndarray
        An array representing an ensemble of `K` images of equal shape.
    fill : float or array-like of floats or 'mean', optional
        Value to fill the padding areas and/or the extra tiles in
        the output array. Has to be `float` for single channel collections.
        For multichannel collections has to be an array-like of shape of
        number of channels. If `mean`, uses the mean value over all images.
    rescale_intensity : bool, optional
        Whether to rescale the intensity of each image to [0, 1].
    grid_shape : tuple, optional
        The desired grid shape for the montage `(ntiles_row, ntiles_column)`.
        The default aspect ratio is square.
    padding_width : int, optional
        The size of the spacing between the tiles and between the tiles and
        the borders. If non-zero, makes the boundaries of individual images
        easier to perceive.
    multichannel : boolean, optional
        If True, the last `arr_in` dimension is threated as a color channel,
        otherwise as spatial.

    Returns
    -------
    arr_out : (K*(M+p)+p, K*(N+p)+p[, C]) ndarray
        Output array with input images glued together (including padding `p`).

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.util import montage
    >>> arr_in = np.arange(3 * 2 * 2).reshape(3, 2, 2)
    >>> arr_in  # doctest: +NORMALIZE_WHITESPACE
    array([[[ 0,  1],
            [ 2,  3]],
           [[ 4,  5],
            [ 6,  7]],
           [[ 8,  9],
            [10, 11]]])
    >>> arr_out = montage(arr_in)
    >>> arr_out.shape
    (4, 4)
    >>> arr_out
    array([[ 0,  1,  4,  5],
           [ 2,  3,  6,  7],
           [ 8,  9,  5,  5],
           [10, 11,  5,  5]])
    >>> arr_in.mean()
    5.5
    >>> arr_out_nonsquare = montage(arr_in, grid_shape=(1, 3))
    >>> arr_out_nonsquare
    array([[ 0,  1,  4,  5,  8,  9],
           [ 2,  3,  6,  7, 10, 11]])
    >>> arr_out_nonsquare.shape
    (2, 6)
    """

    if multichannel:
        arr_in = np.asarray(arr_in)
    else:
        arr_in = np.asarray(arr_in)[..., np.newaxis]

    if arr_in.ndim != 4:
        raise ValueError('Input array has to be either 3- or 4-dimensional')

    n_images, n_rows, n_cols, n_chan = arr_in.shape

    if grid_shape:
        ntiles_row, ntiles_col = [int(s) for s in grid_shape]
    else:
        ntiles_row = ntiles_col = int(np.ceil(np.sqrt(n_images)))

    # Rescale intensity if necessary
    if rescale_intensity:
        for i in range(n_images):
            arr_in[i] = exposure.rescale_intensity(arr_in[i])

    # Calculate the fill value
    if fill == 'mean':
        fill = arr_in.mean(axis=(0, 1, 2))
    fill = np.atleast_1d(fill).astype(arr_in.dtype)

    # Pre-allocate an array with padding for montage
    n_pad = padding_width
    arr_out = np.empty(((n_rows + n_pad) * ntiles_row + n_pad,
                        (n_cols + n_pad) * ntiles_col + n_pad,
                        n_chan), dtype=arr_in.dtype)
    for idx_chan in range(n_chan):
        arr_out[..., idx_chan] = fill[idx_chan]

    slices_row = [slice(n_pad + (n_rows + n_pad) * n,
                        n_pad + (n_rows + n_pad) * n + n_rows)
                  for n in range(ntiles_row)]
    slices_col = [slice(n_pad + (n_cols + n_pad) * n,
                        n_pad + (n_cols + n_pad) * n + n_cols)
                  for n in range(ntiles_col)]

    # Copy the data to the output array
    for idx_image, image in enumerate(arr_in):
        idx_sr = idx_image // ntiles_col
        idx_sc = idx_image % ntiles_col
        arr_out[slices_row[idx_sr], slices_col[idx_sc], :] = image

    if multichannel:
        return arr_out
    else:
        return arr_out[..., 0]


@deprecated('montage', removed_version='0.15')
def montage2d(arr_in, fill='mean', rescale_intensity=False, grid_shape=None,
              padding_width=0):
    """Create a 2-dimensional 'montage' from a 3-dimensional input array
    representing an ensemble of equally shaped 2-dimensional images.

    For example, ``montage2d(arr_in, fill)`` with the following `arr_in`

    +---+---+---+
    | 1 | 2 | 3 |
    +---+---+---+

    will return:

    +---+---+
    | 1 | 2 |
    +---+---+
    | 3 | * |
    +---+---+

    Where the '*' patch will be determined by the `fill` parameter.

    Parameters
    ----------
    arr_in : ndarray, shape=[n_images, height, width]
        3-dimensional input array representing an ensemble of n_images
        of equal shape (i.e. [height, width]).
    fill : float or 'mean', optional
        How to fill the 2-dimensional output array when sqrt(n_images)
        is not an integer. If 'mean' is chosen, then fill = arr_in.mean().
    rescale_intensity : bool, optional
        Whether to rescale the intensity of each image to [0, 1].
    grid_shape : tuple, optional
        The desired grid shape for the montage (tiles_y, tiles_x).
        The default aspect ratio is square.
    padding_width : int, optional
        The size of the spacing between the tiles to make the
        boundaries of individual frames easier to see.

    Returns
    -------
    arr_out : ndarray, shape=[alpha * height, alpha * width]
        Output array where 'alpha' has been determined automatically to
        fit (at least) the `n_images` in `arr_in`.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage.util import montage2d
    >>> arr_in = np.arange(3 * 2 * 2).reshape(3, 2, 2)
    >>> arr_in  # doctest: +NORMALIZE_WHITESPACE
    array([[[ 0,  1],
            [ 2,  3]],
           [[ 4,  5],
            [ 6,  7]],
           [[ 8,  9],
            [10, 11]]])
    >>> arr_out = montage2d(arr_in)
    >>> arr_out.shape
    (4, 4)
    >>> arr_out
    array([[ 0,  1,  4,  5],
           [ 2,  3,  6,  7],
           [ 8,  9,  5,  5],
           [10, 11,  5,  5]])
    >>> arr_in.mean()
    5.5
    >>> arr_out_nonsquare = montage2d(arr_in, grid_shape=(1, 3))
    >>> arr_out_nonsquare
    array([[ 0,  1,  4,  5,  8,  9],
           [ 2,  3,  6,  7, 10, 11]])
    >>> arr_out_nonsquare.shape
    (2, 6)
    """

    assert arr_in.ndim == 3

    # -- fill missing patches (needs to be calculated before border padding)
    if fill == 'mean':
        fill = arr_in.mean()

    # -- add border padding, np.pad does all dimensions
    # so we remove the padding from the first
    if padding_width > 0:
        # only pad after to make the width correct
        bef_aft = (0, padding_width)
        arr_in = np.pad(arr_in, ((0, 0), bef_aft, bef_aft), mode='constant')
    else:
        arr_in = arr_in.copy()

    n_images, height, width = arr_in.shape

    # -- rescale intensity if necessary
    if rescale_intensity:
        for i in range(n_images):
            arr_in[i] = exposure.rescale_intensity(arr_in[i])

    # -- determine alpha
    if grid_shape:
        alpha_y, alpha_x = grid_shape
    else:
        alpha_y = alpha_x = int(np.ceil(np.sqrt(n_images)))

    n_missing = int((alpha_y * alpha_x) - n_images)
    # sometimes the mean returns a float, this ensures the missing
    # has the same type for non-float images
    missing = (np.ones((n_missing, height, width), dtype=arr_in.dtype) *
               fill).astype(arr_in.dtype)
    arr_out = np.vstack((arr_in, missing))

    # -- reshape to 2d montage, step by step
    arr_out = arr_out.reshape(alpha_y, alpha_x, height, width)
    arr_out = arr_out.swapaxes(1, 2)
    arr_out = arr_out.reshape(alpha_y * height, alpha_x * width)

    return arr_out
