import numpy as np

from . import draw


def polygon2mask(image_shape, polygon):
    """Create a binary mask from a polygon.

    Parameters
    ----------
    image_shape : tuple of size 2
        The shape of the mask.
    polygon : (N, 2) array_like
        The polygon coordinates of shape (N, 2) where N is
        the number of points. The coordinates are (row, column).

    Returns
    -------
    mask : 2-D ndarray of type 'bool'
        The binary mask that corresponds to the input polygon.

    See Also
    --------
    polygon:
        Generate coordinates of pixels inside a polygon.

    Notes
    -----
    This function does not do any border checking. Parts of the polygon that
    are outside the coordinate space defined by `image_shape` are not drawn.

    Examples
    --------
    >>> import skimage as ski
    >>> image_shape = (10, 10)
    >>> polygon = np.array([[1, 1], [2, 7], [8, 4]])
    >>> mask = ski.draw.polygon2mask(image_shape, polygon)
    >>> mask.astype(int)
    array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
           [0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    If vertices / points of the `polygon` are outside the coordinate space
    defined by `image_shape`, only a part (or none at all) of the polygon is
    drawn in the mask.

    >>> offset = np.array([[2, -4]])
    >>> ski.draw.polygon2mask(image_shape, polygon - offset).astype(int)
    array([[0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
           [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
           [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
           [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
           [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    """
    polygon = np.asarray(polygon)
    vertex_row_coords, vertex_col_coords = polygon.T
    fill_row_coords, fill_col_coords = draw.polygon(
        vertex_row_coords, vertex_col_coords, image_shape
    )
    mask = np.zeros(image_shape, dtype=bool)
    mask[fill_row_coords, fill_col_coords] = True
    return mask
