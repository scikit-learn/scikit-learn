"""Convex Hull."""
from itertools import product
import numpy as np
from scipy.spatial import ConvexHull
from ..measure.pnpoly import grid_points_in_poly
from ._convex_hull import possible_hull
from ..measure._label import label
from ..util import unique_rows
from .._shared.utils import warn

__all__ = ['convex_hull_image', 'convex_hull_object']


def _offsets_diamond(ndim):
    offsets = np.zeros((2 * ndim, ndim))
    for vertex, (axis, offset) in enumerate(product(range(ndim), (-0.5, 0.5))):
        offsets[vertex, axis] = offset
    return offsets


def convex_hull_image(image, offset_coordinates=True, tolerance=1e-10):
    """Compute the convex hull image of a binary image.

    The convex hull is the set of pixels included in the smallest convex
    polygon that surround all white pixels in the input image.

    Parameters
    ----------
    image : array
        Binary input image. This array is cast to bool before processing.
    offset_coordinates : bool, optional
        If ``True``, a pixel at coordinate, e.g., (4, 7) will be represented
        by coordinates (3.5, 7), (4.5, 7), (4, 6.5), and (4, 7.5). This adds
        some "extent" to a pixel when computing the hull.
    tolerance : float, optional
        Tolerance when determining whether a point is inside the hull. Due
        to numerical floating point errors, a tolerance of 0 can result in
        some points erroneously being classified as being outside the hull.

    Returns
    -------
    hull : (M, N) array of bool
        Binary image with pixels in convex hull set to True.

    References
    ----------
    .. [1] http://blogs.mathworks.com/steve/2011/10/04/binary-image-convex-hull-algorithm-notes/

    """
    ndim = image.ndim
    if np.count_nonzero(image) == 0:
        warn("Input image is entirely zero, no valid convex hull. "
             "Returning empty image", UserWarning)
        return np.zeros(image.shape, dtype=np.bool_)
    # In 2D, we do an optimisation by choosing only pixels that are
    # the starting or ending pixel of a row or column.  This vastly
    # limits the number of coordinates to examine for the virtual hull.
    if ndim == 2:
        coords = possible_hull(image.astype(np.uint8))
    else:
        coords = np.transpose(np.nonzero(image))
        if offset_coordinates:
            # when offsetting, we multiply number of vertices by 2 * ndim.
            # therefore, we reduce the number of coordinates by using a
            # convex hull on the original set, before offsetting.
            hull0 = ConvexHull(coords)
            coords = hull0.points[hull0.vertices]

    # Add a vertex for the middle of each pixel edge
    if offset_coordinates:
        offsets = _offsets_diamond(image.ndim)
        coords = (coords[:, np.newaxis, :] + offsets).reshape(-1, ndim)

    # repeated coordinates can *sometimes* cause problems in
    # scipy.spatial.ConvexHull, so we remove them.
    coords = unique_rows(coords)

    # Find the convex hull
    hull = ConvexHull(coords)
    vertices = hull.points[hull.vertices]

    # If 2D, use fast Cython function to locate convex hull pixels
    if ndim == 2:
        mask = grid_points_in_poly(image.shape, vertices)
    else:
        gridcoords = np.reshape(np.mgrid[tuple(map(slice, image.shape))],
                                (ndim, -1))
        # A point is in the hull if it satisfies all of the hull's inequalities
        coords_in_hull = np.all(hull.equations[:, :ndim].dot(gridcoords) +
                                hull.equations[:, ndim:] < tolerance, axis=0)
        mask = np.reshape(coords_in_hull, image.shape)

    return mask


def convex_hull_object(image, neighbors=8):
    """Compute the convex hull image of individual objects in a binary image.

    The convex hull is the set of pixels included in the smallest convex
    polygon that surround all white pixels in the input image.

    Parameters
    ----------
    image : (M, N) array
        Binary input image.
    neighbors : {4, 8}, int
        Whether to use 4- or 8-connectivity.

    Returns
    -------
    hull : ndarray of bool
        Binary image with pixels in convex hull set to True.

    Notes
    -----
    This function uses skimage.morphology.label to define unique objects,
    finds the convex hull of each using convex_hull_image, and combines
    these regions with logical OR. Be aware the convex hulls of unconnected
    objects may overlap in the result. If this is suspected, consider using
    convex_hull_image separately on each object.

    """
    if image.ndim > 2:
        raise ValueError("Input must be a 2D image")

    if neighbors != 4 and neighbors != 8:
        raise ValueError('Neighbors must be either 4 or 8.')

    labeled_im = label(image, neighbors, background=0)
    convex_obj = np.zeros(image.shape, dtype=bool)
    convex_img = np.zeros(image.shape, dtype=bool)

    for i in range(1, labeled_im.max() + 1):
        convex_obj = convex_hull_image(labeled_im == i)
        convex_img = np.logical_or(convex_img, convex_obj)

    return convex_img
