"""watershed.py - watershed algorithm

This module implements a watershed algorithm that apportions pixels into
marked basins. The algorithm uses a priority queue to hold the pixels
with the metric for the priority queue being pixel value, then the time
of entry into the queue - this settles ties in favor of the closest marker.

Some ideas taken from
Soille, "Automated Basin Delineation from Digital Elevation Models Using
Mathematical Morphology", Signal Processing 20 (1990) 171-182.

The most important insight in the paper is that entry time onto the queue
solves two problems: a pixel should be assigned to the neighbor with the
largest gradient or, if there is no gradient, pixels on a plateau should
be split between markers on opposite sides.

Originally part of CellProfiler, code licensed under both GPL and BSD licenses.
Website: http://www.cellprofiler.org

Copyright (c) 2003-2009 Massachusetts Institute of Technology
Copyright (c) 2009-2011 Broad Institute
All rights reserved.

Original author: Lee Kamentsky
"""

import numpy as np
from scipy import ndimage as ndi

from . import _watershed
from ..util import crop, regular_seeds


def _validate_inputs(image, markers, mask):
    """Ensure that all inputs to watershed have matching shapes and types.

    Parameters
    ----------
    image : array
        The input image.
    markers : int or array of int
        The marker image.
    mask : array, or None
        A boolean mask, True where we want to compute the watershed.

    Returns
    -------
    image, markers, mask : arrays
        The validated and formatted arrays. Image will have dtype float64,
        markers int32, and mask int8. If ``None`` was given for the mask,
        it is a volume of all 1s.

    Raises
    ------
    ValueError
        If the shapes of the given arrays don't match.
    """
    if not isinstance(markers, (np.ndarray, list, tuple)):
        # not array-like, assume int
        markers = regular_seeds(image.shape, markers)
    elif markers.shape != image.shape:
        raise ValueError("`markers` (shape {}) must have same shape "
                         "as `image` (shape {})".format(markers.shape, image.shape))
    if mask is not None and mask.shape != image.shape:
        raise ValueError("`mask` must have same shape as `image`")
    if mask is None:
        # Use a complete `True` mask if none is provided
        mask = np.ones(image.shape, bool)
    return (image.astype(np.float64),
            markers.astype(np.int32),
            mask.astype(np.int8))


def _validate_connectivity(image_dim, connectivity, offset):
    """Convert any valid connectivity to a structuring element and offset.

    Parameters
    ----------
    image_dim : int
        The number of dimensions of the input image.
    connectivity : int, array, or None
        The neighborhood connectivity. An integer is interpreted as in
        ``scipy.ndimage.generate_binary_structure``, as the maximum number
        of orthogonal steps to reach a neighbor. An array is directly
        interpreted as a structuring element and its shape is validated against
        the input image shape. ``None`` is interpreted as a connectivity of 1.
    offset : tuple of int, or None
        The coordinates of the center of the structuring element.

    Returns
    -------
    c_connectivity : array of bool
        The structuring element corresponding to the input `connectivity`.
    offset : array of int
        The offset corresponding to the center of the structuring element.

    Raises
    ------
    ValueError:
        If the image dimension and the connectivity or offset dimensions don't
        match.
    """
    if connectivity is None:
        connectivity = 1

    if np.isscalar(connectivity):
        c_connectivity = ndi.generate_binary_structure(image_dim, connectivity)
    else:
        c_connectivity = np.array(connectivity, bool)
        if c_connectivity.ndim != image_dim:
            raise ValueError("Connectivity dimension must be same as image")

    if offset is None:
        if any([x % 2 == 0 for x in c_connectivity.shape]):
            raise ValueError("Connectivity array must have an unambiguous "
                             "center")

        offset = np.array(c_connectivity.shape) // 2

    return c_connectivity, offset


def _compute_neighbors(image, structure, offset):
    """Compute neighborhood as an array of linear offsets into the image.

    These are sorted according to Euclidean distance from the center (given
    by `offset`), ensuring that immediate neighbors are visited first.
    """
    structure[tuple(offset)] = 0  # ignore the center; it's not a neighbor
    locations = np.transpose(np.nonzero(structure))
    sqdistances = np.sum((locations - offset)**2, axis=1)
    neighborhood = (np.ravel_multi_index(locations.T, image.shape) -
                    np.ravel_multi_index(offset, image.shape))
    sorted_neighborhood = neighborhood[np.argsort(sqdistances)]
    return sorted_neighborhood


def watershed(image, markers, connectivity=1, offset=None, mask=None,
              compactness=0, watershed_line=False):
    """Find watershed basins in `image` flooded from given `markers`.

    Parameters
    ----------
    image: ndarray (2-D, 3-D, ...) of integers
        Data array where the lowest value points are labeled first.
    markers: int, or ndarray of int, same shape as `image`
        The desired number of markers, or an array marking the basins with the
        values to be assigned in the label matrix. Zero means not a marker.
    connectivity: ndarray, optional
        An array with the same number of dimensions as `image` whose
        non-zero elements indicate neighbors for connection.
        Following the scipy convention, default is a one-connected array of
        the dimension of the image.
    offset: array_like of shape image.ndim, optional
        offset of the connectivity (one offset per dimension)
    mask: ndarray of bools or 0s and 1s, optional
        Array of same shape as `image`. Only points at which mask == True
        will be labeled.
    compactness : float, optional
        Use compact watershed [3]_ with given compactness parameter.
        Higher values result in more regularly-shaped watershed basins.
    watershed_line : bool, optional
        If watershed_line is True, a one-pixel wide line separates the regions
        obtained by the watershed algorithm. The line has the label 0.

    Returns
    -------
    out: ndarray
        A labeled matrix of the same type and shape as markers

    See also
    --------
    skimage.segmentation.random_walker: random walker segmentation
        A segmentation algorithm based on anisotropic diffusion, usually
        slower than the watershed but with good results on noisy data and
        boundaries with holes.

    Notes
    -----
    This function implements a watershed algorithm [1]_ [2]_ that apportions
    pixels into marked basins. The algorithm uses a priority queue to hold
    the pixels with the metric for the priority queue being pixel value, then
    the time of entry into the queue - this settles ties in favor of the
    closest marker.

    Some ideas taken from
    Soille, "Automated Basin Delineation from Digital Elevation Models Using
    Mathematical Morphology", Signal Processing 20 (1990) 171-182

    The most important insight in the paper is that entry time onto the queue
    solves two problems: a pixel should be assigned to the neighbor with the
    largest gradient or, if there is no gradient, pixels on a plateau should
    be split between markers on opposite sides.

    This implementation converts all arguments to specific, lowest common
    denominator types, then passes these to a C algorithm.

    Markers can be determined manually, or automatically using for example
    the local minima of the gradient of the image, or the local maxima of the
    distance function to the background for separating overlapping objects
    (see example).

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Watershed_%28image_processing%29

    .. [2] http://cmm.ensmp.fr/~beucher/wtshed.html

    .. [3] Peer Neubert & Peter Protzel (2014). Compact Watershed and
           Preemptive SLIC: On Improving Trade-offs of Superpixel Segmentation
           Algorithms. ICPR 2014, pp 996-1001. DOI:10.1109/ICPR.2014.181
           https://www.tu-chemnitz.de/etit/proaut/forschung/rsrc/cws_pSLIC_ICPR.pdf

    Examples
    --------
    The watershed algorithm is useful to separate overlapping objects.

    We first generate an initial image with two overlapping circles:

    >>> x, y = np.indices((80, 80))
    >>> x1, y1, x2, y2 = 28, 28, 44, 52
    >>> r1, r2 = 16, 20
    >>> mask_circle1 = (x - x1)**2 + (y - y1)**2 < r1**2
    >>> mask_circle2 = (x - x2)**2 + (y - y2)**2 < r2**2
    >>> image = np.logical_or(mask_circle1, mask_circle2)

    Next, we want to separate the two circles. We generate markers at the
    maxima of the distance to the background:

    >>> from scipy import ndimage as ndi
    >>> distance = ndi.distance_transform_edt(image)
    >>> from skimage.feature import peak_local_max
    >>> local_maxi = peak_local_max(distance, labels=image,
    ...                             footprint=np.ones((3, 3)),
    ...                             indices=False)
    >>> markers = ndi.label(local_maxi)[0]

    Finally, we run the watershed on the image and markers:

    >>> labels = watershed(-distance, markers, mask=image)

    The algorithm works also for 3-D images, and can be used for example to
    separate overlapping spheres.
    """
    image, markers, mask = _validate_inputs(image, markers, mask)
    connectivity, offset = _validate_connectivity(image.ndim, connectivity,
                                                  offset)

    # pad the image, markers, and mask so that we can use the mask to
    # keep from running off the edges
    pad_width = [(p, p) for p in offset]
    image = np.pad(image, pad_width, mode='constant')
    mask = np.pad(mask, pad_width, mode='constant').ravel()
    output = np.pad(markers, pad_width, mode='constant')

    flat_neighborhood = _compute_neighbors(image, connectivity, offset)
    marker_locations = np.flatnonzero(output)
    image_strides = np.array(image.strides, dtype=np.intp) // image.itemsize

    _watershed.watershed_raveled(image.ravel(),
                                 marker_locations, flat_neighborhood,
                                 mask, image_strides, compactness,
                                 output.ravel(),
                                 watershed_line)

    output = crop(output, pad_width, copy=True)

    return output
