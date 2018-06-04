from __future__ import division

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import scipy.ndimage as ndi
from scipy.ndimage.filters import laplace
import skimage
from ..measure import label


def _get_neighborhood(nd_idx, radius, nd_shape):
    bounds_lo = (nd_idx - radius).clip(min=0)
    bounds_hi = (nd_idx + radius + 1).clip(max=nd_shape)
    return bounds_lo, bounds_hi


def _inpaint_biharmonic_single_channel(mask, out, limits):
    # Initialize sparse matrices
    matrix_unknown = sparse.lil_matrix((np.sum(mask), out.size))
    matrix_known = sparse.lil_matrix((np.sum(mask), out.size))

    # Find indexes of masked points in flatten array
    mask_i = np.ravel_multi_index(np.where(mask), mask.shape)

    # Find masked points and prepare them to be easily enumerate over
    mask_pts = np.array(np.where(mask)).T

    # Iterate over masked points
    for mask_pt_n, mask_pt_idx in enumerate(mask_pts):
        # Get bounded neighborhood of selected radius
        b_lo, b_hi = _get_neighborhood(mask_pt_idx, 2, out.shape)

        # Create biharmonic coefficients ndarray
        neigh_coef = np.zeros(b_hi - b_lo)
        neigh_coef[tuple(mask_pt_idx - b_lo)] = 1
        neigh_coef = laplace(laplace(neigh_coef))

        # Iterate over masked point's neighborhood
        it_inner = np.nditer(neigh_coef, flags=['multi_index'])
        for coef in it_inner:
            if coef == 0:
                continue
            tmp_pt_idx = np.add(b_lo, it_inner.multi_index)
            tmp_pt_i = np.ravel_multi_index(tmp_pt_idx, mask.shape)

            if mask[tuple(tmp_pt_idx)]:
                matrix_unknown[mask_pt_n, tmp_pt_i] = coef
            else:
                matrix_known[mask_pt_n, tmp_pt_i] = coef

    # Prepare diagonal matrix
    flat_diag_image = sparse.dia_matrix((out.flatten(), np.array([0])),
                                        shape=(out.size, out.size))

    # Calculate right hand side as a sum of known matrix's columns
    matrix_known = matrix_known.tocsr()
    rhs = -(matrix_known * flat_diag_image).sum(axis=1)

    # Solve linear system for masked points
    matrix_unknown = matrix_unknown[:, mask_i]
    matrix_unknown = sparse.csr_matrix(matrix_unknown)
    result = spsolve(matrix_unknown, rhs)

    # Handle enormous values
    result = np.clip(result, *limits)

    result = result.ravel()

    # Substitute masked points with inpainted versions
    for mask_pt_n, mask_pt_idx in enumerate(mask_pts):
        out[tuple(mask_pt_idx)] = result[mask_pt_n]

    return out


def inpaint_biharmonic(image, mask, multichannel=False):
    """Inpaint masked points in image with biharmonic equations.

    Parameters
    ----------
    image : (M[, N[, ..., P]][, C]) ndarray
        Input image.
    mask : (M[, N[, ..., P]]) ndarray
        Array of pixels to be inpainted. Have to be the same shape as one
        of the 'image' channels. Unknown pixels have to be represented with 1,
        known pixels - with 0.
    multichannel : boolean, optional
        If True, the last `image` dimension is considered as a color channel,
        otherwise as spatial.

    Returns
    -------
    out : (M[, N[, ..., P]][, C]) ndarray
        Input image with masked pixels inpainted.

    References
    ----------
    .. [1]  N.S.Hoang, S.B.Damelin, "On surface completion and image inpainting
            by biharmonic functions: numerical aspects",
            https://arxiv.org/abs/1707.06567
    .. [2]  C. K. Chui and H. N. Mhaskar, MRA Contextual-Recovery Extension of
            Smooth Functions on Manifolds, Appl. and Comp. Harmonic Anal.,
            28 (2010), 104-113,
            DOI: 10.1016/j.acha.2009.04.004

    Examples
    --------
    >>> img = np.tile(np.square(np.linspace(0, 1, 5)), (5, 1))
    >>> mask = np.zeros_like(img)
    >>> mask[2, 2:] = 1
    >>> mask[1, 3:] = 1
    >>> mask[0, 4:] = 1
    >>> out = inpaint_biharmonic(img, mask)
    """

    if image.ndim < 1:
        raise ValueError('Input array has to be at least 1D')

    img_baseshape = image.shape[:-1] if multichannel else image.shape
    if img_baseshape != mask.shape:
        raise ValueError('Input arrays have to be the same shape')

    if np.ma.isMaskedArray(image):
        raise TypeError('Masked arrays are not supported')

    image = skimage.img_as_float(image)
    mask = mask.astype(np.bool)

    # Split inpainting mask into independent regions
    kernel = ndi.morphology.generate_binary_structure(mask.ndim, 1)
    mask_dilated = ndi.morphology.binary_dilation(mask, structure=kernel)
    mask_labeled, num_labels = label(mask_dilated, return_num=True)
    mask_labeled *= mask

    if not multichannel:
        image = image[..., np.newaxis]

    out = np.copy(image)

    for idx_channel in range(image.shape[-1]):
        known_points = image[..., idx_channel][~mask]
        limits = (np.min(known_points), np.max(known_points))

        for idx_region in range(1, num_labels+1):
            mask_region = mask_labeled == idx_region
            _inpaint_biharmonic_single_channel(mask_region,
                out[..., idx_channel], limits)

    if not multichannel:
        out = out[..., 0]

    return out
