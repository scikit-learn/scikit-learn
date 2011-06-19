"""
Utilities to extract features from images.
"""

# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
from scipy import sparse
from itertools import product
from ..utils.fixes import in1d
from ..utils import check_random_state
from ..base import BaseEstimator

################################################################################
# From an image to a graph

def _make_edges_3d(n_x, n_y, n_z=1):
    """Returns a list of edges for a 3D image.

    Parameters
    ===========
    n_x: integer
        The size of the grid in the x direction.
    n_y: integer
        The size of the grid in the y direction.
    n_z: integer, optional
        The size of the grid in the z direction, defaults to 1
    """
    vertices = np.arange(n_x*n_y*n_z).reshape((n_x, n_y, n_z))
    edges_deep = np.vstack((vertices[:, :, :-1].ravel(),
                            vertices[:, :, 1:].ravel()))
    edges_right = np.vstack((vertices[:, :-1].ravel(), vertices[:, 1:].ravel()))
    edges_down = np.vstack((vertices[:-1].ravel(), vertices[1:].ravel()))
    edges = np.hstack((edges_deep, edges_right, edges_down))
    return edges


def _compute_gradient_3d(edges, img):
    n_x, n_y, n_z = img.shape
    gradient = np.abs(img[edges[0]/(n_y*n_z),
                                (edges[0] % (n_y*n_z))/n_z,
                                (edges[0] % (n_y*n_z))%n_z] -
                           img[edges[1]/(n_y*n_z),
                                (edges[1] % (n_y*n_z))/n_z,
                                (edges[1] % (n_y*n_z)) % n_z])
    return gradient


# XXX: Why mask the image after computing the weights?

def _mask_edges_weights(mask, edges, weights=None):
    """Apply a mask to edges (weighted or not)"""
    inds = np.arange(mask.size)
    inds = inds[mask.ravel()]
    ind_mask = np.logical_and(in1d(edges[0], inds),
                              in1d(edges[1], inds))
    edges = edges[:, ind_mask]
    if weights is not None:
        weights = weights[ind_mask]
    maxval = edges.max()
    order = np.searchsorted(np.unique(edges.ravel()), np.arange(maxval+1))
    edges = order[edges]
    if weights is None:
        return edges
    else:
        return edges, weights


def _to_graph(n_x, n_y, n_z, mask=None, img=None,
              return_as=sparse.coo_matrix, dtype=None):
    """Auxiliary function for img_to_graph and grid_to_graph
    """
    edges = _make_edges_3d(n_x, n_y, n_z)

    if dtype is None:
        if img is None:
            dtype = np.bool
        else:
            dtype = img.dtype

    if img is not None:
        img = np.atleast_3d(img)
        weights = _compute_gradient_3d(edges, img)
        if mask is not None:
            edges, weights = _mask_edges_weights(mask, edges, weights)
            diag = img.squeeze()[mask]
        else:
            diag = img.ravel()
        n_voxels = diag.size
    else:
        if mask is not None:
            edges = _mask_edges_weights(mask, edges)
            n_voxels = np.sum(mask)
        else:
            n_voxels = n_x * n_y * n_z
        weights = np.ones(edges.shape[1], dtype=dtype)
        diag = np.ones(n_voxels, dtype=dtype)

    diag_idx = np.arange(n_voxels)
    i_idx = np.hstack((edges[0], edges[1]))
    j_idx = np.hstack((edges[1], edges[0]))
    graph = sparse.coo_matrix((np.hstack((weights, weights, diag)),
                              (np.hstack((i_idx, diag_idx)),
                               np.hstack((j_idx, diag_idx)))),
                              (n_voxels, n_voxels),
                              dtype=dtype)
    if return_as is np.ndarray:
        return graph.todense()
    return return_as(graph)


def img_to_graph(img, mask=None, return_as=sparse.coo_matrix, dtype=None):
    """Graph of the pixel-to-pixel gradient connections

    Edges are weighted with the gradient values.

    Parameters
    ===========
    img: ndarray, 2D or 3D
        2D or 3D image
    mask : ndarray of booleans, optional
        An optional mask of the image, to consider only part of the
        pixels.
    return_as: np.ndarray or a sparse matrix class, optional
        The class to use to build the returned adjacency matrix.
    dtype: None or dtype, optional
        The data of the returned sparse matrix. By default it is the
        dtype of img
    """
    img = np.atleast_3d(img)
    n_x, n_y, n_z = img.shape
    return _to_graph(n_x, n_y, n_z, mask, img, return_as, dtype)


def grid_to_graph(n_x, n_y, n_z=1, mask=None, return_as=sparse.coo_matrix,
                  dtype=np.bool):
    """Graph of the pixel-to-pixel connections

    Edges exist if 2 voxels are connected.

    Parameters
    ===========
    n_x: int
        Dimension in x axis
    n_y: int
        Dimension in y axis
    n_z: int, optional, default 1
        Dimension in z axis
    mask : ndarray of booleans, optional
        An optional mask of the image, to consider only part of the
        pixels.
    return_as: np.ndarray or a sparse matrix class, optional
        The class to use to build the returned adjacency matrix.
    dtype: dtype, optional, default bool
        The data of the returned sparse matrix. By default it is bool
    """
    return _to_graph(n_x, n_y, n_z, mask=mask, return_as=return_as, dtype=dtype)


################################################################################
# From an image to a set of small image patches

def extract_patches_2d(image, image_size, patch_size, max_patches=None,
                       seed=None):
    """Reshape a 2D image into a collection of patches

    The resulting patches are allocated in a dedicated array.

    Parameters
    ----------
    image: array with shape (i_h, i_w) or (i_h * i_w)
        the original image data

    image_size: tuple of ints (i_h, i_w)
        the dimensions of the images

    patch_size: tuple of ints (p_h, p_w)
        the dimensions of one patch

    max_patches: integer or float, optional default is None
        The maximum number of patches to extract. If max_patches is a float
        between 0 and 1, it is taken to be a proportion of the total number
        of patches.
    
    seed: int or RandomState
        Seed for the random number generator used in case max_patches is used.

    Returns
    -------
    patches: array
         shape is (n_patches, patch_height, patch_width, n_colors)
         or (n_patches, patch_height, patch_width) if n_colors is 1

    Examples
    --------

    >>> one_image = np.arange(16).reshape((4, 4))
    >>> one_image
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])

    >>> patches = extract_patches_2d(one_image, (4, 4), (2, 2))
    >>> patches.shape
    (9, 2, 2)

    >>> patches[0]
    array([[0, 1],
           [4, 5]])

    >>> patches[1]
    array([[1, 2],
           [5, 6]])

    >>> patches[8]
    array([[10, 11],
           [14, 15]])

    """
    i_h, i_w = image_size[:2]
    p_h, p_w = patch_size

    image = np.atleast_2d(image)

    image = image.reshape((i_h, i_w, -1))
    n_colors = image.shape[-1]

    # compute the dimensions of the patches array
    n_h = i_h - p_h + 1
    n_w = i_w - p_w + 1
    all_patches = n_h * n_w
    
    if max_patches:
        if isinstance(max_patches, int) and max_patches < all_patches:
            n_patches = max_patches
        elif isinstance(max_patches, float) and 0 < max_patches < 1:
            n_patches = max_patches * n_patches
        else:
            raise ValueError("Invalid value for max_patches!")
        
        rng = check_random_state(seed)
        patches = np.zeros((n_patches, p_h, p_w, n_colors), dtype=image.dtype)
        i_s = rng.randint(n_h, size=n_patches)
        j_s = rng.randint(n_w, size=n_patches)
        for offset, (i, j) in enumerate(zip(i_s, j_s)):
            patches[offset] = image[i:i + p_h, j:j + p_w, :] 
    else:
        n_patches = all_patches
        patches = np.zeros((n_patches, p_h, p_w, n_colors), dtype=image.dtype)
        offset = 0
        for i in xrange(n_h):
            for j in xrange(n_w):
                patches[offset] = image[i:i + p_h, j:j + p_w, :]
                offset += 1

    # remove the color dimension if useless
    if patches.shape[-1] == 1:
        return patches.reshape((n_patches, p_h, p_w))
    else:
        return patches

def reconstruct_patches(patches, image_size, patch_size):
    """Reconstruct the image from all of its patches"""

    # XXX: make it work with colour images too!

    i_h, i_w = image_size[:2]
    p_h, p_w = patch_size
    img = np.zeros(image_size)
    # compute the dimensions of the patches array
    n_h = i_h - p_h + 1
    n_w = i_w - p_w + 1
    for offset, (i, j) in enumerate(product(xrange(n_h), xrange(n_w))):
        img[i:i + p_h, j:j + p_w] += patches[offset]
    for i in xrange(i_h):
        for j in xrange(i_w):
            img[i, j] /= float(min(i + 1, p_h, i_h - i) * 
                               min(j + 1, p_w, i_w - j))
    return img

class PatchExtractor(BaseEstimator):
    """Extracts patches from a collection of images

    Parameters
    ----------
    patch_size: tuple of ints (p_h, p_w)
        the dimensions of one patch

    max_patches: integer or float, optional default is None
        The maximum number of patches per image to extract. If max_patches is a 
        float in (0, 1), it is taken to mean a proportion of the total number
        of patches.
    
    seed: int or RandomState
        Seed for the random number generator used in case max_patches is used.
    """

    def __init__(self, patch_size, max_patches=None, seed=None):
        self.patch_size = patch_size
        self.max_patches = max_patches
        self.seed = seed
    
    def fit(self, X, y=None):
        return self

    def transform(self, X):
        patches = np.empty((0,) + patch_size)
        for image in X:
            partial_patches = extract_patches_2d(image, image.shape,
                              self.patch_size, self.max_patches, self.seed)
            patches = np.r_[patches, partial_patches]
        return patches
