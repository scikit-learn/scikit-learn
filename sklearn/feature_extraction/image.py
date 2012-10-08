"""
The :mod:`sklearn.feature_extraction.image` submodule gathers utilities to
extract features from images.
"""

# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Olivier Grisel
#          Vlad Niculae
# License: BSD

from itertools import product
import numpy as np
from scipy import sparse

from ..utils.fixes import in1d
from ..utils import array2d, check_random_state
from ..base import BaseEstimator

__all__ = ['PatchExtractor',
           'extract_patches_2d',
           'grid_to_graph',
           'img_to_graph',
           'reconstruct_from_patches_2d']

###############################################################################
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
    vertices = np.arange(n_x * n_y * n_z).reshape((n_x, n_y, n_z))
    edges_deep = np.vstack((vertices[:, :, :-1].ravel(),
                            vertices[:, :, 1:].ravel()))
    edges_right = np.vstack((vertices[:, :-1].ravel(),
                             vertices[:, 1:].ravel()))
    edges_down = np.vstack((vertices[:-1].ravel(), vertices[1:].ravel()))
    edges = np.hstack((edges_deep, edges_right, edges_down))
    return edges


def _compute_gradient_3d(edges, img):
    n_x, n_y, n_z = img.shape
    gradient = np.abs(img[edges[0] // (n_y * n_z),
                                (edges[0] % (n_y * n_z)) // n_z,
                                (edges[0] % (n_y * n_z)) % n_z] -
                                img[edges[1] // (n_y * n_z),
                                (edges[1] % (n_y * n_z)) // n_z,
                                (edges[1] % (n_y * n_z)) % n_z])
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
    if len(edges.ravel()):
        maxval = edges.max()
    else:
        maxval = 0
    order = np.searchsorted(np.unique(edges.ravel()), np.arange(maxval + 1))
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
            dtype = np.int
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
            mask = mask.astype(np.bool)
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
                  dtype=np.int):
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
    dtype: dtype, optional, default int
        The data of the returned sparse matrix. By default it is int
    """
    return _to_graph(n_x, n_y, n_z, mask=mask, return_as=return_as,
                     dtype=dtype)


###############################################################################
# From an image to a set of small image patches


def extract_patches_2d(image, patch_size, max_patches=None, random_state=None):
    """Reshape a 2D image into a collection of patches

    The resulting patches are allocated in a dedicated array.

    Parameters
    ----------
    image: array, shape = (image_height, image_width) or
        (image_height, image_width, n_channels)
        The original image data. For color images, the last dimension specifies
        the channel: a RGB image would have `n_channels=3`.

    patch_size: tuple of ints (patch_height, patch_width)
        the dimensions of one patch

    max_patches: integer or float, optional default is None
        The maximum number of patches to extract. If max_patches is a float
        between 0 and 1, it is taken to be a proportion of the total number
        of patches.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling to use if
        `max_patches` is not None.

    Returns
    -------
    patches: array, shape = (n_patches, patch_height, patch_width) or
         (n_patches, patch_height, patch_width, n_channels)
         The collection of patches extracted from the image, where `n_patches`
         is either `max_patches` or the total number of patches that can be
         extracted.

    Examples
    --------

    >>> from sklearn.feature_extraction import image
    >>> one_image = np.arange(16).reshape((4, 4))
    >>> one_image
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])
    >>> patches = image.extract_patches_2d(one_image, (2, 2))
    >>> print patches.shape
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
    i_h, i_w = image.shape[:2]
    p_h, p_w = patch_size

    image = array2d(image)

    image = image.reshape((i_h, i_w, -1))
    n_colors = image.shape[-1]

    # compute the dimensions of the patches array
    n_h = i_h - p_h + 1
    n_w = i_w - p_w + 1
    all_patches = n_h * n_w

    if max_patches:
        if (isinstance(max_patches, (int, np.integer))
                and max_patches < all_patches):
            n_patches = max_patches
        elif (isinstance(max_patches, (float, np.floating))
                and 0 < max_patches < 1):
            n_patches = int(max_patches * all_patches)
        else:
            raise ValueError("Invalid value for max_patches: %r" % max_patches)

        rng = check_random_state(random_state)
        patches = np.empty((n_patches, p_h, p_w, n_colors), dtype=image.dtype)
        i_s = rng.randint(n_h, size=n_patches)
        j_s = rng.randint(n_w, size=n_patches)
        for p, i, j in zip(patches, i_s, j_s):
            p[:] = image[i:i + p_h, j:j + p_w, :]
    else:
        n_patches = all_patches
        patches = np.empty((n_patches, p_h, p_w, n_colors), dtype=image.dtype)
        for p, (i, j) in zip(patches, product(xrange(n_h), xrange(n_w))):
            p[:] = image[i:i + p_h, j:j + p_w, :]

    # remove the color dimension if useless
    if patches.shape[-1] == 1:
        return patches.reshape((n_patches, p_h, p_w))
    else:
        return patches


def reconstruct_from_patches_2d(patches, image_size):
    """Reconstruct the image from all of its patches.

    Patches are assumed to overlap and the image is constructed by filling in
    the patches from left to right, top to bottom, averaging the overlapping
    regions.

    Parameters
    ----------
    patches: array, shape = (n_patches, patch_height, patch_width) or
        (n_patches, patch_height, patch_width, n_channels)
        The complete set of patches. If the patches contain colour information,
        channels are indexed along the last dimension: RGB patches would
        have `n_channels=3`.

    image_size: tuple of ints (image_height, image_width) or
        (image_height, image_width, n_channels)
        the size of the image that will be reconstructed

    Returns
    -------
    image: array, shape = image_size
        the reconstructed image

    """
    i_h, i_w = image_size[:2]
    p_h, p_w = patches.shape[1:3]
    img = np.zeros(image_size)
    # compute the dimensions of the patches array
    n_h = i_h - p_h + 1
    n_w = i_w - p_w + 1
    for p, (i, j) in zip(patches, product(xrange(n_h), xrange(n_w))):
        img[i:i + p_h, j:j + p_w] += p

    for i in xrange(i_h):
        for j in xrange(i_w):
            # divide by the amount of overlap
            # XXX: is this the most efficient way? memory-wise yes, cpu wise?
            img[i, j] /= float(min(i + 1, p_h, i_h - i) *
                               min(j + 1, p_w, i_w - j))
    return img


class PatchExtractor(BaseEstimator):
    """Extracts patches from a collection of images

    Parameters
    ----------
    patch_size: tuple of ints (patch_height, patch_width)
        the dimensions of one patch

    max_patches: integer or float, optional default is None
        The maximum number of patches per image to extract. If max_patches is a
        float in (0, 1), it is taken to mean a proportion of the total number
        of patches.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    """
    def __init__(self, patch_size=None, max_patches=None, random_state=None):
        self.patch_size = patch_size
        self.max_patches = max_patches
        self.random_state = random_state

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.
        """
        return self

    def transform(self, X):
        """Transforms the image samples in X into a matrix of patch data.

        Parameters
        ----------
        X : array, shape = (n_samples, image_height, image_width) or
            (n_samples, image_height, image_width, n_channels)
            Array of images from which to extract patches. For color images,
            the last dimension specifies the channel: a RGB image would have
            `n_channels=3`.

        Returns
        -------
        patches: array, shape = (n_patches, patch_height, patch_width) or
             (n_patches, patch_height, patch_width, n_channels)
             The collection of patches extracted from the images, where
             `n_patches` is either `n_samples * max_patches` or the total
             number of patches that can be extracted.

        """
        self.random_state = check_random_state(self.random_state)
        n_images, i_h, i_w = X.shape[:3]
        X = np.reshape(X, (n_images, i_h, i_w, -1))
        n_channels = X.shape[-1]
        if self.patch_size is None:
            patch_size = i_h / 10, i_w / 10
        else:
            patch_size = self.patch_size

        if self.max_patches:
            n_patches = self.max_patches
        else:
            p_h, p_w = patch_size
            n_patches = (i_h - p_h + 1) * (i_w - p_w + 1)
        patches_shape = (n_images * n_patches,) + patch_size
        if n_channels > 1:
            patches_shape += (n_channels,)
        patches = np.empty(patches_shape)
        for ii, image in enumerate(X):
            patches[ii * n_patches:(ii + 1) * n_patches] = extract_patches_2d(
                image, patch_size, self.max_patches, self.random_state)
        return patches
