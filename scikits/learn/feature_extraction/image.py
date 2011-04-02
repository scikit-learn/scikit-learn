"""
Utilities to extract features from images.
"""

# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
from scipy import sparse
from ..utils.fixes import in1d

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


