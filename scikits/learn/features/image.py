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
    """ Returns a list of edges for a 3D image.
    
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
    gradient = np.abs(img[edges[0]/(n_y*n_z), \
                                (edges[0] % (n_y*n_z))/n_z, \
                                (edges[0] % (n_y*n_z))%n_z] - \
                           img[edges[1]/(n_y*n_z), \
                                (edges[1] % (n_y*n_z))/n_z, \
                                (edges[1] % (n_y*n_z)) % n_z])
    return gradient


# XXX: Notes: what it the difference between the normed and the non
# normed versions of the laplacien?

# XXX: Why mask the image after computing the weights?

def _mask_edges_weights(mask, edges, weights):
    """ Given a image mask and the 
    """
    inds = np.arange(mask.size)
    inds = inds[mask.ravel()]
    ind_mask = np.logical_and(in1d(edges[0], inds),
                              in1d(edges[1], inds))
    edges, weights = edges[:, ind_mask], weights[ind_mask]
    maxval = edges.max()
    order = np.searchsorted(np.unique(edges.ravel()), np.arange(maxval+1))
    edges = order[edges]
    return edges, weights


def img_to_graph(img, mask=None,
                    return_as=sparse.csc_matrix):
    """ Create a graph of the pixel-to-pixel connections with the
        gradient of the image as a the edge value.

        Parameters
        ===========
        img: ndarray, 2D or 3D
            2D or 3D image
        mask : ndarray of booleans, optional
            An optional mask of the image, to consider only part of the 
            pixels.
        return_as: np.ndarray or a sparse matrix class, optional
            The class to use to build the returned adjacency matrix.
    """
    img = np.atleast_3d(img)
    n_x, n_y, n_z = img.shape
    edges   = _make_edges_3d(n_x, n_y, n_z)
    weights = _compute_gradient_3d(edges, img)
    if mask is not None:
        edges, weights = _mask_edges_weights(mask, edges, weights)
        img = img.squeeze()[mask]
    else:
        img = img.ravel()
    n_voxels = img.size
    diag_idx = np.arange(n_voxels)
    i_idx = np.hstack((edges[0], edges[1]))
    j_idx = np.hstack((edges[1], edges[0]))
    graph = sparse.coo_matrix((np.hstack((weights, weights, img)),
                              (np.hstack((i_idx, diag_idx)),
                               np.hstack((j_idx, diag_idx)))),
                              shape=(n_voxels, n_voxels))
    if return_as is np.ndarray:
        return graph.todense()
    return return_as(graph)



################################################################################
# Graph laplacien

#def _make_weights_3d(edges, data, beta=130, eps=1.e-6):
#    weights = np.exp(- beta*gradients / (10*data.std())) + eps
#    return weights


def _make_laplacian_sparse(edges, weights):
    """
    Sparse implementation
    """
    n_voxels = len(np.unique(edges.ravel()))
    diag = np.arange(n_voxels)
    i_indices = np.hstack((edges[0], edges[1]))
    j_indices = np.hstack((edges[1], edges[0]))
    data = np.hstack((-weights, -weights))
    lap = sparse.coo_matrix((data, (i_indices, j_indices)), 
                            shape=(n_voxels, n_voxels))
    connect = - np.ravel(lap.sum(axis=1)) 
    lap = sparse.coo_matrix((np.hstack((data, connect)),
                (np.hstack((i_indices,diag)), np.hstack((j_indices, diag)))), 
                shape=(n_voxels, n_voxels))
    return lap.tocsc()


def _make_normed_laplacian(edges, weights):
    """
    Sparse implementation
    """
    n_voxels = len(np.unique(edges.ravel()))
    i_indices = np.hstack((edges[0], edges[1]))
    j_indices = np.hstack((edges[1], edges[0]))
    data = np.hstack((-weights, -weights))
    lap = sparse.coo_matrix((data, (i_indices, j_indices)), 
                            shape=(n_voxels, n_voxels))
    w = -np.ravel(lap.sum(axis=1))
    data *= 1. / (np.sqrt(w[i_indices]*w[j_indices]))
    lap = sparse.coo_matrix((-data, (i_indices, j_indices)),
                            shape=(n_voxels, n_voxels))
    return lap.tocsc(), w

 
def graph_laplacien(graph, normed=False, beta=50):
    if not normed:
        lap =  _make_laplacian_sparse(edges, weights)
        del edges, weights
        return lap
    else:
        lap, w = _make_normed_laplacian(edges, weights)
        del edges, weights
        return lap, w


