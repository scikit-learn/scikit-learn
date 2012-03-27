"""Diffusion map for manifold learning"""

# Author: Satrajit Ghosh  -- <satra@mit.edu>
#         based on matlab code by Andy Sweet and Georg Langs (MIT)
# License: BSD, (c) 2012

import numpy as np
import scipy as sp
from scipy.sparse import spdiags

from ..utils.extmath import safe_sparse_dot

from .locally_linear import null_space


def normalized_laplacian_eigs(S, k, alpha=0.5):
    """Compute the p-rank eigendecomposition of a matrix

    The matrix is essentially equivalent to a normalized Laplacian

    Parameters
    ----------
    S: array [n_points, n_points]
        Matrix of similarities between points

    k: int
        Rank of decomposition to find

    Returns
    -------
    V: array [n_points, k]
       Eigenvectors
    L: array [k]
       Eigenvalues
    """
    n_points  = S.shape[0]
    inv_sqrt_D = spdiags(sp.power(sp.sum(S, axis=1), -alpha), 0, n_points, n_points)
    Ms = safe_sparse_dot(inv_sqrt_D, safe_sparse_dot(S, inv_sqrt_D))
    Ms = (Ms + Ms.T)/2
    return null_space(Ms, k)

def generate_diffusion_map(S, n_eigs, dimensionality=None, diffusion_time=None,
                           delta=None):
    """
    S: array [n_points, n_points]
        Matrix of similarities between points

    n_eigs: int (max n_points)
        Number of eigenvectors to use

    dimensionality : int (max n_eigs)
        Dimensionality of the embedding space

    diffusion_time : float
        amount of time to diffuse

    delta : float
        threshold of ratio between extreme eigenvalues

    Either dimensionality, diffusion_time or delta threshold must
    be left unspecified.

    Returns
    -------
    Gamma : [dimensionality, dimensionality]
        data in embedded space

    determined_param : tuple (name, value)
        which undefined parameter was computed

    Psi : narray
       The right eigenvectors of Markov matrix

    Phi : narray
       The left eigenvectors of Markov matrix
    """

    degrees = sp.sum(S, axis=1).todense()
    n_points = len(degrees)
    V, l = normalized_laplacian_eigs(S, n_eigs)
    if dimensionality is None:
        min_dim = np.nonzero(np.power(l/l[1], diffusion_time) <=delta)[0]
        if len(min_dim):
            dimensionality = min_dim[0]
        else:
            dimensionality = len(l)
        determined_param = ('dimensionality', dimensionality)
    elif diffusion_time is None:
        diffusion_time = np.log(delta) / np.log(l[dimensionality-1]/l[1])
        determined_param = ('diffusion_time', diffusion_time)
    elif delta is None:
        delta = np.power(l[dimensionality-1]/l[1], diffusion_time)
        determined_param = ('delta', delta)
    else:
        raise ValueError(('Either dimensionality, diffusion time or delta '
                          'threshold must be left unspecified.'))

    # this ensures that norm wrt stationary distribution (phi_0) is 1
    norm_factor = np.sqrt(np.sum(d))

    # construct left and right eigenvectors of stochastic matrix from
    inv_sqrt_D = spdiags(sp.power(d, -0.5), 0, n_points, n_points)
    sqrt_D = spdiags(sp.power(d, 0.5), 0, n_points, n_points)
    Phi = safe_sparse_dot(sqrt_D, V[:, :dimensionality] / norm_factor)
    Psi = norm_factor * safe_sparse_dot(inv_sqrt_D, V[:, :dimensionality])

    # first right eigenvector should be constant positive, so if it isn't positive, we need to flip all eigenvectors
    if Psi[0, 0] < 0:
        Phi *= -1
        Psi *= -1

    # construct diffusion map coordinates
    Lt = spdiags(sp.power(l[:dimensionality], diffusion_time), 0,
                 dimensionality, dimensionality)
    Gamma = safe_sparse_dot(Psi, Lt)
    return Gamma, determined_param, Phi, Psi
