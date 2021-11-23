"""
Manifold Regularized Discriminative Non-Negative Matrix Factorization.
"""
# Author: Jason Zhang <jason.zhang@ucb.com>
# License: BSD 3 clause

import numbers
import numpy as np
import time
import warnings
from typing import Callable, Dict, Optional
from scipy.linalg import fractional_matrix_power

from sklearn._config import config_context
from sklearn.base import (
    BaseEstimator,
    TransformerMixin,
    _ClassNamePrefixFeaturesOutMixin,
)
from sklearn.decomposition._nmf import _beta_divergence, _check_init
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import check_array, check_random_state
from sklearn.utils.validation import (
    check_is_fitted,
    check_non_negative,
)

EPSILON = np.finfo(np.float32).eps


def mdnmf_loss(X, W, H, Lg, Lc, e, alpha: float, beta: float, gamma: float):
    """
    Calculate the loss function for the given parameters:

        .. math::

            KL(X, WH)

            + 0.5 * alpha * tr(W e W^T)

            + 0.5 * beta * tr(H H^T)

            + 0.5 * gamma * tr(H (L_c^{-0.5}))^T L_g L_c^{-0.5} H^T)

    Note the matrix dimensions are reversed in this function
    (compared with conventional (n_samples, n_features)
    we have (n_features, n_samples) instead), to be aligned
    with most literature in NMF.

    Parameters
    ----------
    X : array-like of shape (n_features, n_samples)
        Location of the samples.

    W : ndarray of shape (n_features, n_components)
        Factor matrix of features.

    H : ndarray of shape (n_components, n_samples)
        Factor matrix of samples.

    Lg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph Laplacian.

    Lc : ndarray of shape (n_samples, n_samples)
        Different-class nearest neighbor graph Laplacian.

    e : ndarray of shape (n_samples, n_samples)
        e = 1 - I, where 1 is a matrix whose elements
        are all one and I is an identity matrix.

    alpha : float
        Parameter for orthogonality regularization.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    Returns
    -------
    loss : float
        Value of the loss function given the parameter values.
    """
    kl_term = _beta_divergence(X, W, H, beta=1)
    ortho_term = 0.5 * alpha * np.trace(W @ e @ W.T)
    tiknohov_term = 0.5 * beta * np.trace(H @ H.T)
    Lc_power_minus_half = np.nan_to_num(fractional_matrix_power(Lc, -0.5))
    md_term = (
        0.5
        * gamma
        * np.trace(H @ Lc_power_minus_half.T @ Lg @ Lc_power_minus_half @ H.T)
    )
    return kl_term + ortho_term + tiknohov_term + md_term


"""
Nearest neighbors graph calculations for
Manifold-regularized Discriminative NMF
"""


def _calc_same_class_dist(X, y):
    """
    Calculate the pairwise distances between the samples
    with the same label, setting ones with different labels
    and distance to oneself equal to max_dist + 1.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Location of the samples.

    y : array-like of shape (n_samples)
        Label of the samples.

    Returns
    -------
    pairwise_dist : ndarray of shape (n_samples, n_samples)
        Pairwise distances between the samples of the same label.
    """
    pairwise_dist = pairwise_distances(X)
    max_dist = np.max(pairwise_dist)
    for i in range(X.shape[0]):
        for j in range(i, X.shape[0]):
            if y[i] != y[j]:
                pairwise_dist[i][j] = pairwise_dist[j][i] = max_dist + 1
            pairwise_dist[i][i] = max_dist + 1
    return pairwise_dist


def _calc_diff_class_dist(X, y):
    """
    Calculate the pairwise distances between the samples
    with different labels, setting ones with the same label
    (including distance to oneself) equal to max_dist + 1.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Location of the samples.

    y : array-like of shape (n_samples)
        Label of the samples.

    Returns
    -------
    pairwise_dist : ndarray of shape (n_samples, n_samples)
        Pairwise distances between the samples of different labels.
    """
    pairwise_dist = pairwise_distances(X)
    max_dist = np.max(pairwise_dist)
    for i in range(X.shape[0]):
        for j in range(i, X.shape[0]):
            if y[i] == y[j]:
                pairwise_dist[i][j] = pairwise_dist[j][i] = max_dist + 1
    return pairwise_dist


def calc_nn_graph_laplacian(X, y, k1: int, k2: int, return_all: bool = False):
    """
    Compute the nearest neighbors graph for same and different labels respectively.
    Then calculate the graph Laplacian matrix for each graph.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Location of the samples.

    y : array-like of shape (n_samples)
        Label of the samples, if not given the results will be all zeros.

    k1 : integer
        Number of nearest neighbors to consider for same-class graph.

    k2 : integer
        Number of nearest neighbors to consider for different-class graph.

    return_all : bool
        Whether to return all the calculated matrices.

    Returns
    -------
    Sg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph matrix.

    Dg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph diagonal matrix.

    Lg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph Laplacian.

    Sc : ndarray of shape (n_samples, n_samples)
        Returned only if return_all is True.
        Different-class nearest neighbor graph matrix.

    Dc : ndarray of shape (n_samples, n_samples)
        Returned only if return_all is True.
        Different-class nearest neighbor graph diagonal matrix.

    Lc : ndarray of shape (n_samples, n_samples)
        Different-class nearest neighbor graph Laplacian.
    """
    n = X.shape[0]
    if y is None:
        return (np.zeros((n, n)),) * 4 if not return_all else (np.zeros((n, n)),) * 6

    same_class_dist = _calc_same_class_dist(X, y)
    diff_class_dist = _calc_diff_class_dist(X, y)

    nn = NearestNeighbors(n_neighbors=k1, metric="precomputed").fit(same_class_dist)
    Sg = nn.kneighbors_graph(same_class_dist).toarray()
    for i in range(Sg.shape[0]):
        for j in range(i, Sg.shape[0]):
            if y[i] != y[j]:
                Sg[i][j] = Sg[j][i] = 0
            Sg[i][i] = 0
    Dg = np.diag(np.sum(Sg, axis=1))
    Lg = Dg - Sg
    Lg = Lg + 1e-4 * np.trace(Lg) * np.eye(Lg.shape[0])  # make sure L is invertible

    nn = NearestNeighbors(n_neighbors=k2, metric="precomputed").fit(diff_class_dist)
    Sc = nn.kneighbors_graph(diff_class_dist).toarray()
    for i in range(Sc.shape[0]):
        for j in range(i, Sc.shape[0]):
            if y[i] == y[j]:
                Sc[i][j] = Sc[j][i] = 0
    Dc = np.diag(np.sum(Sc, axis=1))
    Lc = Dc - Sc
    Lc = Lc + 1e-4 * np.trace(Lc) * np.eye(Lc.shape[0])  # make sure L is invertible

    return (Sg, Dg, Lg, Lc) if not return_all else (Sg, Dg, Lg, Sc, Dc, Lc)


"""
Multiplicative update procedures for
Manifold-regularized Discriminative NMF.

Note the matrix dimensions are reversed in this section
(compared with conventional (n_samples, n_features)
we have (n_features, n_samples) instead), to be aligned
with most literature in NMF.
"""


def _fill_zeros(X):
    """
    Replace zeros in X with a small non-zero value, EPSILON.

    Parameters
    ----------
    X : array-like of any shape
        Target matrix to be filled.

    Returns
    -------
    X : array-like of any shape
        Target matrix with zeros filled with EPSILON.
    """
    X[X == 0] = EPSILON
    return X


def _grad_w(X, W, H, E, e, alpha: float):
    """
    Calculate the factor matrix to be multiplied by W (the gradient):

        .. math::

            \\frac{\\frac{X}{W H} H^T}{E H^T + \\alpha W e}

    Parameters
    ----------
    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    W : ndarray of shape (n_features, n_components)
        Factor matrix to be updated.

    H : ndarray of shape (n_components, n_samples)
        Factor matrix used for update calculation.

    E : ndarray of shape (n_features, n_samples)
        A rectangle matrix whose elements are all one.

    e : ndarray of shape (n_samples, n_samples)
        e = 1 - I, where 1 is a matrix whose elements
        are all one and I is an identity matrix.

    alpha : float
        Parameter for orthogonality regularization.

    Returns
    -------
    delta_W : ndarray of shape (n_features, n_components)
        Factor matrix to be multiplied by W.
    """
    nominator = (X / (_fill_zeros(W @ H))) @ H.T
    denominator = _fill_zeros(E @ H.T + alpha * (W @ e))
    return np.divide(nominator, denominator)


def _grad_h(X, W, H, S, D, E, beta: float, gamma: float):
    """
    Calculate the factor matrix to be multiplied by H (the gradient):

        .. math::

            \\frac{W^T \\frac{X}{W H} + \\gamma H S}{W^T E + \\beta H + \\gamma H D}

    Parameters
    ----------
    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    W : ndarray of shape (n_features, n_components)
        Factor matrix to be updated.

    H : ndarray of shape (n_components, n_samples)
        Factor matrix used for update calculation.

    S : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph matrix.

    D : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph diagonal matrix.

    E : ndarray of shape (n_features, n_samples)
        A rectangle matrix whose elements are all one.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    Returns
    -------
    delta_H : ndarray of shape (n_components, n_samples)
        Factor matrix to be multiplied by H.
    """
    nominator = W.T @ (X / (_fill_zeros(W @ H))) + gamma * (H @ S)
    denominator = _fill_zeros(W.T @ E + beta * H + gamma * (H @ D))
    return np.divide(nominator, denominator)


def fit_multiplicative_update(
    X,
    W,
    H,
    Sg,
    Dg,
    Lg,
    Lc,
    E,
    e,
    tol: float = 1e-4,
    max_iter: int = 200,
    alpha: float = 1e-2,
    beta: float = 1e-1,
    gamma: float = 100,
    update_W: bool = True,
    update_H: bool = True,
    verbose: int = 0,
):
    """
    Run multiplicative updates on factor matrices
    to minimize the loss function.

    Parameters
    ----------
    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    W : ndarray of shape (n_features, n_components)
        Initial value for the factor matrix.

    H : ndarray of shape (n_components, n_samples)
        Initial value for the factor matrix.

    Sg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph matrix.

    Dg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph diagonal matrix.

    Lg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph Laplacian.

    Lc : ndarray of shape (n_samples, n_samples)
        Different-class nearest neighbor graph Laplacian.

    E : ndarray of shape (n_features, n_samples)
        A rectangle matrix whose elements are all one.

    e : ndarray of shape (n_samples, n_samples)
        e = 1 - I, where 1 is a matrix whose elements
        are all one and I is an identity matrix.

    tol : float
        Minimum relative improvement in the loss function
        over every 10 iterations.

    max_iter : int
        Maximum number of iterations that can be run.

    alpha : float
        Parameter for orthogonality regularization.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    update_W : bool
        Whether to update W during the iterations.

    update_H : bool
        Whether to update H during the iterations.

    verbose : int
        Verbosity level of the iterations.

    Returns
    -------
    W : ndarray of shape (n_features, n_components)
        Updated factor matrix after the iterations.

    H : ndarray of shape (n_components, n_samples)
        Updated factor matrix after the iterations.

    n_iter : int
        Number of iterations that have run.
    """
    start_time = time.time()
    Lc_power_minus_half = np.nan_to_num(fractional_matrix_power(Lc, -0.5))
    S = Lc_power_minus_half.T @ Sg @ Lc_power_minus_half
    D = Lc_power_minus_half.T @ Dg @ Lc_power_minus_half
    loss_init = mdnmf_loss(
        X=X, W=W, H=H, Lg=Lg, Lc=Lc, e=e, alpha=alpha, beta=beta, gamma=gamma
    )
    previous_loss = loss_init

    for n_iter in range(1, max_iter + 1):
        # update W
        if update_W:
            W = W * _grad_w(X, W, H, E, e, alpha)
            W[W < np.finfo(np.float64).eps] = 0.0
            W = np.abs(W)

        # update H
        if update_H:
            H = H * _grad_h(X, W, H, S, D, E, beta, gamma)
            H[H < np.finfo(np.float64).eps] = 0.0
            H = np.abs(H)

        # test convergence criterion every 10 iterations
        if tol > 0 and n_iter % 10 == 0:
            loss = mdnmf_loss(
                X=X, W=W, H=H, Lg=Lg, Lc=Lc, e=e, alpha=alpha, beta=beta, gamma=gamma
            )

            if verbose:
                iter_time = time.time()
                print(
                    f"Epoch {n_iter:02d} reached after "
                    f"{iter_time - start_time:.3f} seconds, loss: {loss:f}"
                )

            if (previous_loss - loss) / loss_init < tol:
                break
            previous_loss = loss

    # do not print if we have already printed in the convergence test
    if verbose and (tol == 0 or n_iter % 10 != 0):
        end_time = time.time()
        print(f"Epoch {n_iter:02d} reached after {end_time - start_time:.3f} seconds.")

    return W, H, n_iter


"""
Coordinate descent update for
Manifold Regularized Discriminative NMF.

Note the matrix dimensions are reversed in this section
(compared with conventional (n_samples, n_features)
we have (n_features, n_samples) instead), to be aligned
with most literature in NMF.
"""


def _has_converged(
    previous_value: Optional[float], current_value: float, tol: float = 1e-4
) -> bool:
    """
    Check whether the relative change is below the tolerance.

    Parameters
    ----------
    previous_value : float
        Previous target value.

    current_value : float
        Current target value.

    tol : float
        Tolerance threshold of the relative change.

    Returns
    -------
    has_converged : bool
        Whether the relative change is below the tolerance.
    """
    if previous_value is None:
        return False
    return abs((current_value - previous_value) / previous_value) < tol


def _calc_lambda(A, delta_A) -> float:
    """
    Find the supremum, i.e. the least upper bound:

        .. math::

            max\\{(A_{ij}) / (|\\Delta_{A_{ij}}|) | \\Delta_{A_{ij}} < 0\\}

    Parameters
    ----------
    A : array-like of shape (n, m)
        Target matrix.

    delta_A : array-like of shape (n, m)
        Gradient of loss function w.r.t. A.

    Returns
    -------
    lambda_ : float
        Supremum value as defined above.
    """
    lambda_ = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if delta_A[i][j] < 0 and A[i][j] / -delta_A[i][j] > lambda_:
                lambda_ = A[i][j] / -delta_A[i][j]
    return lambda_


def _calc_phi_p_w(
    theta: float, X, W, delta_W, W_H, delta_W_H, e_delta_W_T, alpha: float
):
    """
    Calculate the first derivative of phi_W, the perturbed loss function in W,
    with respect to theta:

        .. math::

            \\phi^\\prime_W(\\theta) =
            & (\\alpha tr (\\Delta_W e \\Delta_W^T)) \\theta \\\\
            & - \\sum_{i, j} \\frac{X_{ij} (\\Delta_W H)_{ij}}{(WH)_{ij}
            + (\\Delta_W H)_{ij}\\theta} + \\sum_{i, j} (\\Delta_W H)_{ij} \\\\
            & + \\alpha tr (W e \\Delta_W^T)

    Parameters
    ----------
    theta: float
        Step size in the gradient update.

    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    W : ndarray of shape (n_features, n_components)
        Factor matrix been calculated.

    delta_W : ndarray of shape (n_features, n_components)
        Gradient matrix of the loss function w.r.t. W.

    W_H : ndarray of shape (n_features, n_samples)
        W @ H, passed to avoid repetitive computations.

    delta_W_H : ndarray of shape (n_features, n_samples)
        delta_W @ H, passed to avoid repetitive computations.

    e_delta_W_T : ndarray of shape (n_components, n_features)
        e @ delta_W.T, passed to avoid repetitive computations.

    alpha : float
        Parameter for orthogonality regularization.

    Returns
    -------
    phi_p_w : float
        The first derivative of phi_W, the perturbed loss function in W.
    """
    term_1 = alpha * np.trace(delta_W @ e_delta_W_T) * theta
    term_2 = -np.sum(np.divide(X * delta_W_H, _fill_zeros(W_H + delta_W_H * theta)))
    term_3 = np.sum(delta_W_H)
    term_4 = alpha * np.trace(W @ e_delta_W_T)
    return term_1 + term_2 + term_3 + term_4


def _calc_phi_pp_w(theta: float, X, delta_W, W_H, delta_W_H, e_delta_W_T, alpha: float):
    """
    Calculate the second derivative of phi_W, the perturbed loss function in W,
    with respect to theta:

        .. math::

            \\phi^{\\prime\\prime}_W(\\theta) =
            & \\alpha tr (\\Delta_W e \\Delta_W^T) \\\\
            & + \\sum_{i, j} \\frac{X_{ij} (\\Delta_W H)^2_{ij}}{((WH)_{ij}
            + (\\Delta_W H)_{ij}\\theta)^2}

    Parameters
    ----------
    theta: float
        Step size in the gradient update.

    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    delta_W : ndarray of shape (n_features, n_components)
        Gradient matrix of the loss function w.r.t. W.

    W_H : ndarray of shape (n_features, n_samples)
        W @ H, passed to avoid repetitive computations.

    delta_W_H : ndarray of shape (n_features, n_samples)
        delta_W @ H, passed to avoid repetitive computations.

    e_delta_W_T : ndarray of shape (n_components, n_features)
        e @ delta_W.T, passed to avoid repetitive computations.

    alpha : float
        Parameter for orthogonality regularization.

    Returns
    -------
    phi_pp_w : float
        The second derivative of phi_W, the perturbed loss function in W.
    """
    term_1 = alpha * np.trace(delta_W @ e_delta_W_T)
    term_2 = np.sum(
        np.divide(
            X * np.square(delta_W_H), _fill_zeros(np.square(W_H + delta_W_H * theta))
        )
    )
    return term_1 + term_2


def _calc_phi_p_h(
    theta: float, X, H, delta_H, W_H, W_delta_H, L_delta_H_T, beta: float, gamma: float
):
    """
    Calculate the first derivative of phi_H, the perturbed loss function in H,
    with respect to theta:

        .. math::

            \\phi^\\prime_H(\\theta) =
            & (\\beta tr (\\Delta_H \\Delta_H^T)
            + \\gamma tr (\\Delta_H L \\Delta_H^T)) \\theta \\\\
            & - \\sum_{i, j} \\frac{X_{ij} (W \\Delta_H)_{ij}}{(WH)_{ij}
            + (W \\Delta_H)_{ij}\\theta} + \\sum_{i, j} (W \\Delta_H)_{ij} \\\\
            & + \\beta tr (H \\Delta_H^T) + \\gamma tr (H L \\Delta_H^T)

    Parameters
    ----------
    theta: float
        Step size in the gradient update.

    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    H : ndarray of shape (n_components, n_samples)
        Factor matrix been calculated.

    delta_H : ndarray of shape (n_components, n_samples)
        Gradient matrix of the loss function w.r.t. H.

    W_H : ndarray of shape (n_features, n_samples)
        W @ H, passed to avoid repetitive computations.

    W_delta_H : ndarray of shape (n_features, n_samples)
        W @ delta_H, passed to avoid repetitive computations.

    L_delta_H_T : ndarray of shape (n_samples, n_components)
        L @ delta_H.T, passed to avoid repetitive computations.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    Returns
    -------
    phi_p_h : float
        The first derivative of phi_H, the perturbed loss function in H.
    """
    term_1 = beta * np.trace(delta_H @ delta_H.T) * theta
    term_2 = gamma * np.trace(delta_H @ L_delta_H_T) * theta
    term_3 = -np.sum(np.divide(X * W_delta_H, _fill_zeros(W_H + W_delta_H * theta)))
    term_4 = np.sum(W_delta_H)
    term_5 = beta * np.trace(H @ delta_H.T)
    term_6 = gamma * np.trace(H @ L_delta_H_T)
    return term_1 + term_2 + term_3 + term_4 + term_5 + term_6


def _calc_phi_pp_h(
    theta: float, X, delta_H, W_H, W_delta_H, L_delta_H_T, beta: float, gamma: float
):
    """
    Calculate the second derivative of phi_H, the perturbed loss function in H,
    with respect to theta:

        .. math::

            \\phi^{\\prime\\prime}_H(\\theta) =
            & \\beta tr (\\Delta_H \\Delta_H^T)
            + \\gamma tr (\\Delta_H L \\Delta_H^T) \\\\
            & + \\sum_{i, j} \\frac{X_{ij} (W \\Delta_H)^2_{ij}}{((WH)_{ij}
            + (W \\Delta_H)_{ij}\\theta)^2}

    Parameters
    ----------
    theta: float
        Step size in the gradient update.

    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    delta_H : ndarray of shape (n_components, n_samples)
        Gradient matrix of the loss function w.r.t. H.

    W_H : ndarray of shape (n_features, n_samples)
        W @ H, passed to avoid repetitive computations.

    W_delta_H : ndarray of shape (n_features, n_samples)
        W @ delta_H, passed to avoid repetitive computations.

    L_delta_H_T : ndarray of shape (n_samples, n_components)
        L @ delta_H.T, passed to avoid repetitive computations.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    Returns
    -------
    phi_pp_h : float
        The second derivative of phi_H, the perturbed loss function in H.
    """
    term_1 = beta * np.trace(delta_H @ delta_H.T)
    term_2 = gamma * np.trace(delta_H @ L_delta_H_T)
    term_3 = np.sum(
        np.divide(
            X * np.square(W_delta_H), _fill_zeros(np.square(W_H + W_delta_H * theta))
        )
    )
    return term_1 + term_2 + term_3


def _coordinate_descent(
    _calc_phi_p: Callable[[float], float],
    _calc_phi_pp: Callable[[float], float],
    prev_theta: float,
    theta_tol: float,
    max_cd_iter: float,
    theta_p: float,
    optimizer_config: Dict,
):
    """
    Calculate the step size, theta, for the gradient update.

    Parameters
    ----------
    _calc_phi_p : Callable[[float], float]
        Callable for calculating first derivate of the loss function.

    _calc_phi_p : Callable[[float], float]
        Callable for calculating second derivate of the loss function.

    prev_theta : float
        Previous theta value, used for initial value calculation for this step.

    theta_tol : float
        Tolerance on the relative change in theta convergence.

    max_cd_iter : int
        Maximum number of iterations to be run for coordinate descent.

    theta_p : float
        Upper bound heuristic on theta.

    optimizer_config : Dict
        Configuration for the optimizer used for coordinate descent.

    Returns
    -------
    theta : float
        Step size in the gradient update.
    """
    # parse universal optimizer configurations
    optimizer = optimizer_config.get("name")
    theta_min = optimizer_config.get("theta_min", 0.0)
    theta_max = optimizer_config.get("theta_max", 10.0)

    # parse optimizer-specific configurations
    if optimizer in ("momentum", "nesterov"):
        gamma = optimizer_config.get("gamma", 0.9)
        eta = optimizer_config.get("eta", 0.1)
        v = 0
    elif optimizer in ("adam", "adamax", "nadam", "amsgrad"):
        beta_1 = optimizer_config.get("beta_1", 0.9)
        beta_2 = optimizer_config.get("beta_2", 0.999)
        epsilon = optimizer_config.get("epsilon", 1e-8)
        eta = optimizer_config.get("eta", 0.002)
        m = v = m_hat = v_hat = 0
    else:  # default to use SOR Newton-Ralphson
        omega = optimizer_config.get("omega", 0.1)

    # initialize search variables
    iter_idx = 1
    previous_theta = None
    theta = (theta_min + prev_theta) / 2

    # search for optimal step size
    while iter_idx <= max_cd_iter and not _has_converged(
        previous_theta, theta, theta_tol
    ):
        previous_theta = theta

        if optimizer == "momentum":
            phi_p = _calc_phi_p(theta)
            v = gamma * v + eta * phi_p
            theta -= v
        elif optimizer == "nesterov":
            phi_p = _calc_phi_p(theta - gamma * v)
            v = gamma * v + eta * phi_p
            theta -= v
        elif optimizer in ("adam", "adamax", "nadam", "amsgrad"):
            phi_p = _calc_phi_p(theta)
            m = beta_1 * m + (1 - beta_1) * phi_p
            v = beta_2 * v + (1 - beta_2) * phi_p * phi_p
            if optimizer == "adam":
                m_hat = m / (1 - beta_1 ** iter_idx)
                v_hat = v / (1 - beta_2 ** iter_idx)
                theta -= eta * m_hat / (np.sqrt(v_hat) + epsilon)
            elif optimizer == "adamax":
                m_hat = m / (1 - beta_1 ** iter_idx)
                u = max(beta_2 * v, abs(phi_p))
                theta -= eta * m_hat / u
            elif optimizer == "nadam":
                m_hat = m / (1 - beta_1 ** iter_idx) + (1 - beta_1) * phi_p / (
                    1 - beta_1 ** iter_idx
                )
                v_hat = v / (1 - beta_2 ** iter_idx)
                theta -= eta * m_hat / (np.sqrt(v_hat) + epsilon)
            else:
                v_hat = max(v_hat, v)
                theta -= eta * m / (np.sqrt(v_hat) + epsilon)
        else:  # default to use SOR Newton-Ralphson
            phi_p = _calc_phi_p(theta)
            phi_pp = _calc_phi_pp(theta)
            theta = (1 - omega) * theta + omega * (theta - phi_p / phi_pp)

        iter_idx += 1

    return min(max(min(theta, theta_p), theta_min), theta_max)


def fit_coordinate_descent(
    X,
    W,
    H,
    Sg,
    Dg,
    Lg,
    Lc,
    E,
    e,
    tol: float = 1e-4,
    cd_tol: float = 1e-4,
    max_iter: int = 200,
    max_cd_iter: int = 200,
    optimizer_config: Optional[Dict] = None,
    alpha: float = 1e-2,
    beta: float = 1e-1,
    gamma: float = 100,
    update_W: bool = True,
    update_H: bool = True,
    verbose: int = 0,
):
    """
    Run coordinate descent updates on factor matrices
    to minimize the loss function.

    Parameters
    ----------
    X : array-like of shape (n_features, n_samples)
        Matrix to be factorized.

    W : ndarray of shape (n_features, n_components)
        Initial value for the factor matrix.

    H : ndarray of shape (n_components, n_samples)
        Initial value for the factor matrix.

    Sg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph matrix.

    Dg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph diagonal matrix.

    Lg : ndarray of shape (n_samples, n_samples)
        Same-class nearest neighbor graph Laplacian.

    Lc : ndarray of shape (n_samples, n_samples)
        Different-class nearest neighbor graph Laplacian.

    E : ndarray of shape (n_features, n_samples)
        A rectangle matrix whose elements are all one.

    e : ndarray of shape (n_samples, n_samples)
        e = 1 - I, where 1 is a matrix whose elements
        are all one and I is an identity matrix.

    tol : float
        Minimum relative improvement in the loss function
        over every 10 iterations.

    cd_tol : float
        Tolerance on the relative change in theta convergence.

    max_iter : int
        Maximum number of iterations that can be run.

    max_cd_iter : int
        Maximum number of iterations that can be run
        during each step size search.

    optimizer_config : Dict
        Configuration for the optimizer.

    alpha : float
        Parameter for orthogonality regularization.

    beta : float
        Parameter for Tiknohov regularization.

    gamma : float
        Parameter for discriminative term.

    update_W : bool
        Whether to update W during the iterations.

    update_H : bool
        Whether to update H during the iterations.

    verbose : int
        Verbosity level of the iterations.

    Returns
    -------
    W : ndarray of shape (n_features, n_components)
        Updated factor matrix after the iterations.

    H : ndarray of shape (n_components, n_samples)
        Updated factor matrix after the iterations.

    n_iter : int
        Number of iterations that have run.
    """
    X = check_array(X, accept_sparse="csr")
    optimizer_config = optimizer_config if optimizer_config else {}
    start_time = time.time()
    if update_H:
        Lc_power_minus_half = np.nan_to_num(fractional_matrix_power(Lc, -0.5))
        S = Lc_power_minus_half.T @ Sg @ Lc_power_minus_half
        D = Lc_power_minus_half.T @ Dg @ Lc_power_minus_half
        L = D - S
    loss_init = mdnmf_loss(
        X=X, W=W, H=H, Lg=Lg, Lc=Lc, e=e, alpha=alpha, beta=beta, gamma=gamma
    )
    previous_loss = loss_init
    theta_W = theta_H = 1.0
    W_H = W @ H

    for n_iter in range(1, max_iter + 1):
        # Update W
        if update_W:
            delta_W = W * _grad_w(X, W, H, E, e, alpha) - W
            theta_p = 0.01 + 0.99 * _calc_lambda(W, delta_W)
            delta_W_H = _fill_zeros(delta_W @ H)
            e_delta_W_T = e @ delta_W.T
            theta_W = _coordinate_descent(
                lambda theta: _calc_phi_p_w(
                    theta, X, W, delta_W, W_H, delta_W_H, e_delta_W_T, alpha
                ),
                lambda theta: _calc_phi_pp_w(
                    theta, X, delta_W, W_H, delta_W_H, e_delta_W_T, alpha
                ),
                theta_W,
                cd_tol,
                max_cd_iter,
                theta_p,
                optimizer_config,
            )
            W = W + theta_W * delta_W
            W[W < np.finfo(np.float64).eps] = 0.0
            W = np.abs(W)
            W_H = W @ H

        # Update H
        if update_H:
            delta_H = (
                H * _grad_h(X=X, W=W, H=H, D=D, S=S, E=E, beta=beta, gamma=gamma) - H
            )
            theta_p = 0.01 + 0.99 * _calc_lambda(H, delta_H)
            W_delta_H = _fill_zeros(W @ delta_H)
            L_delta_H_T = L @ delta_H.T
            theta_H = _coordinate_descent(
                lambda theta: _calc_phi_p_h(
                    theta, X, H, delta_H, W_H, W_delta_H, L_delta_H_T, beta, gamma
                ),
                lambda theta: _calc_phi_pp_h(
                    theta, X, delta_H, W_H, W_delta_H, L_delta_H_T, beta, gamma
                ),
                theta_H,
                cd_tol,
                max_cd_iter,
                theta_p,
                optimizer_config,
            )
            H = H + theta_H * delta_H
            H[H < np.finfo(np.float64).eps] = 0.0
            H = np.abs(H)
            W_H = W @ H

        # test convergence criterion every 10 iterations
        if tol > 0 and n_iter % 10 == 0:
            loss = mdnmf_loss(
                X=X, W=W, H=H, Lg=Lg, Lc=Lc, e=e, alpha=alpha, beta=beta, gamma=gamma
            )

            if verbose:
                iter_time = time.time()
                print(
                    f"Epoch {n_iter:02d} reached after "
                    f"{iter_time - start_time:.3f} seconds, loss: {loss:f}"
                )

            if (previous_loss - loss) / loss_init < tol:
                break
            previous_loss = loss

    # do not print if we have already printed in the convergence test
    if verbose and (tol == 0 or n_iter % 10 != 0):
        end_time = time.time()
        print(f"Epoch {n_iter:02d} reached after {end_time - start_time:.3f} seconds.")

    return W, H, n_iter


class MDNMF(_ClassNamePrefixFeaturesOutMixin, TransformerMixin, BaseEstimator):
    """
    Manifold-regularised Discriminative Non-Negative Matrix Factorization (MDNMF).

    Find two non-negative matrices (W, H) whose product approximates the non-
    negative matrix X, while preserving the geometry of the data space with
    manifold regularization and incorporating discriminative information by
    maximizing the inter-class margins. This factorization can be used for
    dimensionality reduction, source separation or topic extraction.

    The objective function is:

        .. math::

            KL(X, WH)

            + 0.5 * alpha * tr(W e W^T)

            + 0.5 * beta * tr(H H^T)

            + 0.5 * gamma * tr(H (L_c^{-0.5}))^T L_g L_c^{-0.5} H^T)

    Where:

    :math:`KL(A, B) = \\sum_{i,j} (A_{ij} \\log \\frac{A_{ij}}{B_{ij}}
        - A_{ij} + B{ij})` (Kullback-Leibler divergence)

    :math:`L_c` is the different-class nearest neighbor graph Laplacian while
    :math:`L_g` is for samples of the same-class.

    The objective function is minimized with an alternating minimization of W
    and H.

    Read more in the :ref:`User Guide <MDNMF>`.

    Parameters
    ----------
    n_components : int, default=None
        Number of components, if n_components is not set all features
        are kept.

    solver : {'cd', 'mu'}, default='cd'
        Numerical solver to use:
        'cd' is a Coordinate Descent solver.
        'mu' is a Multiplicative Update solver.

    tol : float, default=1e-4
        Minimum relative improvement in the loss function
        over every 10 iterations.

    cd_tol : float, default=1e-4
        Tolerance on the relative change in theta convergence.

    max_iter : int, default=200
        Maximum number of iterations before timing out.

    max_cd_iter : int, default=200
        Maximum number of iterations that can be run
        during each step size search.

    optimizer_config : Dict
        Configuration for the optimizer used in coordinate descent.

    random_state : int, RandomState instance or None, default=None
        Used for initialisation (when ``init`` == 'nndsvdar' or
        'random'), and in Coordinate Descent. Pass an int for reproducible
        results across multiple function calls.
        See :term:`Glossary <random_state>`.

    alpha : float, default=1e-2
        Parameter for orthogonality regularization.

    beta : float, default=1e-1
        Parameter for Tiknohov regularization.

    gamma : float, default=100
        Parameter for discriminative term.

    k1 : int, default=1
        Number of nearest neighbors to consider for the same class
        nearest neighbors graph laplacian.

    k2 : int, default=1
        Number of nearest neighbors to consider for the different class
        nearest neighbors graph laplacian.

    verbose : int, default=0
        Verbosity level of the iterations.

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Factorization matrix, sometimes called 'dictionary'.

    n_components_ : int
        The number of components. It is same as the `n_components` parameter
        if it was given. Otherwise, it will be same as the number of
        features.

    reconstruction_err_ : float
        Frobenius norm of the matrix difference, or beta-divergence, between
        the training data ``X`` and the reconstructed data ``WH`` from
        the fitted model.

    n_iter_ : int
        Actual number of iterations.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`, )
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    DictionaryLearning : Find a dictionary that sparsely encodes data.
    MiniBatchSparsePCA : Mini-batch Sparse Principal Components Analysis.
    PCA : Principal component analysis.
    SparseCoder : Find a sparse representation of data from a fixed,
        precomputed dictionary.
    SparsePCA : Sparse Principal Components Analysis.
    TruncatedSVD : Dimensionality reduction using truncated SVD.
    NMF : Non-negative Matrix Factorization without geometric and
        discriminative constraints.

    References
    ----------
    Guan, Naiyang, et al. "Manifold regularized discriminative nonnegative
    matrix factorization with fast gradient descent."
    IEEE Transactions on Image Processing 20.7 (2011): 2030-2048.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> y = [0, 1, 0, 0, 1, 2]
    >>> from sklearn.decomposition import MDNMF
    >>> model = MDNMF(n_components=2, random_state=0)
    >>> W = model.fit_transform(np.c_[X, y])
    >>> H = model.components_
    """

    def __init__(
        self,
        n_components: Optional[int] = None,
        solver="cd",
        tol: float = 1e-4,
        cd_tol: float = 1e-4,
        max_iter: int = 200,
        max_cd_iter: int = 200,
        optimizer_config: Optional[Dict] = None,
        random_state: Optional[int] = None,
        alpha: float = 1e-2,
        beta: float = 1e-1,
        gamma: float = 100,
        k1: int = 1,
        k2: int = 1,
        verbose: Optional[int] = 0,
    ):
        self.n_components = n_components
        self.solver = solver
        self.tol = tol
        self.cd_tol = cd_tol
        self.max_iter = max_iter
        self.max_cd_iter = max_cd_iter
        self.optimizer_config = optimizer_config
        self.random_state = random_state
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.k1 = k1
        self.k2 = k2
        self.verbose = verbose

    def fit(self, X, y=None, W=None, H=None):
        """
        Learn a MDNMF model for the data X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features + 1)
            Training vector, where `n_samples` is the number of samples,
            `n_features` is the number of features and `+1` for labels.

        y : Ignored
            Not used, present for API consistency by convention.

        W : ndarray of shape (n_features, n_components)
            Initial value for the factor matrix.

        H : ndarray of shape (n_components, n_samples)
            Initial value for the factor matrix.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        X, y = self._check_x_y(X)

        with config_context(assume_finite=True):
            W, H, self.n_iter_ = self._mdnmf(X, y, W=W, H=H, update_H=True)

        self.reconstruction_err_ = _beta_divergence(
            X, W, H, "kullback-leibler", square_root=True
        )
        self.n_features_in_ = X.shape[1] + 1
        self.n_components_ = H.shape[0]
        self.components_ = H

        return self

    def transform(self, X):
        """
        Transform the data X according to the fitted MDNMF model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features + 1)
            Training vector, where `n_samples` is the number of samples,
            `n_features` is the number of features and `+1` for labels.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        check_is_fitted(self)
        X, y = self._check_x_y(X, reset=False)

        with config_context(assume_finite=True):
            W, *_ = self._mdnmf(
                X, y, H=self.components_.astype(X.dtype, copy=False), update_H=False
            )
        return W

    def fit_transform(self, X, y=None, W=None, H=None):
        """
        Learn a MDNMF model for the data X and returns the transformed data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features + 1)
            Training vector, where `n_samples` is the number of samples,
            `n_features` is the number of features and `+1` for labels.

        y : Ignored
            Not used, present for API consistency by convention.

        W : ndarray of shape (n_features, n_components)
            Initial value for the factor matrix.

        H : ndarray of shape (n_components, n_samples)
            Initial value for the factor matrix.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.

        H : ndarray of shape (n_components, n_features)
            Factorization matrix, sometimes called 'dictionary'.

        n_iter_ : int
            Actual number of iterations.
        """
        self.fit(X, y, W, H)
        return self.transform(X)

    @property
    def _n_features_out(self):
        """
        Number of transformed output features.
        """
        return self.components_.shape[0]

    def _more_tags(self):
        return {"requires_positive_X": True}

    def _check_params(self, X):
        # n_components
        self._n_components = self.n_components
        if self._n_components is None:
            self._n_components = X.shape[1]
        if (
            not isinstance(self._n_components, numbers.Integral)
            or self._n_components <= 0
        ):
            raise ValueError(
                "Number of components must be a positive integer; got "
                f"(n_components={self._n_components!r})"
            )

        # max_iter
        if not isinstance(self.max_iter, numbers.Integral) or self.max_iter < 0:
            raise ValueError(
                "Maximum number of iterations must be a positive "
                f"integer; got (max_iter={self.max_iter!r})"
            )

        # tol
        if not isinstance(self.tol, numbers.Number) or self.tol < 0:
            raise ValueError(
                "Tolerance for stopping criteria must be positive; got "
                f"(tol={self.tol!r})"
            )

        # solver
        allowed_solver = ("cd", "mu")
        if self.solver not in allowed_solver:
            raise ValueError(
                f"Invalid solver parameter: got {self.solver!r} instead of one of "
                f"{allowed_solver}"
            )

        if self.solver == "cd":
            # cd_tol
            if not isinstance(self.cd_tol, numbers.Number) or self.cd_tol < 0:
                raise ValueError(
                    "Tolerance for coordinate descent stopping criteria "
                    f"must be positive; got (tol={self.cd_tol!r})"
                )

        return self

    def _check_w_h(self, X, W, H, update_H):
        # check W and H, or initialize them
        n_samples, n_features = X.shape
        avg = np.sqrt(X.mean() / self._n_components)
        rng = check_random_state(self.random_state)
        if update_H:
            H = avg * rng.randn(self._n_components, n_features).astype(
                X.dtype, copy=False
            )
        else:
            _check_init(H, (self._n_components, n_features), "NMF (input H)")
            if H.dtype != X.dtype:
                raise TypeError(
                    f"H should have the same dtype as X. Got H.dtype = {H.dtype}."
                )
        W = avg * rng.randn(n_samples, self._n_components).astype(X.dtype, copy=False)
        np.abs(W, out=W)
        return W, H

    def _check_x_y(self, X_y, reset=True):
        X_y = self._validate_data(X_y, dtype=[np.float64, np.float32], reset=reset)
        if X_y.shape[1] == 1:
            raise ValueError(
                "Please provide more than 1 feature(s) as the last feature will be used"
                " as discriminative information."
            )
        return X_y[:, :-1], X_y[:, -1]

    def _mdnmf(self, X, y=None, W=None, H=None, update_H=True):
        """
        Learn a MDNMF model for the data X and returns the transformed data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Data matrix to be decomposed.

        y : Ignored.

        W : array-like of shape (n_samples, n_components)
            Initial guess for the solution.

        H : array-like of shape (n_components, n_features)
            Initial guess for the solution.
            If update_H=False, it is used as a constant, to solve for W only.

        update_H : bool, default=True
            If True, both W and H will be estimated from initial guesses,
            otherwise only W will be estimated.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.

        H : ndarray of shape (n_components, n_features)
            Factorization matrix, sometimes called 'dictionary'.

        n_iter_ : int
            Actual number of iterations.
        """
        check_non_negative(X, "NMF (input X)")

        # check parameters
        self._check_params(X)

        # initialize or check W and H
        W, H = self._check_w_h(X, W, H, update_H)

        # initialize intermediate quantities
        Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, self.k1, self.k2)
        e = np.ones((self._n_components, self._n_components)) - np.eye(
            self._n_components
        )
        X_fit, W_fit, H_fit, update_W_fit = X.T, H.T, W.T, update_H
        E = np.ones(X_fit.shape)

        if self.solver == "cd":
            W_fit, H_fit, n_iter = fit_coordinate_descent(
                X=X_fit,
                W=W_fit,
                H=H_fit,
                Sg=Sg,
                Dg=Dg,
                Lg=Lg,
                Lc=Lc,
                E=E,
                e=e,
                tol=self.tol,
                cd_tol=self.cd_tol,
                max_iter=self.max_iter,
                max_cd_iter=self.max_cd_iter,
                optimizer_config=self.optimizer_config,
                alpha=self.alpha,
                beta=self.beta,
                gamma=self.gamma,
                update_W=update_W_fit,
                verbose=self.verbose,
            )
        elif self.solver == "mu":
            W_fit, H_fit, n_iter = fit_multiplicative_update(
                X=X_fit,
                W=W_fit,
                H=H_fit,
                Sg=Sg,
                Dg=Dg,
                Lg=Lg,
                Lc=Lc,
                E=E,
                e=e,
                tol=self.tol,
                max_iter=self.max_iter,
                alpha=self.alpha,
                beta=self.beta,
                gamma=self.gamma,
                update_W=update_W_fit,
                verbose=self.verbose,
            )
        else:
            raise ValueError("Invalid solver parameter '%s'." % self.solver)
        W, H = H_fit.T, W_fit.T

        if n_iter == self.max_iter and self.tol > 0:
            warnings.warn(
                f"Maximum number of iterations {self.max_iter} reached. Increase "
                "it to improve convergence.",
                ConvergenceWarning,
            )

        return W, H, n_iter
