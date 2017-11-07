"""GraphLasso: sparse inverse covariance estimation with an l1-penalized
estimator.
"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause
# Copyright: INRIA
import numpy as np
from ..utils import deprecated
from .graphical_lasso_ import GraphicalLasso, GraphicalLassoCV, graphical_lasso


# The g-lasso algorithm

@deprecated("The 'graph_lasso' was renamed to 'graphical_lasso' "
            "in version 0.20")
def graph_lasso(emp_cov, alpha, cov_init=None, mode='cd', tol=1e-4,
                enet_tol=1e-4, max_iter=100, verbose=False,
                return_costs=False, eps=np.finfo(np.float64).eps,
                return_n_iter=False):
    """l1-penalized covariance estimator

    .. deprecated:: 0.20
        The `graph_lasso` was renamed to `graphical_lasso` in version 0.20

    Read more in the :ref:`User Guide <sparse_inverse_covariance>`.

    Parameters
    ----------
    emp_cov : 2D ndarray, shape (n_features, n_features)
        Empirical covariance from which to compute the covariance estimate.

    alpha : positive float
        The regularization parameter: the higher alpha, the more
        regularization, the sparser the inverse covariance.

    cov_init : 2D array (n_features, n_features), optional
        The initial guess for the covariance.

    mode : {'cd', 'lars'}
        The Lasso solver to use: coordinate descent or LARS. Use LARS for
        very sparse underlying graphs, where p > n. Elsewhere prefer cd
        which is more numerically stable.

    tol : positive float, optional
        The tolerance to declare convergence: if the dual gap goes below
        this value, iterations are stopped.

    enet_tol : positive float, optional
        The tolerance for the elastic net solver used to calculate the descent
        direction. This parameter controls the accuracy of the search direction
        for a given column update, not of the overall parameter estimate. Only
        used for mode='cd'.

    max_iter : integer, optional
        The maximum number of iterations.

    verbose : boolean, optional
        If verbose is True, the objective function and dual gap are
        printed at each iteration.

    return_costs : boolean, optional
        If return_costs is True, the objective function and dual gap
        at each iteration are returned.

    eps : float, optional
        The machine-precision regularization in the computation of the
        Cholesky diagonal factors. Increase this for very ill-conditioned
        systems.

    return_n_iter : bool, optional
        Whether or not to return the number of iterations.

    Returns
    -------
    covariance : 2D ndarray, shape (n_features, n_features)
        The estimated covariance matrix.

    precision : 2D ndarray, shape (n_features, n_features)
        The estimated (sparse) precision matrix.

    costs : list of (objective, dual_gap) pairs
        The list of values of the objective function and the dual gap at
        each iteration. Returned only if return_costs is True.

    n_iter : int
        Number of iterations. Returned only if `return_n_iter` is set to True.

    See Also
    --------
    GraphLasso, GraphLassoCV

    Notes
    -----
    The algorithm employed to solve this problem is the GLasso algorithm,
    from the Friedman 2008 Biostatistics paper. It is the same algorithm
    as in the R `glasso` package.

    One possible difference with the `glasso` R package is that the
    diagonal coefficients are not penalized.

    """
    return graphical_lasso(emp_cov, alpha, cov_init, mode, tol,
                           enet_tol, max_iter, verbose, return_costs,
                           eps, return_n_iter)


@deprecated("The 'GraphLasso' was renamed to 'GraphicalLasso' "
            "in version 0.20")
class GraphLasso(GraphicalLasso):
    """Sparse inverse covariance estimation with an l1-penalized estimator.

    This class implements the Graphical Lasso algorithm.

    .. deprecated:: 0.20
        The `GraphLasso` was renamed to `GraphicalLasso` in version 0.20

    Read more in the :ref:`User Guide <sparse_inverse_covariance>`.

    Parameters
    ----------
    alpha : positive float, default 0.01
        The regularization parameter: the higher alpha, the more
        regularization, the sparser the inverse covariance.

    mode : {'cd', 'lars'}, default 'cd'
        The Lasso solver to use: coordinate descent or LARS. Use LARS for
        very sparse underlying graphs, where p > n. Elsewhere prefer cd
        which is more numerically stable.

    tol : positive float, default 1e-4
        The tolerance to declare convergence: if the dual gap goes below
        this value, iterations are stopped.

    enet_tol : positive float, optional
        The tolerance for the elastic net solver used to calculate the descent
        direction. This parameter controls the accuracy of the search direction
        for a given column update, not of the overall parameter estimate. Only
        used for mode='cd'.

    max_iter : integer, default 100
        The maximum number of iterations.

    verbose : boolean, default False
        If verbose is True, the objective function and dual gap are
        plotted at each iteration.

    assume_centered : boolean, default False
        If True, data are not centered before computation.
        Useful when working with data whose mean is almost, but not exactly
        zero.
        If False, data are centered before computation.

    Attributes
    ----------
    covariance_ : array-like, shape (n_features, n_features)
        Estimated covariance matrix

    precision_ : array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.

    n_iter_ : int
        Number of iterations run.

    See Also
    --------
    graph_lasso, GraphLassoCV
    """


@deprecated("The 'GraphLassoCV' was renamed to 'GraphicalLassoCV' "
            "in version 0.20")
class GraphLassoCV(GraphicalLassoCV):
    """Sparse inverse covariance w/ cross-validated choice of the l1 penalty

    This class implements the Graphical Lasso algorithm.

    .. deprecated:: 0.20
        The `GraphLassoCV` was renamed to `GraphicalLassoCV` in version 0.20

    Read more in the :ref:`User Guide <sparse_inverse_covariance>`.

    Parameters
    ----------
    alphas : integer, or list positive float, optional
        If an integer is given, it fixes the number of points on the
        grids of alpha to be used. If a list is given, it gives the
        grid to be used. See the notes in the class docstring for
        more details.

    n_refinements : strictly positive integer
        The number of times the grid is refined. Not used if explicit
        values of alphas are passed.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - An object to be used as a cross-validation generator.
        - An iterable yielding train/test splits.

        For integer/None inputs :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    tol : positive float, optional
        The tolerance to declare convergence: if the dual gap goes below
        this value, iterations are stopped.

    enet_tol : positive float, optional
        The tolerance for the elastic net solver used to calculate the descent
        direction. This parameter controls the accuracy of the search direction
        for a given column update, not of the overall parameter estimate. Only
        used for mode='cd'.

    max_iter : integer, optional
        Maximum number of iterations.

    mode : {'cd', 'lars'}
        The Lasso solver to use: coordinate descent or LARS. Use LARS for
        very sparse underlying graphs, where number of features is greater
        than number of samples. Elsewhere prefer cd which is more numerically
        stable.

    n_jobs : int, optional
        number of jobs to run in parallel (default 1).

    verbose : boolean, optional
        If verbose is True, the objective function and duality gap are
        printed at each iteration.

    assume_centered : Boolean
        If True, data are not centered before computation.
        Useful when working with data whose mean is almost, but not exactly
        zero.
        If False, data are centered before computation.

    Attributes
    ----------
    covariance_ : numpy.ndarray, shape (n_features, n_features)
        Estimated covariance matrix.

    precision_ : numpy.ndarray, shape (n_features, n_features)
        Estimated precision matrix (inverse covariance).

    alpha_ : float
        Penalization parameter selected.

    cv_alphas_ : list of float
        All penalization parameters explored.

    grid_scores_ : 2D numpy.ndarray (n_alphas, n_folds)
        Log-likelihood score on left-out data across folds.

    n_iter_ : int
        Number of iterations run for the optimal alpha.

    See Also
    --------
    graph_lasso, GraphLasso

    Notes
    -----
    The search for the optimal penalization parameter (alpha) is done on an
    iteratively refined grid: first the cross-validated scores on a grid are
    computed, then a new refined grid is centered around the maximum, and so
    on.

    One of the challenges which is faced here is that the solvers can
    fail to converge to a well-conditioned estimate. The corresponding
    values of alpha then come out as missing values, but the optimum may
    be close to these missing values.
    """
