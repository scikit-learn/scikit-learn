# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Gael Varoquaux <gael.varoquaux@inria.fr>
#
# License: BSD 3 clause

import sys
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import sparse
from joblib import Parallel, delayed, effective_n_jobs

from ._base import LinearModel, _pre_fit
from ..base import RegressorMixin, MultiOutputMixin
from ._base import _preprocess_data
from ..utils import check_array, check_X_y
from ..utils.validation import check_random_state
from ..model_selection import check_cv
from ..utils.extmath import safe_sparse_dot
from ..utils.fixes import _joblib_parallel_args
from ..utils.validation import check_is_fitted
from ..utils.validation import column_or_1d

from . import _cd_fast as cd_fast


###############################################################################
# Paths functions

def _alpha_grid(X, y, Xy=None, l1_ratio=1.0, fit_intercept=True,
                eps=1e-3, n_alphas=100, normalize=False, copy_X=True):
    """ Compute the grid of alpha values for elastic net parameter search

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data. Pass directly as Fortran-contiguous data to avoid
        unnecessary memory duplication

    y : ndarray, shape (n_samples,)
        Target values

    Xy : array-like, optional
        Xy = np.dot(X.T, y) that can be precomputed.

    l1_ratio : float
        The elastic net mixing parameter, with ``0 < l1_ratio <= 1``.
        For ``l1_ratio = 0`` the penalty is an L2 penalty. (currently not
        supported) ``For l1_ratio = 1`` it is an L1 penalty. For
        ``0 < l1_ratio <1``, the penalty is a combination of L1 and L2.

    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``

    n_alphas : int, optional
        Number of alphas along the regularization path

    fit_intercept : boolean, default True
        Whether to fit an intercept or not

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.
    """
    if l1_ratio == 0:
        raise ValueError("Automatic alpha grid generation is not supported for"
                         " l1_ratio=0. Please supply a grid by providing "
                         "your estimator with the appropriate `alphas=` "
                         "argument.")
    n_samples = len(y)

    sparse_center = False
    if Xy is None:
        X_sparse = sparse.isspmatrix(X)
        sparse_center = X_sparse and (fit_intercept or normalize)
        X = check_array(X, 'csc',
                        copy=(copy_X and fit_intercept and not X_sparse))
        if not X_sparse:
            # X can be touched inplace thanks to the above line
            X, y, _, _, _ = _preprocess_data(X, y, fit_intercept,
                                             normalize, copy=False)
        Xy = safe_sparse_dot(X.T, y, dense_output=True)

        if sparse_center:
            # Workaround to find alpha_max for sparse matrices.
            # since we should not destroy the sparsity of such matrices.
            _, _, X_offset, _, X_scale = _preprocess_data(X, y, fit_intercept,
                                                          normalize,
                                                          return_mean=True)
            mean_dot = X_offset * np.sum(y)

    if Xy.ndim == 1:
        Xy = Xy[:, np.newaxis]

    if sparse_center:
        if fit_intercept:
            Xy -= mean_dot[:, np.newaxis]
        if normalize:
            Xy /= X_scale[:, np.newaxis]

    alpha_max = (np.sqrt(np.sum(Xy ** 2, axis=1)).max() /
                 (n_samples * l1_ratio))

    if alpha_max <= np.finfo(float).resolution:
        alphas = np.empty(n_alphas)
        alphas.fill(np.finfo(float).resolution)
        return alphas

    return np.logspace(np.log10(alpha_max * eps), np.log10(alpha_max),
                       num=n_alphas)[::-1]


def lasso_path(X, y, eps=1e-3, n_alphas=100, alphas=None,
               precompute='auto', Xy=None, copy_X=True, coef_init=None,
               verbose=False, return_n_iter=False, positive=False, **params):
    """Compute Lasso path with coordinate descent

    The Lasso optimization function varies for mono and multi-outputs.

    For mono-output tasks it is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

    For multi-output tasks it is::

        (1 / (2 * n_samples)) * ||Y - XW||^2_Fro + alpha * ||W||_21

    Where::

        ||W||_21 = \\sum_i \\sqrt{\\sum_j w_{ij}^2}

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <lasso>`.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data. Pass directly as Fortran-contiguous data to avoid
        unnecessary memory duplication. If ``y`` is mono-output then ``X``
        can be sparse.

    y : ndarray, shape (n_samples,), or (n_samples, n_outputs)
        Target values

    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : ndarray, optional
        List of alphas where to compute the models.
        If ``None`` alphas are set automatically

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument.

    Xy : array-like, optional
        Xy = np.dot(X.T, y) that can be precomputed. It is useful
        only when the Gram matrix is precomputed.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    coef_init : array, shape (n_features, ) | None
        The initial values of the coefficients.

    verbose : bool or integer
        Amount of verbosity.

    return_n_iter : bool
        whether to return the number of iterations or not.

    positive : bool, default False
        If set to True, forces coefficients to be positive.
        (Only allowed when ``y.ndim == 1``).

    **params : kwargs
        keyword arguments passed to the coordinate descent solver.

    Returns
    -------
    alphas : array, shape (n_alphas,)
        The alphas along the path where models are computed.

    coefs : array, shape (n_features, n_alphas) or \
            (n_outputs, n_features, n_alphas)
        Coefficients along the path.

    dual_gaps : array, shape (n_alphas,)
        The dual gaps at the end of the optimization for each alpha.

    n_iters : array-like, shape (n_alphas,)
        The number of iterations taken by the coordinate descent optimizer to
        reach the specified tolerance for each alpha.

    Notes
    -----
    For an example, see
    :ref:`examples/linear_model/plot_lasso_coordinate_descent_path.py
    <sphx_glr_auto_examples_linear_model_plot_lasso_coordinate_descent_path.py>`.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.

    Note that in certain cases, the Lars solver may be significantly
    faster to implement this functionality. In particular, linear
    interpolation can be used to retrieve model coefficients between the
    values output by lars_path

    Examples
    --------

    Comparing lasso_path and lars_path with interpolation:

    >>> X = np.array([[1, 2, 3.1], [2.3, 5.4, 4.3]]).T
    >>> y = np.array([1, 2, 3.1])
    >>> # Use lasso_path to compute a coefficient path
    >>> _, coef_path, _ = lasso_path(X, y, alphas=[5., 1., .5])
    >>> print(coef_path)
    [[0.         0.         0.46874778]
     [0.2159048  0.4425765  0.23689075]]

    >>> # Now use lars_path and 1D linear interpolation to compute the
    >>> # same path
    >>> from sklearn.linear_model import lars_path
    >>> alphas, active, coef_path_lars = lars_path(X, y, method='lasso')
    >>> from scipy import interpolate
    >>> coef_path_continuous = interpolate.interp1d(alphas[::-1],
    ...                                             coef_path_lars[:, ::-1])
    >>> print(coef_path_continuous([5., 1., .5]))
    [[0.         0.         0.46915237]
     [0.2159048  0.4425765  0.23668876]]


    See also
    --------
    lars_path
    Lasso
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode
    """
    return enet_path(X, y, l1_ratio=1., eps=eps, n_alphas=n_alphas,
                     alphas=alphas, precompute=precompute, Xy=Xy,
                     copy_X=copy_X, coef_init=coef_init, verbose=verbose,
                     positive=positive, return_n_iter=return_n_iter, **params)


def enet_path(X, y, l1_ratio=0.5, eps=1e-3, n_alphas=100, alphas=None,
              precompute='auto', Xy=None, copy_X=True, coef_init=None,
              verbose=False, return_n_iter=False, positive=False,
              check_input=True, **params):
    """
    Compute elastic net path with coordinate descent.

    The elastic net optimization function varies for mono and multi-outputs.

    For mono-output tasks it is::

        1 / (2 * n_samples) * ||y - Xw||^2_2
        + alpha * l1_ratio * ||w||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2

    For multi-output tasks it is::

        (1 / (2 * n_samples)) * ||Y - XW||^Fro_2
        + alpha * l1_ratio * ||W||_21
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2

    Where::

        ||W||_21 = \\sum_i \\sqrt{\\sum_j w_{ij}^2}

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <elastic_net>`.

    Parameters
    ----------
    X : {array-like}, shape (n_samples, n_features)
        Training data. Pass directly as Fortran-contiguous data to avoid
        unnecessary memory duplication. If ``y`` is mono-output then ``X``
        can be sparse.

    y : ndarray, shape (n_samples,) or (n_samples, n_outputs)
        Target values.

    l1_ratio : float, optional
        Number between 0 and 1 passed to elastic net (scaling between
        l1 and l2 penalties). ``l1_ratio=1`` corresponds to the Lasso.

    eps : float
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional
        Number of alphas along the regularization path.

    alphas : ndarray, optional
        List of alphas where to compute the models.
        If None alphas are set automatically.

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument.

    Xy : array-like, optional
        Xy = np.dot(X.T, y) that can be precomputed. It is useful
        only when the Gram matrix is precomputed.

    copy_X : bool, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    coef_init : array, shape (n_features, ) | None
        The initial values of the coefficients.

    verbose : bool or int
        Amount of verbosity.

    return_n_iter : bool
        Whether to return the number of iterations or not.

    positive : bool, default False
        If set to True, forces coefficients to be positive.
        (Only allowed when ``y.ndim == 1``).

    check_input : bool, default True
        Skip input validation checks, including the Gram matrix when provided
        assuming there are handled by the caller when check_input=False.

    **params : kwargs
        Keyword arguments passed to the coordinate descent solver.

    Returns
    -------
    alphas : array, shape (n_alphas,)
        The alphas along the path where models are computed.

    coefs : array, shape (n_features, n_alphas) or \
            (n_outputs, n_features, n_alphas)
        Coefficients along the path.

    dual_gaps : array, shape (n_alphas,)
        The dual gaps at the end of the optimization for each alpha.

    n_iters : array-like, shape (n_alphas,)
        The number of iterations taken by the coordinate descent optimizer to
        reach the specified tolerance for each alpha.
        (Is returned when ``return_n_iter`` is set to True).

    See Also
    --------
    MultiTaskElasticNet
    MultiTaskElasticNetCV
    ElasticNet
    ElasticNetCV

    Notes
    -----
    For an example, see
    :ref:`examples/linear_model/plot_lasso_coordinate_descent_path.py
    <sphx_glr_auto_examples_linear_model_plot_lasso_coordinate_descent_path.py>`.
    """
    # We expect X and y to be already Fortran ordered when bypassing
    # checks
    if check_input:
        X = check_array(X, 'csc', dtype=[np.float64, np.float32],
                        order='F', copy=copy_X)
        y = check_array(y, 'csc', dtype=X.dtype.type, order='F', copy=False,
                        ensure_2d=False)
        if Xy is not None:
            # Xy should be a 1d contiguous array or a 2D C ordered array
            Xy = check_array(Xy, dtype=X.dtype.type, order='C', copy=False,
                             ensure_2d=False)

    n_samples, n_features = X.shape

    multi_output = False
    if y.ndim != 1:
        multi_output = True
        _, n_outputs = y.shape

    if multi_output and positive:
        raise ValueError('positive=True is not allowed for multi-output'
                         ' (y.ndim != 1)')

    # MultiTaskElasticNet does not support sparse matrices
    if not multi_output and sparse.isspmatrix(X):
        if 'X_offset' in params:
            # As sparse matrices are not actually centered we need this
            # to be passed to the CD solver.
            X_sparse_scaling = params['X_offset'] / params['X_scale']
            X_sparse_scaling = np.asarray(X_sparse_scaling, dtype=X.dtype)
        else:
            X_sparse_scaling = np.zeros(n_features, dtype=X.dtype)

    # X should be normalized and fit already if function is called
    # from ElasticNet.fit
    if check_input:
        X, y, X_offset, y_offset, X_scale, precompute, Xy = \
            _pre_fit(X, y, Xy, precompute, normalize=False,
                     fit_intercept=False, copy=False, check_input=check_input)
    if alphas is None:
        # No need to normalize of fit_intercept: it has been done
        # above
        alphas = _alpha_grid(X, y, Xy=Xy, l1_ratio=l1_ratio,
                             fit_intercept=False, eps=eps, n_alphas=n_alphas,
                             normalize=False, copy_X=False)
    else:
        alphas = np.sort(alphas)[::-1]  # make sure alphas are properly ordered

    n_alphas = len(alphas)
    tol = params.get('tol', 1e-4)
    max_iter = params.get('max_iter', 1000)
    dual_gaps = np.empty(n_alphas)
    n_iters = []

    rng = check_random_state(params.get('random_state', None))
    selection = params.get('selection', 'cyclic')
    if selection not in ['random', 'cyclic']:
        raise ValueError("selection should be either random or cyclic.")
    random = (selection == 'random')

    if not multi_output:
        coefs = np.empty((n_features, n_alphas), dtype=X.dtype)
    else:
        coefs = np.empty((n_outputs, n_features, n_alphas),
                         dtype=X.dtype)

    if coef_init is None:
        coef_ = np.zeros(coefs.shape[:-1], dtype=X.dtype, order='F')
    else:
        coef_ = np.asfortranarray(coef_init, dtype=X.dtype)

    for i, alpha in enumerate(alphas):
        l1_reg = alpha * l1_ratio * n_samples
        l2_reg = alpha * (1.0 - l1_ratio) * n_samples
        if not multi_output and sparse.isspmatrix(X):
            model = cd_fast.sparse_enet_coordinate_descent(
                coef_, l1_reg, l2_reg, X.data, X.indices,
                X.indptr, y, X_sparse_scaling,
                max_iter, tol, rng, random, positive)
        elif multi_output:
            model = cd_fast.enet_coordinate_descent_multi_task(
                coef_, l1_reg, l2_reg, X, y, max_iter, tol, rng, random)
        elif isinstance(precompute, np.ndarray):
            # We expect precompute to be already Fortran ordered when bypassing
            # checks
            if check_input:
                precompute = check_array(precompute, dtype=X.dtype.type,
                                         order='C')
            model = cd_fast.enet_coordinate_descent_gram(
                coef_, l1_reg, l2_reg, precompute, Xy, y, max_iter,
                tol, rng, random, positive)
        elif precompute is False:
            model = cd_fast.enet_coordinate_descent(
                coef_, l1_reg, l2_reg, X, y, max_iter, tol, rng, random,
                positive)
        else:
            raise ValueError("Precompute should be one of True, False, "
                             "'auto' or array-like. Got %r" % precompute)
        coef_, dual_gap_, eps_, n_iter_ = model
        coefs[..., i] = coef_
        dual_gaps[i] = dual_gap_
        n_iters.append(n_iter_)

        if verbose:
            if verbose > 2:
                print(model)
            elif verbose > 1:
                print('Path: %03i out of %03i' % (i, n_alphas))
            else:
                sys.stderr.write('.')

    if return_n_iter:
        return alphas, coefs, dual_gaps, n_iters
    return alphas, coefs, dual_gaps


###############################################################################
# ElasticNet model


class ElasticNet(MultiOutputMixin, RegressorMixin, LinearModel):
    """Linear regression with combined L1 and L2 priors as regularizer.

    Minimizes the objective function::

            1 / (2 * n_samples) * ||y - Xw||^2_2
            + alpha * l1_ratio * ||w||_1
            + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2

    If you are interested in controlling the L1 and L2 penalty
    separately, keep in mind that this is equivalent to::

            a * L1 + b * L2

    where::

            alpha = a + b and l1_ratio = a / (a + b)

    The parameter l1_ratio corresponds to alpha in the glmnet R package while
    alpha corresponds to the lambda parameter in glmnet. Specifically, l1_ratio
    = 1 is the lasso penalty. Currently, l1_ratio <= 0.01 is not reliable,
    unless you supply your own sequence of alpha.

    Read more in the :ref:`User Guide <elastic_net>`.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the penalty terms. Defaults to 1.0.
        See the notes for the exact mathematical meaning of this
        parameter. ``alpha = 0`` is equivalent to an ordinary least square,
        solved by the :class:`LinearRegression` object. For numerical
        reasons, using ``alpha = 0`` with the ``Lasso`` object is not advised.
        Given this, you should use the :class:`LinearRegression` object.

    l1_ratio : float
        The ElasticNet mixing parameter, with ``0 <= l1_ratio <= 1``. For
        ``l1_ratio = 0`` the penalty is an L2 penalty. ``For l1_ratio = 1`` it
        is an L1 penalty.  For ``0 < l1_ratio < 1``, the penalty is a
        combination of L1 and L2.

    fit_intercept : bool
        Whether the intercept should be estimated or not. If ``False``, the
        data is assumed to be already centered.

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    precompute : True | False | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. The Gram matrix can also be passed as argument.
        For sparse input this option is always ``True`` to preserve sparsity.

    max_iter : int, optional
        The maximum number of iterations

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    warm_start : bool, optional
        When set to ``True``, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

    positive : bool, optional
        When set to ``True``, forces the coefficients to be positive.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    coef_ : array, shape (n_features,) | (n_targets, n_features)
        parameter vector (w in the cost function formula)

    sparse_coef_ : scipy.sparse matrix, shape (n_features, 1) | \
            (n_targets, n_features)
        ``sparse_coef_`` is a readonly property derived from ``coef_``

    intercept_ : float | array, shape (n_targets,)
        independent term in decision function.

    n_iter_ : array-like, shape (n_targets,)
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance.

    Examples
    --------
    >>> from sklearn.linear_model import ElasticNet
    >>> from sklearn.datasets import make_regression

    >>> X, y = make_regression(n_features=2, random_state=0)
    >>> regr = ElasticNet(random_state=0)
    >>> regr.fit(X, y)
    ElasticNet(random_state=0)
    >>> print(regr.coef_)
    [18.83816048 64.55968825]
    >>> print(regr.intercept_)
    1.451...
    >>> print(regr.predict([[0, 0]]))
    [1.451...]


    Notes
    -----
    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.

    See also
    --------
    ElasticNetCV : Elastic net model with best model selection by
        cross-validation.
    SGDRegressor: implements elastic net regression with incremental training.
    SGDClassifier: implements logistic regression with elastic net penalty
        (``SGDClassifier(loss="log", penalty="elasticnet")``).
    """
    path = staticmethod(enet_path)

    def __init__(self, alpha=1.0, l1_ratio=0.5, fit_intercept=True,
                 normalize=False, precompute=False, max_iter=1000,
                 copy_X=True, tol=1e-4, warm_start=False, positive=False,
                 random_state=None, selection='cyclic'):
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.positive = positive
        self.random_state = random_state
        self.selection = selection

    def fit(self, X, y, check_input=True):
        """Fit model with coordinate descent.

        Parameters
        ----------
        X : ndarray or scipy.sparse matrix, (n_samples, n_features)
            Data

        y : ndarray, shape (n_samples,) or (n_samples, n_targets)
            Target. Will be cast to X's dtype if necessary

        check_input : boolean, (default=True)
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.

        Notes
        -----

        Coordinate descent is an algorithm that considers each column of
        data at a time hence it will automatically convert the X input
        as a Fortran-contiguous numpy array if necessary.

        To avoid memory re-allocation it is advised to allocate the
        initial data in memory directly using that format.
        """

        if self.alpha == 0:
            warnings.warn("With alpha=0, this algorithm does not converge "
                          "well. You are advised to use the LinearRegression "
                          "estimator", stacklevel=2)

        if isinstance(self.precompute, str):
            raise ValueError('precompute should be one of True, False or'
                             ' array-like. Got %r' % self.precompute)

        # Remember if X is copied
        X_copied = False
        # We expect X and y to be float64 or float32 Fortran ordered arrays
        # when bypassing checks
        if check_input:
            X_copied = self.copy_X and self.fit_intercept
            X, y = check_X_y(X, y, accept_sparse='csc',
                             order='F', dtype=[np.float64, np.float32],
                             copy=X_copied, multi_output=True, y_numeric=True)
            y = check_array(y, order='F', copy=False, dtype=X.dtype.type,
                            ensure_2d=False)

        # Ensure copying happens only once, don't do it again if done above
        should_copy = self.copy_X and not X_copied
        X, y, X_offset, y_offset, X_scale, precompute, Xy = \
            _pre_fit(X, y, None, self.precompute, self.normalize,
                     self.fit_intercept, copy=should_copy,
                     check_input=check_input)
        if y.ndim == 1:
            y = y[:, np.newaxis]
        if Xy is not None and Xy.ndim == 1:
            Xy = Xy[:, np.newaxis]

        n_samples, n_features = X.shape
        n_targets = y.shape[1]

        if self.selection not in ['cyclic', 'random']:
            raise ValueError("selection should be either random or cyclic.")

        if not self.warm_start or not hasattr(self, "coef_"):
            coef_ = np.zeros((n_targets, n_features), dtype=X.dtype,
                             order='F')
        else:
            coef_ = self.coef_
            if coef_.ndim == 1:
                coef_ = coef_[np.newaxis, :]

        dual_gaps_ = np.zeros(n_targets, dtype=X.dtype)
        self.n_iter_ = []

        for k in range(n_targets):
            if Xy is not None:
                this_Xy = Xy[:, k]
            else:
                this_Xy = None
            _, this_coef, this_dual_gap, this_iter = \
                self.path(X, y[:, k],
                          l1_ratio=self.l1_ratio, eps=None,
                          n_alphas=None, alphas=[self.alpha],
                          precompute=precompute, Xy=this_Xy,
                          fit_intercept=False, normalize=False, copy_X=True,
                          verbose=False, tol=self.tol, positive=self.positive,
                          X_offset=X_offset, X_scale=X_scale,
                          return_n_iter=True, coef_init=coef_[k],
                          max_iter=self.max_iter,
                          random_state=self.random_state,
                          selection=self.selection,
                          check_input=False)
            coef_[k] = this_coef[:, 0]
            dual_gaps_[k] = this_dual_gap[0]
            self.n_iter_.append(this_iter[0])

        if n_targets == 1:
            self.n_iter_ = self.n_iter_[0]
            self.coef_ = coef_[0]
            self.dual_gap_ = dual_gaps_[0]
        else:
            self.coef_ = coef_
            self.dual_gap_ = dual_gaps_

        self._set_intercept(X_offset, y_offset, X_scale)

        # workaround since _set_intercept will cast self.coef_ into X.dtype
        self.coef_ = np.asarray(self.coef_, dtype=X.dtype)

        # return self for chaining fit and predict calls
        return self

    @property
    def sparse_coef_(self):
        """ sparse representation of the fitted ``coef_`` """
        return sparse.csr_matrix(self.coef_)

    def _decision_function(self, X):
        """Decision function of the linear model

        Parameters
        ----------
        X : numpy array or scipy.sparse matrix of shape (n_samples, n_features)

        Returns
        -------
        T : array, shape (n_samples,)
            The predicted decision function
        """
        check_is_fitted(self)
        if sparse.isspmatrix(X):
            return safe_sparse_dot(X, self.coef_.T,
                                   dense_output=True) + self.intercept_
        else:
            return super()._decision_function(X)


###############################################################################
# Lasso model

class Lasso(ElasticNet):
    """Linear Model trained with L1 prior as regularizer (aka the Lasso)

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

    Technically the Lasso model is optimizing the same objective function as
    the Elastic Net with ``l1_ratio=1.0`` (no L2 penalty).

    Read more in the :ref:`User Guide <lasso>`.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1 term. Defaults to 1.0.
        ``alpha = 0`` is equivalent to an ordinary least square, solved
        by the :class:`LinearRegression` object. For numerical
        reasons, using ``alpha = 0`` with the ``Lasso`` object is not advised.
        Given this, you should use the :class:`LinearRegression` object.

    fit_intercept : boolean, optional, default True
        Whether to calculate the intercept for this model. If set
        to False, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    precompute : True | False | array-like, default=False
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument. For sparse input
        this option is always ``True`` to preserve sparsity.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

    positive : bool, optional
        When set to ``True``, forces the coefficients to be positive.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    coef_ : array, shape (n_features,) | (n_targets, n_features)
        parameter vector (w in the cost function formula)

    sparse_coef_ : scipy.sparse matrix, shape (n_features, 1) | \
            (n_targets, n_features)
        ``sparse_coef_`` is a readonly property derived from ``coef_``

    intercept_ : float | array, shape (n_targets,)
        independent term in decision function.

    n_iter_ : int | array-like, shape (n_targets,)
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.Lasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso(alpha=0.1)
    >>> print(clf.coef_)
    [0.85 0.  ]
    >>> print(clf.intercept_)
    0.15...

    See also
    --------
    lars_path
    lasso_path
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.
    """
    path = staticmethod(enet_path)

    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 precompute=False, copy_X=True, max_iter=1000,
                 tol=1e-4, warm_start=False, positive=False,
                 random_state=None, selection='cyclic'):
        super().__init__(
            alpha=alpha, l1_ratio=1.0, fit_intercept=fit_intercept,
            normalize=normalize, precompute=precompute, copy_X=copy_X,
            max_iter=max_iter, tol=tol, warm_start=warm_start,
            positive=positive, random_state=random_state,
            selection=selection)


###############################################################################
# Functions for CV with paths functions

def _path_residuals(X, y, train, test, path, path_params, alphas=None,
                    l1_ratio=1, X_order=None, dtype=None):
    """Returns the MSE for the models computed by 'path'

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Training data.

    y : array-like, shape (n_samples,) or (n_samples, n_targets)
        Target values

    train : list of indices
        The indices of the train set

    test : list of indices
        The indices of the test set

    path : callable
        function returning a list of models on the path. See
        enet_path for an example of signature

    path_params : dictionary
        Parameters passed to the path function

    alphas : array-like, optional
        Array of float that is used for cross-validation. If not
        provided, computed using 'path'.

    l1_ratio : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). For ``l1_ratio = 0`` the penalty is an
        L2 penalty. For ``l1_ratio = 1`` it is an L1 penalty. For ``0
        < l1_ratio < 1``, the penalty is a combination of L1 and L2.

    X_order : {'F', 'C', or None}, optional
        The order of the arrays expected by the path function to
        avoid memory copies

    dtype : a numpy dtype or None
        The dtype of the arrays expected by the path function to
        avoid memory copies
    """
    X_train = X[train]
    y_train = y[train]
    X_test = X[test]
    y_test = y[test]
    fit_intercept = path_params['fit_intercept']
    normalize = path_params['normalize']

    if y.ndim == 1:
        precompute = path_params['precompute']
    else:
        # No Gram variant of multi-task exists right now.
        # Fall back to default enet_multitask
        precompute = False

    X_train, y_train, X_offset, y_offset, X_scale, precompute, Xy = \
        _pre_fit(X_train, y_train, None, precompute, normalize, fit_intercept,
                 copy=False)

    path_params = path_params.copy()
    path_params['Xy'] = Xy
    path_params['X_offset'] = X_offset
    path_params['X_scale'] = X_scale
    path_params['precompute'] = precompute
    path_params['copy_X'] = False
    path_params['alphas'] = alphas

    if 'l1_ratio' in path_params:
        path_params['l1_ratio'] = l1_ratio

    # Do the ordering and type casting here, as if it is done in the path,
    # X is copied and a reference is kept here
    X_train = check_array(X_train, 'csc', dtype=dtype, order=X_order)
    alphas, coefs, _ = path(X_train, y_train, **path_params)
    del X_train, y_train

    if y.ndim == 1:
        # Doing this so that it becomes coherent with multioutput.
        coefs = coefs[np.newaxis, :, :]
        y_offset = np.atleast_1d(y_offset)
        y_test = y_test[:, np.newaxis]

    if normalize:
        nonzeros = np.flatnonzero(X_scale)
        coefs[:, nonzeros] /= X_scale[nonzeros][:, np.newaxis]

    intercepts = y_offset[:, np.newaxis] - np.dot(X_offset, coefs)
    X_test_coefs = safe_sparse_dot(X_test, coefs)
    residues = X_test_coefs - y_test[:, :, np.newaxis]
    residues += intercepts
    this_mses = ((residues ** 2).mean(axis=0)).mean(axis=0)

    return this_mses


class LinearModelCV(MultiOutputMixin, LinearModel, metaclass=ABCMeta):
    """Base class for iterative model fitting along a regularization path"""

    @abstractmethod
    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
                 normalize=False, precompute='auto', max_iter=1000, tol=1e-4,
                 copy_X=True, cv=None, verbose=False, n_jobs=None,
                 positive=False, random_state=None, selection='cyclic'):
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.tol = tol
        self.copy_X = copy_X
        self.cv = cv
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.positive = positive
        self.random_state = random_state
        self.selection = selection

    def fit(self, X, y):
        """Fit linear model with coordinate descent

        Fit is on grid of alphas and best alpha estimated by cross-validation.

        Parameters
        ----------
        X : {array-like}, shape (n_samples, n_features)
            Training data. Pass directly as Fortran-contiguous data
            to avoid unnecessary memory duplication. If y is mono-output,
            X can be sparse.

        y : array-like, shape (n_samples,) or (n_samples, n_targets)
            Target values
        """
        y = check_array(y, copy=False, dtype=[np.float64, np.float32],
                        ensure_2d=False)
        if y.shape[0] == 0:
            raise ValueError("y has 0 samples: %r" % y)

        if hasattr(self, 'l1_ratio'):
            model_str = 'ElasticNet'
        else:
            model_str = 'Lasso'

        if isinstance(self, ElasticNetCV) or isinstance(self, LassoCV):
            if model_str == 'ElasticNet':
                model = ElasticNet()
            else:
                model = Lasso()
            if y.ndim > 1 and y.shape[1] > 1:
                raise ValueError("For multi-task outputs, use "
                                 "MultiTask%sCV" % (model_str))
            y = column_or_1d(y, warn=True)
        else:
            if sparse.isspmatrix(X):
                raise TypeError("X should be dense but a sparse matrix was"
                                "passed")
            elif y.ndim == 1:
                raise ValueError("For mono-task outputs, use "
                                 "%sCV" % (model_str))
            if model_str == 'ElasticNet':
                model = MultiTaskElasticNet()
            else:
                model = MultiTaskLasso()

        if self.selection not in ["random", "cyclic"]:
            raise ValueError("selection should be either random or cyclic.")

        # This makes sure that there is no duplication in memory.
        # Dealing right with copy_X is important in the following:
        # Multiple functions touch X and subsamples of X and can induce a
        # lot of duplication of memory
        copy_X = self.copy_X and self.fit_intercept

        if isinstance(X, np.ndarray) or sparse.isspmatrix(X):
            # Keep a reference to X
            reference_to_old_X = X
            # Let us not impose fortran ordering so far: it is
            # not useful for the cross-validation loop and will be done
            # by the model fitting itself
            X = check_array(X, 'csc', dtype=[np.float64, np.float32],
                            copy=False)
            if sparse.isspmatrix(X):
                if (hasattr(reference_to_old_X, "data") and
                   not np.may_share_memory(reference_to_old_X.data, X.data)):
                    # X is a sparse matrix and has been copied
                    copy_X = False
            elif not np.may_share_memory(reference_to_old_X, X):
                # X has been copied
                copy_X = False
            del reference_to_old_X
        else:
            X = check_array(X, 'csc', dtype=[np.float64, np.float32],
                            order='F', copy=copy_X)
            copy_X = False

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have inconsistent dimensions (%d != %d)"
                             % (X.shape[0], y.shape[0]))

        # All LinearModelCV parameters except 'cv' are acceptable
        path_params = self.get_params()
        if 'l1_ratio' in path_params:
            l1_ratios = np.atleast_1d(path_params['l1_ratio'])
            # For the first path, we need to set l1_ratio
            path_params['l1_ratio'] = l1_ratios[0]
        else:
            l1_ratios = [1, ]
        path_params.pop('cv', None)
        path_params.pop('n_jobs', None)

        alphas = self.alphas
        n_l1_ratio = len(l1_ratios)
        if alphas is None:
            alphas = [_alpha_grid(X, y, l1_ratio=l1_ratio,
                                  fit_intercept=self.fit_intercept,
                                  eps=self.eps, n_alphas=self.n_alphas,
                                  normalize=self.normalize, copy_X=self.copy_X)
                      for l1_ratio in l1_ratios]
        else:
            # Making sure alphas is properly ordered.
            alphas = np.tile(np.sort(alphas)[::-1], (n_l1_ratio, 1))
        # We want n_alphas to be the number of alphas used for each l1_ratio.
        n_alphas = len(alphas[0])
        path_params.update({'n_alphas': n_alphas})

        path_params['copy_X'] = copy_X
        # We are not computing in parallel, we can modify X
        # inplace in the folds
        if effective_n_jobs(self.n_jobs) > 1:
            path_params['copy_X'] = False

        # init cross-validation generator
        cv = check_cv(self.cv)

        # Compute path for all folds and compute MSE to get the best alpha
        folds = list(cv.split(X, y))
        best_mse = np.inf

        # We do a double for loop folded in one, in order to be able to
        # iterate in parallel on l1_ratio and folds
        jobs = (delayed(_path_residuals)(X, y, train, test, self.path,
                                         path_params, alphas=this_alphas,
                                         l1_ratio=this_l1_ratio, X_order='F',
                                         dtype=X.dtype.type)
                for this_l1_ratio, this_alphas in zip(l1_ratios, alphas)
                for train, test in folds)
        mse_paths = Parallel(n_jobs=self.n_jobs, verbose=self.verbose,
                             **_joblib_parallel_args(prefer="threads"))(jobs)
        mse_paths = np.reshape(mse_paths, (n_l1_ratio, len(folds), -1))
        mean_mse = np.mean(mse_paths, axis=1)
        self.mse_path_ = np.squeeze(np.rollaxis(mse_paths, 2, 1))
        for l1_ratio, l1_alphas, mse_alphas in zip(l1_ratios, alphas,
                                                   mean_mse):
            i_best_alpha = np.argmin(mse_alphas)
            this_best_mse = mse_alphas[i_best_alpha]
            if this_best_mse < best_mse:
                best_alpha = l1_alphas[i_best_alpha]
                best_l1_ratio = l1_ratio
                best_mse = this_best_mse

        self.l1_ratio_ = best_l1_ratio
        self.alpha_ = best_alpha
        if self.alphas is None:
            self.alphas_ = np.asarray(alphas)
            if n_l1_ratio == 1:
                self.alphas_ = self.alphas_[0]
        # Remove duplicate alphas in case alphas is provided.
        else:
            self.alphas_ = np.asarray(alphas[0])

        # Refit the model with the parameters selected
        common_params = {name: value
                         for name, value in self.get_params().items()
                         if name in model.get_params()}
        model.set_params(**common_params)
        model.alpha = best_alpha
        model.l1_ratio = best_l1_ratio
        model.copy_X = copy_X
        precompute = getattr(self, "precompute", None)
        if isinstance(precompute, str) and precompute == "auto":
            model.precompute = False
        model.fit(X, y)
        if not hasattr(self, 'l1_ratio'):
            del self.l1_ratio_
        self.coef_ = model.coef_
        self.intercept_ = model.intercept_
        self.dual_gap_ = model.dual_gap_
        self.n_iter_ = model.n_iter_
        return self


class LassoCV(RegressorMixin, LinearModelCV):
    """Lasso linear model with iterative fitting along a regularization path.

    See glossary entry for :term:`cross-validation estimator`.

    The best model is selected by cross-validation.

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

    Read more in the :ref:`User Guide <lasso>`.

    Parameters
    ----------
    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If ``None`` alphas are set automatically

    fit_intercept : boolean, default True
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    verbose : bool or integer
        Amount of verbosity.

    n_jobs : int or None, optional (default=None)
        Number of CPUs to use during the cross validation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    positive : bool, optional
        If positive, restrict regression coefficients to be positive

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    alpha_ : float
        The amount of penalization chosen by cross validation

    coef_ : array, shape (n_features,) | (n_targets, n_features)
        parameter vector (w in the cost function formula)

    intercept_ : float | array, shape (n_targets,)
        independent term in decision function.

    mse_path_ : array, shape (n_alphas, n_folds)
        mean square error for the test set on each fold, varying alpha

    alphas_ : numpy array, shape (n_alphas,)
        The grid of alphas used for fitting

    dual_gap_ : ndarray, shape ()
        The dual gap at the end of the optimization for the optimal alpha
        (``alpha_``).

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance for the optimal alpha.

    Examples
    --------
    >>> from sklearn.linear_model import LassoCV
    >>> from sklearn.datasets import make_regression
    >>> X, y = make_regression(noise=4, random_state=0)
    >>> reg = LassoCV(cv=5, random_state=0).fit(X, y)
    >>> reg.score(X, y)
    0.9993...
    >>> reg.predict(X[:1,])
    array([-78.4951...])

    Notes
    -----
    For an example, see
    :ref:`examples/linear_model/plot_lasso_model_selection.py
    <sphx_glr_auto_examples_linear_model_plot_lasso_model_selection.py>`.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.

    See also
    --------
    lars_path
    lasso_path
    LassoLars
    Lasso
    LassoLarsCV
    """
    path = staticmethod(lasso_path)

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
                 normalize=False, precompute='auto', max_iter=1000, tol=1e-4,
                 copy_X=True, cv=None, verbose=False, n_jobs=None,
                 positive=False, random_state=None, selection='cyclic'):
        super().__init__(
            eps=eps, n_alphas=n_alphas, alphas=alphas,
            fit_intercept=fit_intercept, normalize=normalize,
            precompute=precompute, max_iter=max_iter, tol=tol, copy_X=copy_X,
            cv=cv, verbose=verbose, n_jobs=n_jobs, positive=positive,
            random_state=random_state, selection=selection)

    def _more_tags(self):
        return {'multioutput': False}

class ElasticNetCV(RegressorMixin, LinearModelCV):
    """Elastic Net model with iterative fitting along a regularization path.

    See glossary entry for :term:`cross-validation estimator`.

    Read more in the :ref:`User Guide <elastic_net>`.

    Parameters
    ----------
    l1_ratio : float or array of floats, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). For ``l1_ratio = 0``
        the penalty is an L2 penalty. For ``l1_ratio = 1`` it is an L1 penalty.
        For ``0 < l1_ratio < 1``, the penalty is a combination of L1 and L2
        This parameter can be a list, in which case the different
        values are tested by cross-validation and the one giving the best
        prediction score is used. Note that a good choice of list of
        values for l1_ratio is often to put more values close to 1
        (i.e. Lasso) and less close to 0 (i.e. Ridge), as in ``[.1, .5, .7,
        .9, .95, .99, 1]``

    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional
        Number of alphas along the regularization path, used for each l1_ratio.

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : bool or integer
        Amount of verbosity.

    n_jobs : int or None, optional (default=None)
        Number of CPUs to use during the cross validation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    positive : bool, optional
        When set to ``True``, forces the coefficients to be positive.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    alpha_ : float
        The amount of penalization chosen by cross validation

    l1_ratio_ : float
        The compromise between l1 and l2 penalization chosen by
        cross validation

    coef_ : array, shape (n_features,) | (n_targets, n_features)
        Parameter vector (w in the cost function formula),

    intercept_ : float | array, shape (n_targets, n_features)
        Independent term in the decision function.

    mse_path_ : array, shape (n_l1_ratio, n_alpha, n_folds)
        Mean square error for the test set on each fold, varying l1_ratio and
        alpha.

    alphas_ : numpy array, shape (n_alphas,) or (n_l1_ratio, n_alphas)
        The grid of alphas used for fitting, for each l1_ratio.

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance for the optimal alpha.

    Examples
    --------
    >>> from sklearn.linear_model import ElasticNetCV
    >>> from sklearn.datasets import make_regression

    >>> X, y = make_regression(n_features=2, random_state=0)
    >>> regr = ElasticNetCV(cv=5, random_state=0)
    >>> regr.fit(X, y)
    ElasticNetCV(cv=5, random_state=0)
    >>> print(regr.alpha_)
    0.199...
    >>> print(regr.intercept_)
    0.398...
    >>> print(regr.predict([[0, 0]]))
    [0.398...]


    Notes
    -----
    For an example, see
    :ref:`examples/linear_model/plot_lasso_model_selection.py
    <sphx_glr_auto_examples_linear_model_plot_lasso_model_selection.py>`.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.

    The parameter l1_ratio corresponds to alpha in the glmnet R package
    while alpha corresponds to the lambda parameter in glmnet.
    More specifically, the optimization objective is::

        1 / (2 * n_samples) * ||y - Xw||^2_2
        + alpha * l1_ratio * ||w||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2

    If you are interested in controlling the L1 and L2 penalty
    separately, keep in mind that this is equivalent to::

        a * L1 + b * L2

    for::

        alpha = a + b and l1_ratio = a / (a + b).

    See also
    --------
    enet_path
    ElasticNet

    """
    path = staticmethod(enet_path)

    def __init__(self, l1_ratio=0.5, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True, normalize=False, precompute='auto',
                 max_iter=1000, tol=1e-4, cv=None, copy_X=True,
                 verbose=0, n_jobs=None, positive=False, random_state=None,
                 selection='cyclic'):
        self.l1_ratio = l1_ratio
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.tol = tol
        self.cv = cv
        self.copy_X = copy_X
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.positive = positive
        self.random_state = random_state
        self.selection = selection

    def _more_tags(self):
        return {'multioutput': False}

###############################################################################
# Multi Task ElasticNet and Lasso models (with joint feature selection)


class MultiTaskElasticNet(Lasso):
    """Multi-task ElasticNet model trained with L1/L2 mixed-norm as regularizer

    The optimization objective for MultiTaskElasticNet is::

        (1 / (2 * n_samples)) * ||Y - XW||_Fro^2
        + alpha * l1_ratio * ||W||_21
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2

    Where::

        ||W||_21 = sum_i sqrt(sum_j w_ij ^ 2)

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <multi_task_elastic_net>`.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1/L2 term. Defaults to 1.0

    l1_ratio : float
        The ElasticNet mixing parameter, with 0 < l1_ratio <= 1.
        For l1_ratio = 1 the penalty is an L1/L2 penalty. For l1_ratio = 0 it
        is an L2 penalty.
        For ``0 < l1_ratio < 1``, the penalty is a combination of L1/L2 and L2.

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    warm_start : bool, optional
        When set to ``True``, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    intercept_ : array, shape (n_tasks,)
        Independent term in decision function.

    coef_ : array, shape (n_tasks, n_features)
        Parameter vector (W in the cost function formula). If a 1D y is
        passed in at fit (non multi-task usage), ``coef_`` is then a 1D array.
        Note that ``coef_`` stores the transpose of ``W``, ``W.T``.

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.MultiTaskElasticNet(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [[0, 0], [1, 1], [2, 2]])
    MultiTaskElasticNet(alpha=0.1)
    >>> print(clf.coef_)
    [[0.45663524 0.45612256]
     [0.45663524 0.45612256]]
    >>> print(clf.intercept_)
    [0.0872422 0.0872422]

    See also
    --------
    MultiTaskElasticNet : Multi-task L1/L2 ElasticNet with built-in
        cross-validation.
    ElasticNet
    MultiTaskLasso

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.
    """
    def __init__(self, alpha=1.0, l1_ratio=0.5, fit_intercept=True,
                 normalize=False, copy_X=True, max_iter=1000, tol=1e-4,
                 warm_start=False, random_state=None, selection='cyclic'):
        self.l1_ratio = l1_ratio
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.random_state = random_state
        self.selection = selection

    def fit(self, X, y):
        """Fit MultiTaskElasticNet model with coordinate descent

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            Data
        y : ndarray, shape (n_samples, n_tasks)
            Target. Will be cast to X's dtype if necessary

        Notes
        -----

        Coordinate descent is an algorithm that considers each column of
        data at a time hence it will automatically convert the X input
        as a Fortran-contiguous numpy array if necessary.

        To avoid memory re-allocation it is advised to allocate the
        initial data in memory directly using that format.
        """
        X = check_array(X, dtype=[np.float64, np.float32], order='F',
                        copy=self.copy_X and self.fit_intercept)
        y = check_array(y, dtype=X.dtype.type, ensure_2d=False)

        if hasattr(self, 'l1_ratio'):
            model_str = 'ElasticNet'
        else:
            model_str = 'Lasso'
        if y.ndim == 1:
            raise ValueError("For mono-task outputs, use %s" % model_str)

        n_samples, n_features = X.shape
        _, n_tasks = y.shape

        if n_samples != y.shape[0]:
            raise ValueError("X and y have inconsistent dimensions (%d != %d)"
                             % (n_samples, y.shape[0]))

        X, y, X_offset, y_offset, X_scale = _preprocess_data(
            X, y, self.fit_intercept, self.normalize, copy=False)

        if not self.warm_start or not hasattr(self, "coef_"):
            self.coef_ = np.zeros((n_tasks, n_features), dtype=X.dtype.type,
                                  order='F')

        l1_reg = self.alpha * self.l1_ratio * n_samples
        l2_reg = self.alpha * (1.0 - self.l1_ratio) * n_samples

        self.coef_ = np.asfortranarray(self.coef_)  # coef contiguous in memory

        if self.selection not in ['random', 'cyclic']:
            raise ValueError("selection should be either random or cyclic.")
        random = (self.selection == 'random')

        self.coef_, self.dual_gap_, self.eps_, self.n_iter_ = \
            cd_fast.enet_coordinate_descent_multi_task(
                self.coef_, l1_reg, l2_reg, X, y, self.max_iter, self.tol,
                check_random_state(self.random_state), random)

        self._set_intercept(X_offset, y_offset, X_scale)

        # return self for chaining fit and predict calls
        return self

    def _more_tags(self):
        return {'multioutput_only': True}


class MultiTaskLasso(MultiTaskElasticNet):
    """Multi-task Lasso model trained with L1/L2 mixed-norm as regularizer.

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||Y - XW||^2_Fro + alpha * ||W||_21

    Where::

        ||W||_21 = \\sum_i \\sqrt{\\sum_j w_{ij}^2}

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <multi_task_lasso>`.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1/L2 term. Defaults to 1.0

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    warm_start : bool, optional
        When set to ``True``, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4

    Attributes
    ----------
    coef_ : array, shape (n_tasks, n_features)
        Parameter vector (W in the cost function formula).
        Note that ``coef_`` stores the transpose of ``W``, ``W.T``.

    intercept_ : array, shape (n_tasks,)
        independent term in decision function.

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.MultiTaskLasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [[0, 0], [1, 1], [2, 2]])
    MultiTaskLasso(alpha=0.1)
    >>> print(clf.coef_)
    [[0.89393398 0.        ]
     [0.89393398 0.        ]]
    >>> print(clf.intercept_)
    [0.10606602 0.10606602]

    See also
    --------
    MultiTaskLasso : Multi-task L1/L2 Lasso with built-in cross-validation
    Lasso
    MultiTaskElasticNet

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.
    """
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=1000, tol=1e-4, warm_start=False,
                 random_state=None, selection='cyclic'):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.l1_ratio = 1.0
        self.random_state = random_state
        self.selection = selection


class MultiTaskElasticNetCV(RegressorMixin, LinearModelCV):
    """Multi-task L1/L2 ElasticNet with built-in cross-validation.

    See glossary entry for :term:`cross-validation estimator`.

    The optimization objective for MultiTaskElasticNet is::

        (1 / (2 * n_samples)) * ||Y - XW||^Fro_2
        + alpha * l1_ratio * ||W||_21
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2

    Where::

        ||W||_21 = \\sum_i \\sqrt{\\sum_j w_{ij}^2}

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <multi_task_elastic_net>`.

    .. versionadded:: 0.15

    Parameters
    ----------
    l1_ratio : float or array of floats
        The ElasticNet mixing parameter, with 0 < l1_ratio <= 1.
        For l1_ratio = 1 the penalty is an L1/L2 penalty. For l1_ratio = 0 it
        is an L2 penalty.
        For ``0 < l1_ratio < 1``, the penalty is a combination of L1/L2 and L2.
        This parameter can be a list, in which case the different
        values are tested by cross-validation and the one giving the best
        prediction score is used. Note that a good choice of list of
        values for l1_ratio is often to put more values close to 1
        (i.e. Lasso) and less close to 0 (i.e. Ridge), as in ``[.1, .5, .7,
        .9, .95, .99, 1]``

    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : array-like, optional
        List of alphas where to compute the models.
        If not provided, set automatically.

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : bool or integer
        Amount of verbosity.

    n_jobs : int or None, optional (default=None)
        Number of CPUs to use during the cross validation. Note that this is
        used only if multiple values for l1_ratio are given.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    intercept_ : array, shape (n_tasks,)
        Independent term in decision function.

    coef_ : array, shape (n_tasks, n_features)
        Parameter vector (W in the cost function formula).
        Note that ``coef_`` stores the transpose of ``W``, ``W.T``.

    alpha_ : float
        The amount of penalization chosen by cross validation

    mse_path_ : array, shape (n_alphas, n_folds) or \
                (n_l1_ratio, n_alphas, n_folds)
        mean square error for the test set on each fold, varying alpha

    alphas_ : numpy array, shape (n_alphas,) or (n_l1_ratio, n_alphas)
        The grid of alphas used for fitting, for each l1_ratio

    l1_ratio_ : float
        best l1_ratio obtained by cross-validation.

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance for the optimal alpha.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.MultiTaskElasticNetCV(cv=3)
    >>> clf.fit([[0,0], [1, 1], [2, 2]],
    ...         [[0, 0], [1, 1], [2, 2]])
    MultiTaskElasticNetCV(cv=3)
    >>> print(clf.coef_)
    [[0.52875032 0.46958558]
     [0.52875032 0.46958558]]
    >>> print(clf.intercept_)
    [0.00166409 0.00166409]

    See also
    --------
    MultiTaskElasticNet
    ElasticNetCV
    MultiTaskLassoCV

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.
    """
    path = staticmethod(enet_path)

    def __init__(self, l1_ratio=0.5, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True, normalize=False,
                 max_iter=1000, tol=1e-4, cv=None, copy_X=True,
                 verbose=0, n_jobs=None, random_state=None,
                 selection='cyclic'):
        self.l1_ratio = l1_ratio
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.max_iter = max_iter
        self.tol = tol
        self.cv = cv
        self.copy_X = copy_X
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.selection = selection

    def _more_tags(self):
        return {'multioutput_only': True}


class MultiTaskLassoCV(RegressorMixin, LinearModelCV):
    """Multi-task Lasso model trained with L1/L2 mixed-norm as regularizer.

    See glossary entry for :term:`cross-validation estimator`.

    The optimization objective for MultiTaskLasso is::

        (1 / (2 * n_samples)) * ||Y - XW||^Fro_2 + alpha * ||W||_21

    Where::

        ||W||_21 = \\sum_i \\sqrt{\\sum_j w_{ij}^2}

    i.e. the sum of norm of each row.

    Read more in the :ref:`User Guide <multi_task_lasso>`.

    .. versionadded:: 0.15

    Parameters
    ----------
    eps : float, optional
        Length of the path. ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : array-like, optional
        List of alphas where to compute the models.
        If not provided, set automatically.

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    max_iter : int, optional
        The maximum number of iterations.

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than ``tol``, the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than ``tol``.

    copy_X : boolean, optional, default True
        If ``True``, X will be copied; else, it may be overwritten.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    verbose : bool or integer
        Amount of verbosity.

    n_jobs : int or None, optional (default=None)
        Number of CPUs to use during the cross validation. Note that this is
        used only if multiple values for l1_ratio are given.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'

    selection : str, default 'cyclic'
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially by default. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    Attributes
    ----------
    intercept_ : array, shape (n_tasks,)
        Independent term in decision function.

    coef_ : array, shape (n_tasks, n_features)
        Parameter vector (W in the cost function formula).
        Note that ``coef_`` stores the transpose of ``W``, ``W.T``.

    alpha_ : float
        The amount of penalization chosen by cross validation

    mse_path_ : array, shape (n_alphas, n_folds)
        mean square error for the test set on each fold, varying alpha

    alphas_ : numpy array, shape (n_alphas,)
        The grid of alphas used for fitting.

    n_iter_ : int
        number of iterations run by the coordinate descent solver to reach
        the specified tolerance for the optimal alpha.

    Examples
    --------
    >>> from sklearn.linear_model import MultiTaskLassoCV
    >>> from sklearn.datasets import make_regression
    >>> from sklearn.metrics import r2_score
    >>> X, y = make_regression(n_targets=2, noise=4, random_state=0)
    >>> reg = MultiTaskLassoCV(cv=5, random_state=0).fit(X, y)
    >>> r2_score(y, reg.predict(X))
    0.9994...
    >>> reg.alpha_
    0.5713...
    >>> reg.predict(X[:1,])
    array([[153.7971...,  94.9015...]])

    See also
    --------
    MultiTaskElasticNet
    ElasticNetCV
    MultiTaskElasticNetCV

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a Fortran-contiguous numpy array.
    """
    path = staticmethod(lasso_path)

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
                 normalize=False, max_iter=1000, tol=1e-4, copy_X=True,
                 cv=None, verbose=False, n_jobs=None, random_state=None,
                 selection='cyclic'):
        super().__init__(
            eps=eps, n_alphas=n_alphas, alphas=alphas,
            fit_intercept=fit_intercept, normalize=normalize,
            max_iter=max_iter, tol=tol, copy_X=copy_X,
            cv=cv, verbose=verbose, n_jobs=n_jobs, random_state=random_state,
            selection=selection)

    def _more_tags(self):
        return {'multioutput_only': True}
