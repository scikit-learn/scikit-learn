# coding: utf-8
"""
Large Margin Nearest Neighbor Classification
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
#
# License: BSD 3 clause

from __future__ import print_function
from warnings import warn

import numpy as np
import time
import sys
from scipy.optimize import fmin_l_bfgs_b
from scipy.sparse import csr_matrix, csc_matrix, spdiags

from ..base import BaseEstimator, TransformerMixin
from ..neighbors import NearestNeighbors
from ..decomposition import PCA
from ..utils import gen_batches
from ..utils.extmath import row_norms, safe_sparse_dot
from ..utils.random import check_random_state
from ..utils.multiclass import check_classification_targets
from ..utils.validation import check_is_fitted, check_array, check_X_y
from ..exceptions import DataDimensionalityWarning


class LargeMarginNearestNeighbor(BaseEstimator, TransformerMixin):
    """Distance Metric Learning for Large Margin Classification

    Parameters
    ----------
    init_transformation : array, shape (n_features_out, n_features_in),
    optional (default=None)
        An initial linear transformation. If None (default), the initial
        transformation is set to the identity, except if ``warm_start`` or
        ``init_pca`` is True.

    init_pca : bool, optional (default=True)
        Whether to use PCA to initialize the linear transformation.
        If False, the identity will be used, except if ``warm_start`` is True.

    warm_start : bool, optional (default=False)
        If True and :meth:`fit` has been called before, the solution of the
        previous call to :meth:`fit` is used as the initial linear
        transformation.

    n_features_out : int, optional (default=None)
        Preferred dimensionality of the inputs after the transformation.
        If None it is inferred from ``init_pca`` and ``init_transformation``.

    n_neighbors : int, optional (default=3)
        Number of neighbors to use as target neighbors for each sample.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the target neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    max_constraints : int, optional (default=500000)
        Maximum number of constraints to enforce per iteration.

    use_sparse : bool, optional (default=True)
        Whether to use a sparse matrix (default) or a dense matrix for the
        impostor-pairs storage. Using a sparse matrix, the distance to
        impostors is computed twice, but it is somewhat faster for larger
        data sets than using a dense matrix. With a dense matrix, the unique
        impostor pairs have to be identified explicitly.

    max_iter : int, optional (default=50)
        Maximum number of iterations in the optimization.

    tol : float, optional (default=1e-5)
        Convergence tolerance for the optimization.

    max_corrections : int, optional (default=100)
        The maximum number of variable metric corrections
        used to define the limited memory matrix. (The limited memory BFGS
        method does not store the full hessian but uses this many terms in an
        approximation to it.)

    callback : callable, optional (default=None)
        If not None, this function is called after every iteration of the
        optimizer taking as arguments the current solution and the number of
        iterations.


    verbose : int, optional (default=0)
        If 0, no progress messages will be printed.
        If 1, progress messages will be printed to stdout.
        If >1, progress messages will be printed and the ``iprint``
        parameter of :meth:`fmin_l_bfgs_b` of `scipy.optimize` will be set to
        verbose - 2.

    random_state : int or numpy.RandomState or None, optional (default=None)
        A pseudo random number generator used to sample from the constraints.

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Doesn't affect :meth:`fit` method.

    Attributes
    ----------
    transformation_ : array, shape (n_features_out, n_features_in).
        The linear transformation learned during fitting.

    n_neighbors_ : int
        The provided n_neighbors is decreased if it is greater than or equal
        to  min(number of elements in each class).

    n_features_out_ : int
        The dimensionality of a sample after applying to it the
        linear transformation.

    classes_non_singleton_ : array-like, shape (n_classes_non_singleton,)
        The appearing classes that have more than one sample.

    n_funcalls_ : int
        The number of times the optimizer computes the loss and the gradient.

    n_iter_ : int
        The number of iterations of the optimizer. Falls back to
        `n_funcalls` if the version of :meth:`fmin_l_bfgs_b` of
        `scipy.optimize` (< 0.12.0) does not store the number of iterations.

    details_ : dict
        A dictionary of information created by the L-BFGS optimizer during
        fitting.


    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import LargeMarginNearestNeighbor
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> lmnn = LargeMarginNearestNeighbor(n_neighbors=1)
    >>> lmnn.fit(X, y) # doctest: +ELLIPSIS
    LargeMarginNearestNeighbor(...)
    >>> print(lmnn.transform(X))
    [[ 0.        ]
     [ 0.52704628]
     [ 1.05409255]
     [ 1.58113883]]
    >>> knn = KNeighborsClassifier(n_neighbors=lmnn.n_neighbors_)
    >>> knn.fit(lmnn.transform(X), y) # doctest: +ELLIPSIS
    KNeighborsClassifier(...)
    >>> print(knn.score(lmnn.transform(X), y))
    1.0


    Notes
    -----
    Large margin nearest neighbor (LMNN) is a machine learning
    algorithm for metric learning. It learns a (pseudo-)metric in a
    supervised fashion to improve the classification accuracy of the k-nearest
    neighbor rule.
    The main intuition behind LMNN is to learn a pseudometric under which all
    data instances in the training set are surrounded by at least k instances
    that share the same class label. If this is achieved, the leave-one-out
    error is minimized.
    This implementation follows closely Kilian Weinberger's MATLAB code found
    at <https://bitbucket.org/mlcircus/lmnn> which solves the unconstrained
    problem, finding a linear transformation with L-BFGS instead of solving the
    constrained problem that finds the globally optimal metric. Different from
    the paper, the problem solved by this implementation is with the squared
    hinge loss (to make the problem differentiable).


    References
    ----------
    .. [1] Weinberger, Kilian Q., and Lawrence K. Saul. "Distance Metric
    Learning for Large Margin Nearest Neighbor Classification."
    Journal of Machine Learning Research, Vol. 10, Feb. 2009, pp. 207-244.
    (http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf)

    .. [2] Wikipedia entry on Large Margin Nearest Neighbor
    (https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor)

    """

    def __init__(self, init_transformation=None, init_pca=True,
                 warm_start=False, n_features_out=None, n_neighbors=3,
                 algorithm='auto', max_constraints=500000, use_sparse=True,
                 max_iter=50, tol=1e-5, max_corrections=100, callback=None,
                 verbose=0, random_state=None, n_jobs=1):

        # Parameters
        self.init_transformation = init_transformation
        self.init_pca = init_pca
        self.warm_start = warm_start
        self.n_features_out = n_features_out
        self.n_neighbors = n_neighbors
        self.algorithm = algorithm
        self.max_constraints = max_constraints
        self.use_sparse = use_sparse
        self.max_iter = max_iter
        self.tol = tol
        self.max_corrections = max_corrections
        self.callback = callback
        self.verbose = verbose
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Find a linear transformation by optimization of the unconstrained
        problem, such that the k-nearest neighbor classification accuracy
        improves.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features_in)
            The training samples.

        y : array-like, shape (n_samples,)
            The corresponding training labels.

        Returns
        -------
        self : returns a trained LargeMarginNearestNeighbor model.
        """

        # Check training data
        X, y = check_X_y(X, y, ensure_min_samples=2)
        check_classification_targets(y)

        # Check inputs consistency
        X_valid, y_valid = self._validate_params(X, y)

        self.random_state_ = check_random_state(self.random_state)

        # Initialize transformer
        transformation, self.n_features_out_ = self._init_transformer(X_valid)

        # Find the target neighbors
        targets = self._select_target_neighbors(X_valid, y_valid)

        # Compute gradient component of target neighbors
        grad_static = self._compute_grad_static(X_valid, targets, self.verbose)

        # Initialize number of optimizer iterations and objective funcalls
        self.n_iter_ = 0
        self.n_funcalls_ = 0
        iprint = self.verbose - 2 if self.verbose > 1 else -1

        # Create parameters dict for optimizer
        optimizer_params = {'func': self._loss_grad,
                            'x0': transformation,
                            'args': (X_valid, y_valid, targets, grad_static),
                            'm': self.max_corrections,
                            'pgtol': self.tol,
                            'iprint': iprint,
                            'maxiter': self.max_iter,
                            'callback': self._lbfgs_callback
                            }

        if self.verbose:
            print('\n{:>10} {:>10} {:>15} {:>10}'.
                  format('Iteration', 'Func.Call', 'Func.Value', 'Time(s)'))
            print('-' * 48)
            print('{:>10}'.format(self.n_iter_ + 1))

        # Call optimizer
        transformation, loss, info = fmin_l_bfgs_b(**optimizer_params)

        # Reshape result from optimizer
        d = transformation.size // self.n_features_out_
        self.transformation_ = transformation.reshape(self.n_features_out_, d)

        # Store information dict from the optimizer
        self.details_ = info
        self.details_['loss'] = loss
        self.n_iter_ = info['nit']

        if self.verbose:
            termination_reason = info['warnflag']
            n_funcalls = info['funcalls']
            if termination_reason == 0:
                print('Converged after {} function calls.'.format(n_funcalls))
            elif termination_reason == 1:
                print('Too many function evaluations ({}).'.format(n_funcalls))
            elif termination_reason == 2:
                print('Optimization stopped: {}'.format(info['task']))

        return self

    def transform(self, X, check_input=True):
        """Applies the learned transformation to the given data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features_in)
            Data samples.

        check_input: bool, optional (default=True)
            Whether to validate ``X``.

        Returns
        -------
        Lx: array-like, shape (n_samples, n_features_out)
            The data samples transformed.

        Raises
        ------
        NotFittedError
            If :meth:`fit` has not been called before.
        """

        check_is_fitted(self, ['transformation_'])

        if check_input:
            X = check_array(X)

        return X.dot(self.transformation_.T)

    def _validate_params(self, X, y):
        """Validate input parameters as soon as :meth:`fit` is called.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features_in)
            The training samples.

        y : array-like, shape (n_samples,)
            The corresponding training labels.

        Returns
        -------
        X : array, shape (n_samples, n_features_in)
            The validated training samples.

        y : array, shape (n_samples,)
            The validated corresponding training labels.

        Raises
        -------
        TypeError
            If a parameter's type does not match the desired type.

        ValueError
            If a parameter's value violates its legal value range or if the
            combination of two or more given parameters is incompatible.
        """

        # Find the appearing classes and the class index for each sample
        classes, y_inverse = np.unique(y, return_inverse=True)
        classes = np.arange(len(classes))

        # Ignore classes that have less than 2 samples
        class_sizes = np.bincount(y_inverse)
        is_class_singleton = np.array(np.equal(class_sizes, 1))
        singleton_classes, = np.where(is_class_singleton)
        if len(singleton_classes):
            warn('There are {} singleton classes that will be ignored during '
                 'training. A copy of the inputs will be made.'
                 .format(len(singleton_classes)))
            is_sample_single = np.asarray([yi in singleton_classes for yi in
                                           y_inverse])
            # -1 is used by semi-supervised algorithms
            # y_inverse[is_sample_single] = -2
            X = X[~is_sample_single].copy()
            y_inverse = y_inverse[~is_sample_single].copy()

        # Check number of non-singleton classes > 1
        n_classes_non_singleton = len(classes) - len(singleton_classes)
        if n_classes_non_singleton < 2:
            raise ValueError("LargeMarginNearestNeighbor needs at least 2 "
                             "non-singleton classes, got {}."
                             .format(n_classes_non_singleton))

        n_features = len(X[0])
        check_scalar(self.warm_start, 'warm_start', bool)
        if self.warm_start and hasattr(self, 'transformation_'):
            if set(classes) != set(self.classes_non_singleton_):
                raise ValueError("warm_start can only be used where `y` has "
                                 "the same classes as in the previous call "
                                 "to fit. Previously got {}, `y` has {}"
                                 .format(self.classes_non_singleton_, classes))

            if len(self.transformation_[0]) != n_features:
                raise ValueError('The new inputs dimensionality ({}) does not '
                                 'match the previously learned transformation '
                                 'input dimensionality ({}).'
                                 .format(len(self.transformation_[0]),
                                         n_features))

        self.classes_non_singleton_ = classes[~is_class_singleton]

        if self.n_features_out is not None:
            check_scalar(self.n_features_out, 'n_features_out', int, 1)
        check_scalar(self.n_neighbors, 'n_neighbors', int, 1, len(X) - 1)
        check_scalar(self.max_iter, 'max_iter', int, 1)
        check_scalar(self.max_constraints, 'max_constraints', int, 1)
        check_scalar(self.max_corrections, 'max_corrections', int, 1)
        check_scalar(self.n_jobs, 'n_jobs', int, -1)
        check_scalar(self.tol, 'tol', float, 0.)
        check_scalar(self.init_pca, 'init_pca', bool)
        check_scalar(self.use_sparse, 'use_sparse', bool)
        check_scalar(self.verbose, 'verbose', int, 0)

        if self.callback is not None:
            if not callable(self.callback):
                raise ValueError('callback is not callable.')

        # Check linear transformation dimensions
        if self.init_transformation is not None:
            check_array(self.init_transformation)
            if len(self.init_transformation[0]) != n_features:
                raise ValueError('Transformation input dimensionality ({}) '
                                 'must match the inputs dimensionality ({}).'
                                 .format(len(self.init_transformation[0]),
                                         n_features))

            if len(self.init_transformation) > \
                    len(self.init_transformation[0]):
                raise ValueError('Transformation output dimensionality ({}) '
                                 'cannot be greater than the '
                                 'transformation input dimensionality ({}).'.
                                 format(len(self.init_transformation),
                                        len(self.init_transformation[0])))

        # Check preferred output dimensionality
        if self.n_features_out is not None:
            if self.init_transformation is not None:
                if self.n_features_out != len(self.init_transformation):
                    raise ValueError('Preferred outputs dimensionality ({}) '
                                     'does not match the given linear '
                                     'transformation {}!'.format(
                                        self.n_features_out,
                                        len(self.init_transformation)))

            elif self.n_features_out > n_features:
                raise ValueError('Preferred outputs dimensionality ({}) '
                                 'cannot be greater than the given data '
                                 'dimensionality {}!'.format(
                                    self.n_features_out, n_features))

        # Check preferred number of neighbors
        min_non_singleton_size = class_sizes[~is_class_singleton].min()
        if self.n_neighbors >= min_non_singleton_size:
            warn("n_neighbors (={}) is not less than the number of samples in "
                 "the smallest non-singleton class (={}). n_neighbors will be "
                 "set to (min_non_singleton_size - 1) for estimation."
                 .format(self.n_neighbors, min_non_singleton_size))

        self.n_neighbors_ = min(self.n_neighbors, min_non_singleton_size - 1)

        return X, y_inverse

    def _init_transformer(self, X):
        """Initialize the linear transformation by setting to user specified
        parameter, loading from a file, applying PCA or setting to identity.

        Parameters
        ----------
        X : array, shape (n_samples, n_features_in)
            Data samples.

        Returns
        -------
        transformation : array, shape (n_features_out, n_features_in)
            The initial linear transformation.
        """

        if self.init_transformation is not None:
            transformation = np.asarray(self.init_transformation)
        elif self.warm_start and hasattr(self, 'transformation_'):
            transformation = self.transformation_
        elif self.init_pca and X.shape[1] > 1:
            pca = PCA(random_state=self.random_state_)
            if self.verbose:
                print('Finding principal components... ', end='')
                sys.stdout.flush()
                t = time.time()

            pca.fit(X)

            if self.verbose:
                print('done in {:5.2f}s'.format(time.time() - t))

            transformation = pca.components_
        else:
            transformation = np.eye(X.shape[1])

        if self.n_features_out is None:
            n_features_out = transformation.shape[0]
        else:
            n_features_out = self.n_features_out
            if transformation.shape[0] > n_features_out:
                warn('Decreasing the initial linear transformation output '
                     'dimensionality ({}) to the preferred output '
                     'dimensionality ({}).'.format(transformation.shape[0],
                                                   n_features_out),
                     DataDimensionalityWarning)
                transformation = transformation[:n_features_out]

        return transformation, n_features_out

    def _select_target_neighbors(self, X, y):
        """Find the target neighbors of each sample, that stay fixed during
        training.

        Parameters
        ----------
        X : array, shape (n_samples, n_features_in)
            The training samples.

        y : array, shape (n_samples,)
            The corresponding training labels indices.

        Returns
        -------
        target_neighbors: array, shape (n_samples, n_neighbors)
            An array of neighbors indices for each sample.
        """

        if self.verbose:
            print('Finding the target neighbors... ', end='')
            sys.stdout.flush()
            t = time.time()

        target_neighbors = np.empty((X.shape[0], self.n_neighbors_), dtype=int)

        nn = NearestNeighbors(n_neighbors=self.n_neighbors_,
                              algorithm=self.algorithm, n_jobs=self.n_jobs)

        for class_id in self.classes_non_singleton_:
            ind_class, = np.where(np.equal(y, class_id))
            nn.fit(X[ind_class])
            neigh_ind = nn.kneighbors(return_distance=False)
            target_neighbors[ind_class] = ind_class[neigh_ind]

        if self.verbose:
            print('done in {:5.2f}s'.format(time.time() - t))

        return target_neighbors

    @staticmethod
    def _compute_grad_static(X, targets, verbose=False):
        """Compute the gradient component due to the target neighbors that
        stays fixed throughout training

        Parameters
        ----------
        X : array, shape (n_samples, n_features_in)
            The training samples.

        targets : array, shape (n_samples, n_neighbors)
            The k nearest neighbors of each sample from the same class.

        verbose : bool, optional (default=False)
            Whether to print progress info.

        Returns
        -------
        array, shape (n_features_in, n_features_in)
            An array with the sum of all weighted outer products.
        """
        if verbose:
            print('Computing static part of the gradient...')

        n_samples, n_neighbors = targets.shape
        row = np.repeat(range(n_samples), n_neighbors)
        col = targets.ravel()
        targets_sparse = csr_matrix((np.ones(targets.size), (row, col)),
                                    shape=(n_samples, n_samples))

        return sum_outer_products(X, targets_sparse)

    def _lbfgs_callback(self, transformation):
        self.n_iter_ += 1
        if self.verbose:
            print('{:>10}'.format(self.n_iter_ + 1))
        if self.callback is not None:
            self.callback(transformation, self.n_iter_)

    def _loss_grad(self, transformation, X, y, targets, grad_static):
        """Compute the loss under a given ``transformation`` and the
        loss gradient w.r.t. ``transformation``.

        Parameters
        ----------
        transformation : array, shape (n_features_out * n_features_in,)
            The current (flattened) linear transformation.

        X : array-like, shape (n_samples, n_features_in)
            The training samples.

        y : array, shape (n_samples,)
            The corresponding training labels.

        targets : array, shape (n_samples, n_neighbors)
            The target neighbors of each sample.

        grad_static : array, shape (n_features_in, n_features_in)
            The gradient component caused by target neighbors, that stays
            fixed throughout the algorithm.

        Returns
        -------
        loss: float
            The new loss.
        grad: array, shape (n_features_out * n_features_in,)
            The new (flattened) gradient.
        """

        n_samples, n_features_in = X.shape
        self.transformation_ = transformation.reshape(self.n_features_out_,
                                                      n_features_in)

        tic = time.time()
        Lx = self.transform(X, check_input=False)

        # Compute distances to target neighbors (plus margin)
        dist_tn = np.zeros((n_samples, self.n_neighbors_))
        for k in range(self.n_neighbors_):
            dist_tn[:, k] = row_norms(Lx - Lx[targets[:, k]], True) + 1

        # Compute distances to impostors under the current transformation
        margin_radii = dist_tn[:, -1] + 1
        imp_row, imp_col, dist_imp = \
            self._find_impostors(Lx, y, margin_radii, self.use_sparse)

        loss = 0
        shape = (n_samples, n_samples)
        A0 = csr_matrix(shape)
        for k in reversed(range(self.n_neighbors_)):
            loss1 = np.maximum(dist_tn[imp_row, k] - dist_imp, 0)
            ac, = np.where(loss1 > 0)
            A1 = csr_matrix((2*loss1[ac], (imp_row[ac], imp_col[ac])), shape)

            loss2 = np.maximum(dist_tn[imp_col, k] - dist_imp, 0)
            ac, = np.where(loss2 > 0)
            A2 = csc_matrix((2*loss2[ac], (imp_row[ac], imp_col[ac])), shape)

            values = np.squeeze(np.asarray(A1.sum(1).ravel() + A2.sum(0)))
            A0 = A0 - A1 - A2 + csr_matrix((values, (range(n_samples),
                                                     targets[:, k])), shape)
            loss += loss1.dot(loss1) + loss2.dot(loss2)

        grad_new = sum_outer_products(X, A0)
        grad = self.transformation_.dot(grad_static + grad_new)
        grad *= 2
        metric = self.transformation_.T.dot(self.transformation_)
        loss = loss + (grad_static * metric).sum()

        toc = time.time()
        self.n_funcalls_ += 1
        if self.verbose:
            print('{:10} {:>10} {:>15.6e} {:>10.2f}'
                  .format('', self.n_funcalls_, loss, toc - tic))
            sys.stdout.flush()

        return loss, grad.ravel()

    def _find_impostors(self, Lx, y, margin_radii, use_sparse=True):
        """Compute all impostor pairs exactly.

        Parameters
        ----------
        Lx : array, shape (n_samples, n_features_out)
            An array of transformed samples.

        y : array, shape (n_samples,)
            The corresponding class labels.

        margin_radii : array, shape (n_samples,)
            Distances to the farthest target neighbors + margin.

        use_sparse : bool, optional (default=True)
            Whether to use a sparse matrix for storing the impostor pairs.

        Returns
        -------
        imp_row : array, shape (n_impostors,)
            Sample indices.
        imp_col : array, shape (n_impostors,)
            Corresponding sample indices that violate a margin.
        dist : array, shape (n_impostors,)
            dist[i] is the distance between samples imp1[i] and imp2[i].
        """
        n_samples = Lx.shape[0]

        if use_sparse:
            # Initialize impostors matrix
            impostors_sp = csr_matrix((n_samples, n_samples), dtype=np.int8)

            for class_id in self.classes_non_singleton_[:-1]:
                ind_in, = np.where(np.equal(y, class_id))
                ind_out, = np.where(np.greater(y, class_id))

                # Subdivide ind_out x ind_in to chunks of a size that is
                # fitting in memory
                ii, jj = _find_impostors_batch(Lx[ind_out], Lx[ind_in],
                                               margin_radii[ind_out],
                                               margin_radii[ind_in])

                if len(ii):
                    # sample constraints if they are too many
                    if len(ii) > self.max_constraints:
                        dims = (len(ind_out), len(ind_in))
                        ind = np.ravel_multi_index((ii, jj), dims=dims)
                        ind_sampled = self.random_state_.choice(
                            ind, self.max_constraints, replace=False)
                        ii, jj = np.unravel_index(ind_sampled, dims=dims)

                    imp_row = ind_out[ii]
                    imp_col = ind_in[jj]
                    new_imp = csr_matrix((np.ones(len(imp_row), dtype=np.int8),
                                          (imp_row, imp_col)), dtype=np.int8,
                                         shape=(n_samples, n_samples))
                    impostors_sp = impostors_sp + new_imp

            impostors_sp = impostors_sp.tocoo(copy=False)
            imp_row = impostors_sp.row
            imp_col = impostors_sp.col
            dist = paired_distances_batch(Lx, imp_row, imp_col)
        else:
            # Initialize impostors vectors
            imp_row, imp_col, dist = [], [], []
            for class_id in self.classes_non_singleton_[:-1]:
                ind_in, = np.where(np.equal(y, class_id))
                ind_out, = np.where(np.greater(y, class_id))

                # Subdivide ind_out x ind_in to chunks of a size that is
                # fitting in memory
                ii, jj, dd = _find_impostors_batch(
                    Lx[ind_out], Lx[ind_in], margin_radii[ind_out],
                    margin_radii[ind_in], return_distance=True)

                if len(ii):
                    # sample constraints if they are too many
                    if len(ii) > self.max_constraints:
                        dims = (len(ind_out), len(ind_in))
                        ind = np.ravel_multi_index((ii, jj), dims=dims)
                        ind_sampled = self.random_state_.choice(
                            len(ind), self.max_constraints, replace=False)
                        dd = np.asarray(dd)[ind_sampled]
                        ind_sampled = ind[ind_sampled]
                        ii, jj = np.unravel_index(ind_sampled, dims=dims)

                    imp_row.extend(ind_out[ii])
                    imp_col.extend(ind_in[jj])
                    dist.extend(dd)

            imp_row, imp_col = np.asarray(imp_row), np.asarray(imp_col)
            dist = np.asarray(dist)

        return imp_row, imp_col, dist


##########################
# Some helper functions #
#########################


def _find_impostors_batch(X_out, X_in, margin_radii_out, margin_radii_in,
                          return_distance=False, mem_budget=int(1e7)):
    """Find impostor pairs in chunks to avoid large memory usage

    Parameters
    ----------
    X_out : array, shape (n_samples_out, n_features_out)
        An array of transformed data samples from multiple classes.

    X_in : array, shape (n_samples_in, n_features_out)
        Transformed data samples from one class not present in X_out,
        so probably n_samples_in < n_samples_out.

    margin_radii_out : array, shape (n_samples_out,)
        Distances of the samples in ``X_out`` to their margins.

    margin_radii_in : array, shape (n_samples_in,)
        Distances of the samples in ``X_in`` to their margins.

    mem_budget : int, optional (default=int(1e7))
        Memory budget (in bytes) for computing distances.

    return_distance : bool, optional (default=False)
        Whether to return the distances to the impostors.

    Returns
    -------
    imp_row : array, shape (n_impostors,)
        Sample indices.
    imp_col : array, shape (n_impostors,)
        Corresponding sample indices that violate a margin.
    dist : array, shape (n_impostors,), optional
        dist[i] is the distance between samples imp_row[i] and imp_col[i].
    """

    n_samples_out = X_out.shape[0]
    bytes_per_row = X_in.shape[0] * X_in.itemsize
    batch_size = int(mem_budget // bytes_per_row)

    imp_row, imp_col, dist = [], [], []

    # X_in squared norm stays constant, so pre-compute it to get a speed-up
    X_in_norm_squared = row_norms(X_in, squared=True)
    for chunk in gen_batches(n_samples_out, batch_size):

        #dist_out_in = euclidean_distances(X_out[chunk], X_in, squared=True,
        #                                  Y_norm_squared=X_in_norm_squared)
        # check_input in every chunk would add an extra ~8% time of computation

        X = X_out[chunk]
        XX = row_norms(X, squared=True)[:, np.newaxis]
        YY = X_in_norm_squared
        dist_out_in = safe_sparse_dot(X, X_in.T, dense_output=True)
        dist_out_in *= -2
        dist_out_in += XX
        dist_out_in += YY

        ind1, = np.where((dist_out_in < margin_radii_out[chunk, None]).ravel())
        ind2, = np.where((dist_out_in < margin_radii_in[None, :]).ravel())
        ind = np.unique(np.concatenate((ind1, ind2)))

        if len(ind):
            ii, jj = np.unravel_index(ind, dist_out_in.shape)
            imp_row.extend(ii + chunk.start)
            imp_col.extend(jj)
            if return_distance:
                # This np.maximum would add another ~8% time of computation
                np.maximum(dist_out_in, 0, out=dist_out_in)
                dist.extend(dist_out_in[ii, jj])

    if return_distance:
        return imp_row, imp_col, dist
    else:
        return imp_row, imp_col


def check_scalar(x, name, dtype, min_val=None, max_val=None):
    """Validates scalar parameters by checking if their datatype matches and
    if their values are within a valid given range.

    Parameters
    ----------
    x : object
        The scalar parameter to validate.

    name : str
        The name of the parameter to be printed in error messages.

    dtype : type
        The desired datatype for the parameter.

    min_val : float or int, optional (default=None)
        The minimum value value the parameter can take. If None (default) it
        is implied that the parameter does not have a lower bound.

    max_val: float or int, optional (default=None)
        The maximum valid value the parameter can take. If None (default) it
        is implied that the parameter does not have an upper bound.

    Raises
    -------
    TypeError
        If the parameter's type does not match the desired type.

    ValueError
        If the parameter's value violates the given bounds.
    """

    if type(x) is not dtype:
        raise TypeError('{} must be {}.'.format(name, dtype))

    if min_val is not None and x < min_val:
        raise ValueError('{} must be >= {}.'.format(name, min_val))

    if max_val is not None and x > max_val:
        raise ValueError('{} must be <= {}.'.format(name, max_val))


def sum_outer_products(X, weights):
    """Computes the sum of weighted outer products using a sparse weights
    matrix

    Parameters
    ----------
    X : array, shape (n_samples, n_features_in)
        An array of data samples.

    weights : csr_matrix, shape (n_samples, n_samples)
        A sparse weights matrix (indicating target neighbors).


    Returns
    -------
    sum_outer_prods : array, shape (n_features_in, n_features_in)
        The sum of all weighted outer products.
    """

    weights_sym = weights + weights.T
    diag = spdiags(weights_sym.sum(1).ravel(), 0, *weights_sym.shape)
    laplacian = diag - weights_sym
    sum_outer_prods = X.T.dot(laplacian.dot(X))

    return sum_outer_prods


def paired_distances_batch(X, ind_a, ind_b, mem_budget=int(1e7)):
    """Equivalent to  row_norms(X[ind_a] - X[ind_b], True)

    Parameters
    ----------
    X : array, shape (n_samples, n_features_in)
        An array of data samples.
    ind_a : array, shape (n_indices,)
        An array of sample indices.
    ind_b : array, shape (n_indices,)
        Another array of sample indices.
    mem_budget : int, optional (default=int(1e7))
        Memory budget (in bytes) for computing distances.
    Returns
    -------
    dist: array, shape (n_indices,)
        An array of pairwise distances.
    """

    bytes_per_row = X.shape[1] * X.itemsize
    batch_size = int(mem_budget // bytes_per_row)

    n_pairs = len(ind_a)
    dist = np.zeros(n_pairs)
    for chunk in gen_batches(n_pairs, batch_size):
        dist[chunk] = row_norms(X[ind_a[chunk]] - X[ind_b[chunk]], True)

    return dist
