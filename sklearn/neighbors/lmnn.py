# coding: utf-8
"""
Large Margin Nearest Neighbor Classification
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
# License: BSD 3 clause

from __future__ import print_function
from warnings import warn

import sys
import time
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, spdiags

from ..base import BaseEstimator, TransformerMixin
from ..neighbors import NearestNeighbors
from ..decomposition import PCA
from ..utils import gen_batches
from ..utils.extmath import row_norms, safe_sparse_dot
from ..utils.random import check_random_state
from ..utils.multiclass import check_classification_targets
from ..utils.validation import check_is_fitted, check_array, check_X_y


class LargeMarginNearestNeighbor(BaseEstimator, TransformerMixin):
    """Distance Metric Learning for Large Margin Classification.

    Parameters
    ----------
    n_features_out : int, optional (default=None)
        Preferred dimensionality of the embedding.
        If None it is inferred from ``init``.

    init : string or numpy array, optional (default='pca')
        Initialization of linear transformation. Possible options are 'pca',
        'identity' and a numpy array of shape (n_features_out, n_features).

    warm_start : bool, optional, (default=False)
        If True and :meth:`fit` has been called before, the solution of the
        previous call to :meth:`fit` is used as the initial linear
        transformation (``init`` is ignored).

    n_neighbors : int, optional (default=3)
        Number of neighbors to use as target neighbors for each sample.

    neighbors_algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the target neighbors, passed to
        neighbors.NearestNeighbors instance

    max_constraints : int, optional (default=500000)
        Maximum number of constraints to enforce per iteration.

    imp_store : {'auto', 'list', 'sparse'}, optional
        Data structure used to store the impostors:

        - 'list' will use 3 lists to store the indices of samples,
           their impostors and the distance between them.
        - 'sparse' will use a sparse indicator matrix to store the (sample,
          impostor) pairs. The distances to the impostors will be computed
          twice, but this option tends to be more efficient than 'list' as the
          data set size increases.
        - 'auto' will attempt to decide the most appropriate approach
          based on the values passed to :meth:`fit`.

    max_iter : int, optional (default=50)
        Maximum number of iterations in the optimization.

    tol : float, optional (default=1e-5)
        Convergence tolerance for the optimization.

    callback : callable, optional (default=None)
        If not None, this function is called after every iteration of the
        optimizer taking as arguments the current solution and the number of
        iterations. This might be useful in case one wants to examine or store
        the transformation found after each iteration.

    store_opt_result : bool, optional (default=False)
        If True, the OptimizeResult object returned by :meth:`minimize` of
        scipy.optimize will be stored in the attribute ``opt_result_``.

    verbose : int, optional (default=0)
        If 0, no progress messages will be printed.
        If 1, progress messages will be printed to stdout.
        If > 1, progress messages will be printed and the ``iprint``
        parameter of :meth:`_minimize_lbfgsb` of scipy.optimize will be set
        to ``verbose - 2``.

    random_state : int or numpy.RandomState or None, optional (default=None)
        A pseudo random number generator object or a seed for it if int.

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Doesn't affect :meth:`fit` method.

    Attributes
    ----------
    transformation_ : array, shape (n_features_out, n_features).
        The linear transformation learned during fitting.

    n_neighbors_ : int
        The provided n_neighbors is decreased if it is greater than or equal
        to  min(number of elements in each class).

    classes_inverse_non_singleton_ : array, shape (n_classes_non_singleton,)
        The appearing classes that have more than one sample, encoded as
        integers within the range(0, n_classes).

    n_iter_ : int
        Counts the number of iterations performed by the optimizer.

    opt_result_ : OptimizeResult (optional)
        If ``store_opt_result`` is True, this will be a dictionary of
        information representing the optimization result.

    Examples
    --------
    >>> from sklearn.neighbors import LargeMarginNearestNeighbor
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = load_iris(return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ... stratify=y, test_size=0.7, random_state=42)
    >>> lmnn = LargeMarginNearestNeighbor(n_neighbors=3, random_state=42)
    >>> lmnn.fit(X_train, y_train) # doctest: +ELLIPSIS
    LargeMarginNearestNeighbor(...)
    >>> knn = KNeighborsClassifier(n_neighbors=lmnn.n_neighbors_)
    >>> knn.fit(X_train, y_train) # doctest: +ELLIPSIS
    KNeighborsClassifier(...)
    >>> print(knn.score(X_test, y_test))
    0.933333333333
    >>> knn.fit(lmnn.transform(X_train), y_train) # doctest: +ELLIPSIS
    KNeighborsClassifier(...)
    >>> print(knn.score(lmnn.transform(X_test), y_test))
    0.971428571429

    Notes
    -----
    Large margin nearest neighbor (LMNN) is a machine learning algorithm for
    metric learning. It learns a (pseudo-)metric in a supervised fashion to
    improve the classification accuracy of the k-nearest neighbor rule.
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

    .. warning::

        At least for 32bit systems, one cannot expect precise reproducibility
        of PCA and therefore the transformations in 2 identical runs can
        diverge even before the first iteration of LargeMarginNearestNeighbor.
        Therefore, one should not expect any reproducibility of the
        transformations found by `LargeMarginNearestNeighbor` when
        initialization with PCA is used (``init``='pca').

    References
    ----------
    .. [1] Weinberger, Kilian Q., and Lawrence K. Saul. "Distance Metric
    Learning for Large Margin Nearest Neighbor Classification."
    Journal of Machine Learning Research, Vol. 10, Feb. 2009, pp. 207-244.
    (http://jmlr.csail.mit.edu/papers/volume10/weinberger09a/weinberger09a.pdf)

    .. [2] Wikipedia entry on Large Margin Nearest Neighbor
    (https://en.wikipedia.org/wiki/Large_margin_nearest_neighbor)

    """

    def __init__(self, n_features_out=None, init='pca', warm_start=False,
                 n_neighbors=3, neighbors_algorithm='auto',
                 max_constraints=500000, imp_store='auto', max_iter=50,
                 tol=1e-5, callback=None, store_opt_result=False, verbose=0,
                 random_state=None, n_jobs=1):

        # Parameters
        self.n_features_out = n_features_out
        self.init = init
        self.warm_start = warm_start
        self.n_neighbors = n_neighbors
        self.neighbors_algorithm = neighbors_algorithm
        self.max_constraints = max_constraints
        self.imp_store = imp_store
        self.max_iter = max_iter
        self.tol = tol
        self.callback = callback
        self.store_opt_result = store_opt_result
        self.verbose = verbose
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The training samples.

        y : array-like, shape (n_samples,)
            The corresponding training labels.

        Returns
        -------
        self : object
            returns a trained LargeMarginNearestNeighbor model.
        """

        # Validate the inputs
        X, y = check_X_y(X, y, ensure_min_samples=2)
        check_classification_targets(y)

        # Check that the inputs are consistent with the parameters
        X_valid, y_valid, init = self._validate_params(X, y)

        # Initialize the random generator
        self.random_state_ = check_random_state(self.random_state)

        # Measure the total training time
        t_train = time.time()

        # Initialize the linear transformation
        transformation = self._initialize(X_valid, init)

        # Find the target neighbors
        targets = _select_target_neighbors(X_valid, y_valid, self.n_neighbors_,
                                           verbose=self.verbose,
                                           n_jobs=self.n_jobs,
                                           algorithm=self.neighbors_algorithm)

        # Compute the gradient part contributed by the target neighbors
        grad_static = _compute_grad_static(X_valid, targets, self.verbose)

        # Decide how to store the impostors
        if self.imp_store == 'sparse':
            use_sparse = True
        elif self.imp_store == 'list':
            use_sparse = False
        else:
            # auto: Use a heuristic based on the data set size
            use_sparse = X_valid.shape[0] > 6500

        # Create a dictionary of parameters to be passed to the optimizer
        disp = self.verbose - 2 if self.verbose > 1 else -1
        optimizer_params = {'method': 'L-BFGS-B',
                            'fun': self._loss_grad_lbfgs,
                            'jac': True,
                            'args': (X_valid, y_valid, targets, grad_static,
                                     use_sparse),
                            'x0': transformation,
                            'tol': self.tol,
                            'options': dict(maxiter=self.max_iter, disp=disp),
                            'callback': self._callback
                            }

        # Call the optimizer
        self.n_iter_ = 0
        opt_result = minimize(**optimizer_params)

        # Reshape the solution found by the optimizer
        self.transformation_ = opt_result.x.reshape(-1, X_valid.shape[1])

        # Stop timer
        t_train = time.time() - t_train
        if self.verbose:
            print('Training took {:8.2f}s.'.format(t_train))

        # Optionally store information returned by the optimizer
        if self.store_opt_result:
            self.opt_result_ = opt_result

        return self

    def transform(self, X):
        """Applies the learned transformation to the given data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data samples.

        Returns
        -------
        X_embedded: array-like, shape (n_samples, n_features_out)
            The data samples transformed.

        Raises
        ------
        NotFittedError
            If :meth:`fit` has not been called before.
        """

        return self._transform(X)

    def _transform(self, X, check_input=True):
        """Applies the learned transformation to the given data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data samples.

        check_input: bool, optional (default=True)
            Whether to validate ``X``.

        Returns
        -------
        X_embedded: array-like, shape (n_samples, n_features_out)
            The data samples transformed.

        Raises
        ------
        NotFittedError
            If :meth:`fit` has not been called before.
        """

        check_is_fitted(self, ['transformation_'])

        if check_input:
            X = check_array(X)

        return np.dot(X, self.transformation_.T)

    def _validate_params(self, X, y):
        """Validate parameters as soon as :meth:`fit` is called.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The training samples.

        y : array-like, shape (n_samples,)
            The corresponding training labels.

        Returns
        -------
        X : array, shape (n_samples, n_features)
            The validated training samples.

        y_inverse : array, shape (n_samples,)
            The validated training labels, encoded to be integers in
            the range(0, n_classes).

        init : string or numpy array
            The validated initialization of the linear transformation.

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
        classes_inverse = np.arange(len(classes))

        # Ignore classes that have less than 2 samples (singleton classes)
        class_sizes = np.bincount(y_inverse)
        mask_singleton_class = class_sizes == 1
        singleton_classes, = np.where(mask_singleton_class)
        if len(singleton_classes):
            warn('There are {} singleton classes that will be ignored during '
                 'training. A copy of the inputs `X` and `y` will be made.'
                 .format(len(singleton_classes)))
            mask_singleton_sample = np.asarray([yi in singleton_classes for
                                                yi in y_inverse])
            X = X[~mask_singleton_sample].copy()
            y_inverse = y_inverse[~mask_singleton_sample].copy()

        # Check that there are at least 2 non-singleton classes
        n_classes_non_singleton = len(classes) - len(singleton_classes)
        if n_classes_non_singleton < 2:
            raise ValueError('LargeMarginNearestNeighbor needs at least 2 '
                             'non-singleton classes, got {}.'
                             .format(n_classes_non_singleton))

        self.classes_inverse_non_singleton_ = \
            classes_inverse[~mask_singleton_class]

        # Check the preferred embedding dimensionality
        if self.n_features_out is not None:
            _check_scalar(self.n_features_out, 'n_features_out', int, 1)

            if self.n_features_out > X.shape[1]:
                raise ValueError('The preferred embedding dimensionality '
                                 '`n_features_out` ({}) cannot be greater '
                                 'than the given data dimensionality ({})!'
                                 .format(self.n_features_out, X.shape[1]))

        # If warm_start is enabled, check that the inputs are consistent
        _check_scalar(self.warm_start, 'warm_start', bool)
        if self.warm_start and hasattr(self, 'transformation_'):
            if self.transformation_.shape[1] != X.shape[1]:
                raise ValueError('The new inputs dimensionality ({}) does not '
                                 'match the input dimensionality of the '
                                 'previously learned transformation ({}).'
                                 .format(self.transformation_.shape[1],
                                         X.shape[1]))

        _check_scalar(self.n_neighbors, 'n_neighbors', int, 1, X.shape[0] - 1)
        _check_scalar(self.max_iter, 'max_iter', int, 1)
        _check_scalar(self.tol, 'tol', float, 0.)
        _check_scalar(self.max_constraints, 'max_constraints', int, 1)
        _check_scalar(self.imp_store, 'imp_store', str)
        _check_scalar(self.n_jobs, 'n_jobs', int)
        _check_scalar(self.verbose, 'verbose', int, 0)

        if self.imp_store not in ['auto', 'sparse', 'list']:
            raise ValueError("`imp_store` must be 'auto', 'sparse' or 'list'.")

        if self.callback is not None:
            if not callable(self.callback):
                raise ValueError('`callback` is not callable.')

        # Check how the linear transformation should be initialized
        init = self.init
        if isinstance(init, np.ndarray):
            init = check_array(init)

            # Assert that init.shape[1] = X.shape[1]
            if init.shape[1] != X.shape[1]:
                raise ValueError('The input dimensionality ({}) of the given '
                                 'linear transformation `init` must match the '
                                 'dimensionality of the given inputs `X` ({}).'
                                 .format(init.shape[1], X.shape[1]))

            # Assert that init.shape[0] <= init.shape[1]
            if init.shape[0] > init.shape[1]:
                raise ValueError('The output dimensionality ({}) of the given '
                                 'linear transformation `init` cannot be '
                                 'greater than its input dimensionality ({}).'
                                 .format(init.shape[0], init.shape[1]))

            if self.n_features_out is not None:
                # Assert that self.n_features_out = init.shape[0]
                if self.n_features_out != init.shape[0]:
                    raise ValueError('The preferred embedding dimensionality '
                                     '`n_features_out` ({}) does not match '
                                     'the output dimensionality of the given '
                                     'linear transformation `init` ({})!'
                                     .format(self.n_features_out,
                                             init.shape[0]))
        elif init in ['pca', 'identity']:
            pass
        else:
            raise ValueError("`init` must be 'pca', 'identity', or a numpy "
                             "array of shape (n_features_out, n_features).")

        # Check the preferred number of neighbors
        min_non_singleton_size = class_sizes[~mask_singleton_class].min()
        if self.n_neighbors >= min_non_singleton_size:
            warn('`n_neighbors` (={}) is not less than the number of '
                 'samples in the smallest non-singleton class (={}). '
                 '`n_neighbors_` will be set to {} for estimation.'
                 .format(self.n_neighbors, min_non_singleton_size,
                         min_non_singleton_size-1))

        self.n_neighbors_ = min(self.n_neighbors, min_non_singleton_size-1)

        return X, y_inverse, init

    def _initialize(self, X, init):
        """

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The training samples.

        init : string or numpy array
            The initialization of the linear transformation.

        Returns
        -------
        transformation : array, shape (n_features_out, n_features)
            The initialized linear transformation.

        """

        transformation = init
        if self.warm_start and hasattr(self, 'transformation_'):
            transformation = self.transformation_

        elif isinstance(init, np.ndarray):
            pass

        elif init == 'pca':
            pca = PCA(n_components=self.n_features_out,
                      random_state=self.random_state_)
            t_pca = time.time()
            if self.verbose:
                print('Finding principal components... ', end='')
                sys.stdout.flush()

            pca.fit(X)
            if self.verbose:
                print('done in {:5.2f}s.'.format(time.time() - t_pca))

            transformation = pca.components_

        elif init == 'identity':
            if self.n_features_out is None:
                transformation = np.eye(X.shape[1])
            else:
                transformation = np.eye(self.n_features_out, X.shape[1])

        return transformation

    def _callback(self, transformation):
        """Called after each iteration of the optimizer.

        Parameters
        ----------
        transformation : array, shape(n_features_out, n_features)
            The solution computed by the optimizer in this iteration.
        """
        if self.callback is not None:
            self.callback(transformation, self.n_iter_)

        self.n_iter_ += 1

    def _loss_grad_lbfgs(self, transformation, X, y, targets, grad_static,
                         use_sparse):
        """Compute the loss and the loss gradient w.r.t. ``transformation``.

        Parameters
        ----------
        transformation : array, shape (n_features_out * n_features,)
            The current (flattened) linear transformation.

        X : array-like, shape (n_samples, n_features)
            The training samples.

        y : array, shape (n_samples,)
            The corresponding training labels.

        targets : array, shape (n_samples, n_neighbors)
            The target neighbors of each sample.

        grad_static : array, shape (n_features, n_features)
            The gradient component caused by target neighbors, that stays
            fixed throughout the algorithm.

        use_sparse : bool
            Whether to use a sparse matrix to store the impostors.

        Returns
        -------
        loss: float
            The loss based on the given transformation.

        grad: array, shape (n_features_out * n_features,)
            The new (flattened) gradient of the loss.
        """

        n_samples, n_features = X.shape
        transformation = transformation.reshape(-1, n_features)
        self.transformation_ = transformation

        if self.n_iter_ == 0:
            self.n_iter_ += 1
            if self.verbose:
                header_fields = ['Iteration', 'Objective Value',
                                 '#Constraints', 'Time(s)']
                header_fmt = '{:>10} {:>20} {:>15} {:>10}'
                header = header_fmt.format(*header_fields)
                print('\n{}\n{}'.format(header, '-' * len(header)))

        t_start = time.time()
        X_embedded = self._transform(X, check_input=False)

        # Compute squared distances to target neighbors (plus margin)
        n_neighbors = targets.shape[1]
        dist_tn = np.zeros((n_samples, n_neighbors))
        for k in range(n_neighbors):
            dist_tn[:, k] = row_norms(X_embedded - X_embedded[targets[:, k]],
                                      squared=True) + 1

        # Find the impostors and compute (squared) distances to them
        impostors_graph = \
            self._find_impostors(X_embedded, y, dist_tn[:, -1], use_sparse)

        loss, grad_new, n_constraints = \
            _compute_push_loss(X, targets, dist_tn, impostors_graph)
        grad = np.dot(transformation, grad_static + grad_new)
        grad *= 2
        metric = np.dot(transformation.T, transformation)
        loss += np.dot(grad_static.ravel(), metric.ravel())

        t = time.time() - t_start
        if self.verbose:
            values_fmt = '{:>10} {:>20.6e} {:>15,} {:>10.2f}'
            print(values_fmt.format(self.n_iter_, loss, n_constraints, t))
            sys.stdout.flush()

        return loss, grad.ravel()

    def _find_impostors(self, X_embedded, y, margin_radii, use_sparse=True):
        """Compute the (sample, impostor) pairs exactly.

        Parameters
        ----------
        X_embedded : array, shape (n_samples, n_features_out)
            An array of transformed samples.

        y : array, shape (n_samples,)
            The corresponding (possibly encoded) class labels.

        margin_radii : array, shape (n_samples,)
            Distances to the farthest target neighbors + margin.

        use_sparse : bool, optional (default=True)
            Whether to use a sparse matrix to store the pairs.

        Returns
        -------
        impostors_graph : coo_matrix, shape (n_samples, n_samples)
            Element (i, j) is the distance between samples i and j if j is an
            impostor to i, otherwise zero.

        """
        n_samples = X_embedded.shape[0]

        if use_sparse:
            # Initialize a sparse (indicator) matrix for impostors storage
            impostors_sp = csr_matrix((n_samples, n_samples), dtype=np.int8)
            for class_id in self.classes_inverse_non_singleton_[:-1]:
                ind_in, = np.where(y == class_id)
                ind_out, = np.where(y > class_id)

                # Split ind_out x ind_in into chunks of a size that fits
                # in memory
                imp_ind = _find_impostors_blockwise(X_embedded[ind_out],
                                                    X_embedded[ind_in],
                                                    margin_radii[ind_out],
                                                    margin_radii[ind_in])

                if len(imp_ind):
                    # sample constraints if they are too many
                    if len(imp_ind) > self.max_constraints:
                        imp_ind = self.random_state_.choice(
                            imp_ind, self.max_constraints, replace=False)

                    dims = (len(ind_out), len(ind_in))
                    ii, jj = np.unravel_index(imp_ind, dims=dims)
                    # Convert indices to refer to the original data matrix
                    imp_row = ind_out[ii]
                    imp_col = ind_in[jj]
                    new_imp = csr_matrix((np.ones(len(imp_row), dtype=np.int8),
                                          (imp_row, imp_col)), dtype=np.int8,
                                         shape=(n_samples, n_samples))
                    impostors_sp = impostors_sp + new_imp

            impostors_sp = impostors_sp.tocoo(copy=False)
            imp_row = impostors_sp.row
            imp_col = impostors_sp.col
            imp_dist = _paired_distances_blockwise(X_embedded, imp_row,
                                                   imp_col)
        else:
            # Initialize lists for impostors storage
            imp_row, imp_col, imp_dist = [], [], []
            for class_id in self.classes_inverse_non_singleton_[:-1]:
                ind_in, = np.where(y == class_id)
                ind_out, = np.where(y > class_id)

                # Split ind_out x ind_in into chunks of a size that fits in
                # memory
                imp_ind, dist_batch = _find_impostors_blockwise(
                    X_embedded[ind_out], X_embedded[ind_in],
                    margin_radii[ind_out], margin_radii[ind_in],
                    return_distance=True)

                if len(imp_ind):
                    # sample constraints if they are too many
                    if len(imp_ind) > self.max_constraints:
                        ind_sampled = self.random_state_.choice(
                            len(imp_ind), self.max_constraints, replace=False)
                        imp_ind = imp_ind[ind_sampled]
                        dist_batch = dist_batch[ind_sampled]

                    dims = (len(ind_out), len(ind_in))
                    ii, jj = np.unravel_index(imp_ind, dims=dims)
                    # Convert indices to refer to the original data matrix
                    ind_out_batch = ind_out[ii]
                    ind_in_batch = ind_in[jj]
                    try:
                        imp_row.extend(ind_out_batch)
                        imp_col.extend(ind_in_batch)
                        imp_dist.extend(dist_batch)
                    except TypeError:
                        imp_row.append(ind_out_batch)
                        imp_col.append(ind_in_batch)
                        imp_dist.append(dist_batch)

            imp_row = np.asarray(imp_row, dtype=int)
            imp_col = np.asarray(imp_col, dtype=int)
            imp_dist = np.asarray(imp_dist)

        impostors_graph = coo_matrix((imp_dist, (imp_row, imp_col)),
                                     shape=(n_samples, n_samples))

        return impostors_graph


##########################
# Some helper functions #
#########################


def _select_target_neighbors(X, y, n_neighbors, verbose=0, **kwargs):
    """Find the target neighbors of each data sample.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The training samples.

    y : array, shape (n_samples,)
        The corresponding training labels indices.

    n_neighbors : int
        The number of target neighbors to select for each sample in X.

    verbose : int, optional (default=0)
        Whether to print progress info.

    Returns
    -------
    target_neighbors: array, shape (n_samples, n_neighbors)
        An array of neighbors indices for each sample.
    """

    t_targets = time.time()
    if verbose:
        print('Finding the target neighbors... ', end='')
        sys.stdout.flush()

    target_neighbors = np.zeros((X.shape[0], n_neighbors), dtype=int)

    nn = NearestNeighbors(n_neighbors=n_neighbors, **kwargs)

    classes = np.unique(y)
    for class_id in classes:
        ind_class, = np.where(y == class_id)
        nn.fit(X[ind_class])
        neigh_ind = nn.kneighbors(return_distance=False)
        target_neighbors[ind_class] = ind_class[neigh_ind]

    if verbose:
        print('done in {:5.2f}s.'.format(time.time() - t_targets))

    return target_neighbors


def _compute_grad_static(X, targets, verbose=0):
    """Compute the gradient contributed by the target neighbors.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The training samples.

    targets : array, shape (n_samples, n_neighbors)
        The k nearest neighbors of each sample from the same class.

    verbose : int, optional (default=0)
        Whether to print progress info.

    Returns
    -------
    grad_targets, shape (n_features, n_features)
        An array with the sum of all outer products of samples-targets.
    """

    if verbose:
        print('Computing static part of the gradient... ', end='')

    n_samples, n_neighbors = targets.shape
    row = np.repeat(range(n_samples), n_neighbors)
    col = targets.ravel()
    targets_sparse = csr_matrix((np.ones(targets.size), (row, col)),
                                shape=(n_samples, n_samples))
    grad_targets = _sum_weighted_outer_differences(X, targets_sparse)

    if verbose:
        print('done.')

    return grad_targets


def _find_impostors_blockwise(X_a, X_b, radii_a, radii_b,
                              return_distance=False, block_size=8):
    """Find (sample, impostor) pairs in blocks to avoid large memory usage.

    Parameters
    ----------
    X_a : array, shape (n_samples_a, n_features_out)
        Transformed data samples from class A.

    X_b : array, shape (n_samples_b, n_features_out)
        Transformed data samples from class B.

    radii_a : array, shape (n_samples_a,)
        Squared distances of the samples in ``X_a`` to their margins.

    radii_b : array, shape (n_samples_b,)
        Squared distances of the samples in ``X_b`` to their margins.

    block_size : int, optional (default=8)
        The maximum number of mebibytes (MiB) of memory to use at a time for
        calculating paired distances.

    return_distance : bool, optional (default=False)
        Whether to return the squared distances to the impostors.

    Returns
    -------
    imp_indices : array, shape (n_impostors,)
        Unraveled indices of (sample, impostor) pairs referring to a matrix
        of shape (n_samples_a, n_samples_b).

    imp_distances : array, shape (n_impostors,), optional
        imp_distances[i] is the squared distance between samples imp_row[i] and
        imp_col[i], where
        imp_row, imp_col = np.unravel_index(imp_indices, dims=(n_samples_a,
        n_samples_b))
    """

    n_samples_a = X_a.shape[0]
    bytes_per_row = X_b.shape[0] * X_b.itemsize
    block_n_rows = int(block_size*1024*1024 // bytes_per_row)

    imp_indices, imp_distances = [], []

    # X_b squared norm stays constant, so pre-compute it to get a speed-up
    X_b_norm_squared = row_norms(X_b, squared=True)[np.newaxis, :]
    for chunk in gen_batches(n_samples_a, block_n_rows):
        # from sklearn.metrics.pairwise import euclidean_distances
        # distances_ab = euclidean_distances(X_a[chunk], X_b, squared=True,
        #                                    Y_norm_squared=X_b_norm_squared)
        # check_input in every chunk would add an extra ~8% time of computation

        X_a_chunk = X_a[chunk]
        X_a_norm_squared = row_norms(X_a_chunk, squared=True)[:, np.newaxis]
        distances_ab = safe_sparse_dot(X_a_chunk, X_b.T, dense_output=True)
        distances_ab *= -2
        distances_ab += X_a_norm_squared
        distances_ab += X_b_norm_squared

        ind_b, = np.where((distances_ab < radii_a[chunk, None]).ravel())
        ind_a, = np.where((distances_ab < radii_b[None, :]).ravel())
        ind = np.unique(np.concatenate((ind_a, ind_b)))

        if len(ind):
            ind_plus_offset = ind + chunk.start * X_b.shape[0]
            try:
                imp_indices.extend(ind_plus_offset)
            except TypeError:
                imp_indices.append(ind_plus_offset)

            if return_distance:
                # This np.maximum would add another ~8% time of computation
                np.maximum(distances_ab, 0, out=distances_ab)
                distances_chunk = distances_ab.ravel()[ind]
                try:
                    imp_distances.extend(distances_chunk)
                except TypeError:
                    imp_distances.append(distances_chunk)

    imp_indices = np.asarray(imp_indices)

    if return_distance:
        return imp_indices, np.asarray(imp_distances)
    else:
        return imp_indices


def _paired_distances_blockwise(X, ind_a, ind_b, squared=True, block_size=8):
    """Equivalent to row_norms(X[ind_a] - X[ind_b], squared=squared).

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        An array of data samples.

    ind_a : array, shape (n_indices,)
        An array of sample indices.

    ind_b : array, shape (n_indices,)
        Another array of sample indices.

    squared : bool (default=True)
        Whether to return the squared distances.

    block_size : int, optional (default=8)
        The maximum number of mebibytes (MiB) of memory to use at a time for
        calculating paired distances.

    Returns
    -------
    distances: array, shape (n_indices,)
        An array of pairwise, optionally squared, distances.
    """

    bytes_per_row = X.shape[1] * X.itemsize
    batch_size = int(block_size*1024*1024 // bytes_per_row)

    n_pairs = len(ind_a)
    distances = np.zeros(n_pairs)
    for chunk in gen_batches(n_pairs, batch_size):
        distances[chunk] = row_norms(X[ind_a[chunk]] - X[ind_b[chunk]], True)

    if not squared:
        np.sqrt(distances, distances)

    return distances


def _compute_push_loss(X, targets, dist_targets, impostors_graph):
    """

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        The training input samples.

    targets : array, shape (n_samples, n_neighbors)
        Indices of target neighbors of each sample.

    dist_targets : array, shape (n_samples, n_neighbors)
        Distances of samples to their target neighbors.

    impostors_graph : coo_matrix, shape (n_samples, n_samples)
        Element (i, j) is the distance between sample i and j if j is an
        impostor to i, otherwise zero.

    Returns
    -------
    loss : float
        The push loss caused by the given targets and impostors.

    grad : array, shape (n_features, n_features)
        The gradient of the push loss.

    n_constraints : int
        The number of active triplet constraints.

    """

    n_samples, n_neighbors = dist_targets.shape
    imp_row = impostors_graph.row
    imp_col = impostors_graph.col
    dist_impostors = impostors_graph.data

    loss = 0
    shape = (n_samples, n_samples)
    A0 = csr_matrix(shape)
    n_constraints = 0
    for k in range(n_neighbors-1, -1, -1):
        loss1 = np.maximum(dist_targets[imp_row, k] - dist_impostors, 0)
        ac, = np.where(loss1 > 0)
        n_constraints += len(ac)
        A1 = csr_matrix((2*loss1[ac], (imp_row[ac], imp_col[ac])), shape)

        loss2 = np.maximum(dist_targets[imp_col, k] - dist_impostors, 0)
        ac, = np.where(loss2 > 0)
        n_constraints += len(ac)
        A2 = csc_matrix((2*loss2[ac], (imp_row[ac], imp_col[ac])), shape)

        values = (A1.sum(1).ravel() + A2.sum(0)).getA1()
        A3 = csr_matrix((values, (range(n_samples), targets[:, k])), shape)
        A0 = A0 - A1 - A2 + A3
        loss += np.dot(loss1, loss1) + np.dot(loss2, loss2)

    grad = _sum_weighted_outer_differences(X, A0)

    return loss, grad, n_constraints


def _sum_weighted_outer_differences(X, weights):
    """Compute the sum of weighted outer pairwise differences.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        An array of data samples.

    weights : csr_matrix, shape (n_samples, n_samples)
        A sparse weights matrix.


    Returns
    -------
    sum_weighted_outer_diffs : array, shape (n_features, n_features)
        The sum of all outer weighted differences.
    """

    weights_sym = weights + weights.T
    diag = spdiags(weights_sym.sum(1).ravel(), 0, *weights_sym.shape)
    laplacian = diag - weights_sym

    return np.dot(X.T, safe_sparse_dot(laplacian, X, dense_output=True))


def _check_scalar(x, name, dtype, min_val=None, max_val=None):
    """Validate scalar parameters type and value.

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
        raise TypeError('`{}` must be {}.'.format(name, dtype))

    if min_val is not None and x < min_val:
        raise ValueError('`{}` must be >= {}.'.format(name, min_val))

    if max_val is not None and x > max_val:
        raise ValueError('`{}` must be <= {}.'.format(name, max_val))
