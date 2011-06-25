""" Dictionary learning
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD


import numpy as np

from .sparse_pca import dict_learning, sparse_pca, _update_V_parallel
from ..base import BaseEstimator, TransformerMixin
from ..linear_model import orthogonal_mp


class DictionaryLearning(BaseEstimator, TransformerMixin):
    """ Dictionary learning

    Finds a dictionary (a set of atoms) that can best be used to represent data
    using a sparse code.

    Solves the optimization problem:
    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                 (U,V)
                with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    Parameters
    ----------
    n_atoms: int,
        number of dictionary elements to extract

    alpha: int,
        sparsity controlling parameter

    max_iter: int,
        maximum number of iterations to perform

    tol: float,
        tolerance for numerical error

    method: 'batch'
        algorithm to use

    coding_method: 'lars' | 'cd',
        method to use for solving the lasso problem

    n_jobs: int,
        number of parallel jobs to run

    U_init: array of shape (n_samples, n_atoms),
    V_init: array of shape (n_atoms, n_features),
        initial values for the decomposition for warm restart scenarios

    verbose:
        degree of verbosity of the printed output

    Attributes
    ----------
    components_: array, [n_components, n_features]
        sparse components extracted from the data

    error_: array
        vector of errors at each iteration

    References
    ----------
        XXX: Todo

    See also
    --------
    `scikits.learn.decomposition.SparsePCA`

    """
    def __init__(self, n_atoms, alpha=1, max_iter=1000, tol=1e-8,
                 method='batch', coding_method='lars', n_jobs=1, U_init=None,
                 V_init=None, verbose=False):
        self.n_atoms = n_atoms
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.coding_method = coding_method
        self.n_jobs = n_jobs
        self.U_init = V_init.T if V_init else None
        self.V_init = U_init.T if U_init else None
        self.verbose = verbose

    def fit(self, X, y=None, **params):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._set_params(**params)
        X = np.asanyarray(X)
        if self.method == 'batch':
            U, _, E = sparse_pca(X.T, self.n_atoms, self.alpha,
                                 tol=self.tol, max_iter=self.max_iter,
                                 method=self.coding_method, n_jobs=self.n_jobs,
                                 U_init=self.U_init, V_init=self.V_init,
                                 verbose=self.verbose)
            self.components_ = U.T
            self.error_ = E
        else:
            raise NotImplemented('Online version will be integrated soon')

        return self

    def transform(self, X, y=None, method='omp', **kwargs):
        """Apply the projection onto the learned sparse components
        to new data.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        method: 'omp' | 'lars' | 'cd'
            Sparse coding method to use. Additional parameters are passed
            to the corresponding solver.

        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data
        """
        # XXX : kwargs is not documented

        # XXX: parameters should be made explicit so we can have defaults
        if method == 'omp':
            return orthogonal_mp(self.components_, X.T, **kwargs)
        else:
            return _update_V_parallel(self.components_.T, X.T, **kwargs)
