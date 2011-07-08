""" Dictionary learning
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD


import numpy as np

from .sparse_pca import dict_learning, dict_learning_online, \
                        _update_code_parallel
from ..base import BaseEstimator, TransformerMixin
from ..linear_model import orthogonal_mp
from ..metrics.pairwise import euclidean_distances


class BaseDictionaryLearning(BaseEstimator, TransformerMixin):
    """ Dictionary learning base class
    """
    def __init__(self, n_atoms, transform_method='omp', split_sign=False):
        self.n_atoms = n_atoms
        self.transform_method = transform_method
        self.split_sign = split_sign

    def transform(self, X, y=None, **kwargs):
        """Encode the data as a sparse combination of the learned dictionary
        atoms.

        Coding method is determined by the object parameter `transform_method`.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        TODO: document kwargs for each possible coding method

        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data
        """
        # XXX : kwargs is not documented

        # XXX: parameters should be made explicit so we can have defaults
        if self.transform_method == 'omp':
            code = orthogonal_mp(self.components_.T, X.T, **kwargs).T
        elif self.transform_method in ('lasso_cd', 'lasso_lars'):
            code = _update_code_parallel(self.components_.T, X.T, **kwargs).T

        # XXX: threshold and triangle are not verified to be correct
        elif self.transform_method == 'threshold':
            alpha = float(kwargs['alpha'])
            code = np.dot(X, self.components_.T)
            code = np.sign(code) * np.maximum(np.abs(code) - alpha, 0)
        elif self.transform_method == 'triangle':
            distances = euclidean_distances(X, self.components_)
            distance_means = distances.mean(axis=1)[:, np.newaxis]
            code = np.maximum(0, distance_means - distances)
        else:
            raise NotImplemented('Coding method %s is not implemented' %
                                 self.transform_method)

        if self.split_sign:
            # feature vector is split into a positive and negative side
            n_samples, n_features = code.shape
            split_code = np.empty((n_samples, 2 * n_features))
            split_code[:, :n_features] = np.maximum(code, 0)
            split_code[:, n_features:] = -np.minimum(code, 0)
            code = split_code

        return code


class DictionaryLearning(BaseDictionaryLearning):
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

    transform_method: 'lasso_lars' | 'lasso_cd' | 'omp' | 'threshold' |
                      'triangle'
        method to use for transforming the data after the dictionary has been
        learned

    coding_method: 'lars' | 'cd',
        method to use for solving the lasso problem

    n_jobs: int,
        number of parallel jobs to run

    code_init: array of shape (n_samples, n_atoms),
        initial value for the code, for warm restart

    dict_init: array of shape (n_atoms, n_features),
        initial values for the dictionary, for warm restart

    verbose:
        degree of verbosity of the printed output

    Attributes
    ----------
    components_: array, [n_atoms, n_features]
        dictionary atoms extracted from the data

    error_: array
        vector of errors at each iteration

    References
    ----------
    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)


    See also
    --------
    `scikits.learn.decomposition.SparsePCA`

    """
    def __init__(self, n_atoms, alpha=1, max_iter=1000, tol=1e-8,
                 transform_method='omp', coding_method='lars', n_jobs=1,
                 code_init=None, dict_init=None, verbose=False,
                 split_sign=False):
        self.n_atoms = n_atoms
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.transform_method = transform_method
        self.coding_method = coding_method
        self.n_jobs = n_jobs
        self.code_init = code_init
        self.dict_init = dict_init
        self.verbose = verbose
        self.split_sign = split_sign

    def fit_transform(self, X, y=None, **params):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        code: array-like, shape (n_samples, n_atoms)
            The sparse encoding of the data, coded using the same method used
            in the fit. To transform data using a different sparse coding
            technique such as `OMP`, see the `transform` method.
        """
        self._set_params(**params)
        X = np.asanyarray(X)
        V, U, E = dict_learning(X, self.n_atoms, self.alpha,
                                tol=self.tol, max_iter=self.max_iter,
                                method=self.coding_method,
                                n_jobs=self.n_jobs,
                                code_init=self.code_init,
                                dict_init=self.dict_init,
                                verbose=self.verbose)
        self.components_ = U
        self.error_ = E
        return V

    def fit(self, X, y=None, **params):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self: object
            Returns the object itself
        """
        self.fit_transform(X, y, **params)
        return self


class DictionaryLearningOnline(BaseDictionaryLearning):
    """ Online dictionary learning

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

    n_iter: int,
        total number of iterations to perform

    transform_method: 'lasso_lars' | 'lasso_cd' | 'omp' | 'threshold' |
                      'triangle'
        method to use for transforming the data after the dictionary has been
        learned

    coding_method: 'lars' | 'cd',
        method to use for solving the lasso problem

    n_jobs: int,
        number of parallel jobs to run

    dict_init: array of shape (n_atoms, n_features),
        initial value of the dictionary for warm restart scenarios

    verbose:
        degree of verbosity of the printed output

    chunk_size: int,
        number of samples in each mini-batch

    shuffle: bool,
        whether to shuffle the samples before forming batches

    Attributes
    ----------
    components_: array, [n_atoms, n_features]
        components extracted from the data

    References
    ----------
    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)


    See also
    --------
    `scikits.learn.decomposition.SparsePCA`

    """
    def __init__(self, n_atoms, alpha=1, n_iter=1000, coding_method='lars',
                 n_jobs=1, chunk_size=3, shuffle=True, dict_init=None,
                 transform_method='omp', verbose=False, split_sign=False):
        self.n_atoms = n_atoms
        self.alpha = alpha
        self.n_iter = n_iter
        self.coding_method = coding_method
        self.transform_method = transform_method
        self.n_jobs = n_jobs
        self.dict_init = dict_init
        self.verbose = verbose
        self.shuffle = shuffle
        self.chunk_size = chunk_size
        self.split_sign = split_sign

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
        U = dict_learning_online(X, self.n_atoms, self.alpha,
                                 n_iter=self.n_iter,
                                 coding_method=self.coding_method,
                                 n_jobs=self.n_jobs, dict_init=self.dict_init,
                                 chunk_size=self.chunk_size,
                                 shuffle=self.shuffle, verbose=self.verbose,
                                 return_code=False)
        self.components_ = U
        return self

    def partial_fit(self, X, y=None, iter_offset=0, **params):
        """Updates the model using the data in X as a mini-batch.

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
        X = np.atleast_2d(X)
        if hasattr(self, 'components_'):
            dict_init = self.components_
        else:
            dict_init = self.dict_init
        U = dict_learning_online(X, self.n_atoms, self.alpha,
                                 n_iter=self.n_iter,
                                 coding_method=self.coding_method,
                                 n_jobs=self.n_jobs,
                                 dict_init=dict_init,
                                 chunk_size=len(X), shuffle=False,
                                 verbose=self.verbose, return_code=False,
                                 iter_offset=iter_offset)
        self.components_ = U
        return self
