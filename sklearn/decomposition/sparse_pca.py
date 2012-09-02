"""Matrix factorization with Sparse PCA"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import numpy as np

from ..utils import check_random_state, array2d
from ..linear_model import ridge_regression
from ..base import BaseEstimator, TransformerMixin
from .dict_learning import dict_learning, dict_learning_online


class SparsePCA(BaseEstimator, TransformerMixin):
    """Sparse Principal Components Analysis (SparsePCA)

    Finds the set of sparse components that can optimally reconstruct
    the data.  The amount of sparseness is controllable by the coefficient
    of the L1 penalty, given by the parameter alpha.

    Parameters
    ----------
    n_components : int,
        Number of sparse atoms to extract.

    alpha : float,
        Sparsity controlling parameter. Higher values lead to sparser
        components.

    ridge_alpha : float,
        Amount of ridge shrinkage to apply in order to improve
        conditioning when calling the transform method.

    max_iter : int,
        Maximum number of iterations to perform.

    tol : float,
        Tolerance for the stopping condition.

    method : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    n_jobs : int,
        Number of parallel jobs to run.

    U_init : array of shape (n_samples, n_atoms),
        Initial values for the loadings for warm restart scenarios.

    V_init : array of shape (n_atoms, n_features),
        Initial values for the components for warm restart scenarios.

    verbose :
        Degree of verbosity of the printed output.

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Sparse components extracted from the data.

    `error_` : array
        Vector of errors at each iteration.

    See also
    --------
    PCA
    MiniBatchSparsePCA
    DictionaryLearning
    """
    def __init__(self, n_components=None, alpha=1, ridge_alpha=0.01,
            max_iter=1000, tol=1e-8, method='lars', n_jobs=1, U_init=None,
            V_init=None, verbose=False, random_state=None):
        self.n_components = n_components
        self.alpha = alpha
        self.ridge_alpha = ridge_alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.n_jobs = n_jobs
        self.U_init = U_init
        self.V_init = V_init
        self.verbose = verbose
        self.random_state = random_state

    def fit(self, X, y=None):
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
        self.random_state = check_random_state(self.random_state)
        X = array2d(X)
        if self.n_components is None:
            n_components = X.shape[1]
        else:
            n_components = self.n_components
        code_init = self.V_init.T if self.V_init is not None else None
        dict_init = self.U_init.T if self.U_init is not None else None
        Vt, _, E = dict_learning(X.T, n_components, self.alpha,
                                 tol=self.tol, max_iter=self.max_iter,
                                 method=self.method, n_jobs=self.n_jobs,
                                 verbose=self.verbose,
                                 random_state=self.random_state,
                                 code_init=code_init,
                                 dict_init=dict_init)
        self.components_ = Vt.T
        self.error_ = E
        return self

    def transform(self, X, ridge_alpha=None):
        """Least Squares projection of the data onto the sparse components.

        To avoid instability issues in case the system is under-determined,
        regularization can be applied (Ridge regression) via the
        `ridge_alpha` parameter.

        Note that Sparse PCA components orthogonality is not enforced as in PCA
        hence one cannot use a simple linear projection.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        ridge_alpha: float, default: 0.01
            Amount of ridge shrinkage to apply in order to improve
            conditioning.

        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data.
        """
        ridge_alpha = self.ridge_alpha if ridge_alpha is None else ridge_alpha
        U = ridge_regression(self.components_.T, X.T, ridge_alpha,
                             solver='dense_cholesky')
        s = np.sqrt((U ** 2).sum(axis=0))
        s[s == 0] = 1
        U /= s
        return U


class MiniBatchSparsePCA(SparsePCA):
    """Mini-batch Sparse Principal Components Analysis

    Finds the set of sparse components that can optimally reconstruct
    the data.  The amount of sparseness is controllable by the coefficient
    of the L1 penalty, given by the parameter alpha.

    Parameters
    ----------
    n_components : int,
        number of sparse atoms to extract

    alpha : int,
        Sparsity controlling parameter. Higher values lead to sparser
        components.

    ridge_alpha : float,
        Amount of ridge shrinkage to apply in order to improve
        conditioning when calling the transform method.

    n_iter : int,
        number of iterations to perform for each mini batch

    callback : callable,
        callable that gets invoked every five iterations

    chunk_size : int,
        the number of features to take in each mini batch

    verbose :
        degree of output the procedure will print

    shuffle : boolean,
        whether to shuffle the data before splitting it in batches

    n_jobs : int,
        number of parallel jobs to run, or -1 to autodetect.

    method : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Sparse components extracted from the data.

    `error_` : array
        Vector of errors at each iteration.

    See also
    --------
    PCA
    SparsePCA
    DictionaryLearning
    """
    def __init__(self, n_components=None, alpha=1, ridge_alpha=0.01,
            n_iter=100, callback=None, chunk_size=3, verbose=False,
            shuffle=True, n_jobs=1, method='lars', random_state=None):
        self.n_components = n_components
        self.alpha = alpha
        self.ridge_alpha = ridge_alpha
        self.n_iter = n_iter
        self.callback = callback
        self.chunk_size = chunk_size
        self.verbose = verbose
        self.shuffle = shuffle
        self.n_jobs = n_jobs
        self.method = method
        self.random_state = random_state

    def fit(self, X, y=None):
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
        self.random_state = check_random_state(self.random_state)
        X = array2d(X)
        if self.n_components is None:
            n_components = X.shape[1]
        else:
            n_components = self.n_components
        Vt, _ = dict_learning_online(X.T, n_components, alpha=self.alpha,
                                     n_iter=self.n_iter, return_code=True,
                                     dict_init=None, verbose=self.verbose,
                                     callback=self.callback,
                                     chunk_size=self.chunk_size,
                                     shuffle=self.shuffle,
                                     n_jobs=self.n_jobs, method=self.method,
                                     random_state=self.random_state)
        self.components_ = Vt.T
        return self
