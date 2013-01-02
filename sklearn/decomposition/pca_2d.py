""" 2 Dimensional Principal Component Analysis
"""

# Author :  Aristide TOSSOU <yedtoss@gmail.com>
#
# License:  BSD Style

import numpy as np
from scipy import linalg
from scipy.sparse import issparse
from ..base import BaseEstimator, TransformerMixin
from ..utils import as_float_array, assert_all_finite


class PCA2D(BaseEstimator, TransformerMixin):
    """ Two directionnal two  dimensional principal component analysis
    This technique is based on 2D matrices as opposed to standard PCA which
    is based on 1D vectors. It considers simultaneously the row and column
    directions.2D PCA is computed by examining the covariance of the data.
    It can easily get higher accuracy than 1D PCA.

    This implements the method described in `Daoqiang Zhang, Zhi-Hua Zhou, :
    Two-directional two-dimensional PCA for efficient face representation and
    recognition, Neurocomputing, Volume 69, Issues 1-3, December 2005,
    Pages 224-231`

    This implementation uses the scipy.linalg implementation of singular
    value decomposition. It works on dense matrices. The time complexity is :
    O(n_row**3 + n_column**3 + n_row**2 * n_samples + n_column**2 * n_samples)
    where n_samples is the number of examples, n_row and n_column are the
    dimension of the original 2D matrices. In practice it means that it
    can be used if n_row, n_column, n_samples are all less than 300.
    More formally, it can be used as far as one of the term in the
    complexity is not too large (less than 10**10).

    In most case, it is then significantly faster than standard 1D PCA
    which has a complexity of O( n_row**3 * n_column**3).


    Parameters
    ----------

    n_row_components : int, None, or string
        Number of components in the row direction to keep
        if n_row_components is not set or is 0, all components are kept::
            n_row_components = n_row
        if ``0 < n_row_components < 1``, select the number of components
        such that a proportion (defined by  n_row_components)
        of the variance is retained

    n_column_components : int, None, or string
        Number of components in the column direction to keep
        if n_column_components is not set or is 0, all components are kept::
            n_column_components = n_column
        if ``0 < n_column_components < 1``, select the number of components
        such that a proportion (defined by  n_column_components)
        of the variance is retained

    row_whiten : bool, optional
        When True (False by default) the `row_components` vectors are divided
        by the square root of singular values to ensure uncorrelated outputs
        with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making there data respect some hard-wired assumptions.


    column_whiten : bool, optional
        When True (False by default) the `column_components` vectors are
        divided by the square root of singular values to ensure
        uncorrelated outputs with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making there data respect some hard-wired assumptions.

    epsilon : float, optional
        (default value 1e-5). This is a regularization parameter when
        whitening is enable.Whitening involves division by the square
        root of singular values. When the singular values are close to
        zero it could lead to numerical instability. This parameter is
        added to the singular value to prevent it.

    copy : bool
        If False (True by default), data passed to fit are overwritten

    Attributes
    ----------

    `row_components_` : array, [n_row, n_row_components]
        components with maximum variance on the row direction

    `column_components_` : array, [n_column, n_column_components]
        components with maximum variance on the column direction

    See also
    --------
    PCA
    ProbabilisticPCA
    RandomizedPCA
    KernelPCA
    SparsePCA



    """
    def __init__(self, n_row_components=None, n_column_components=None,
                 row_whiten=False, column_whiten=False, epsilon=1e-5,
                 copy=True):
        self.n_row_components = n_row_components
        self.n_column_components = n_column_components
        self.row_whiten = row_whiten
        self.column_whiten = column_whiten
        self.epsilon = epsilon
        self.copy = copy

    def fit(self, X, y=None):
        """ Fit the model from data in X
        Parameters
        ----------

        X: array-like, shape (n_samples, n_row, n_column)
            Training matrix where n_samples is the number of samples
            and n_row x n_column is the dimension of the features

        Returns
        -------
        self : object
        Returns the instance itself.
        """

        # Checking if X is not sparse
        if issparse(X):
            raise TypeError("sparse matrices are not currently supported")

        # Converting the data to 3 dimensions array
        X = np.asarray(X, np.float64)
        assert_all_finite(X)
        X = np.atleast_3d(X)

        # Copy the data if necessary
        if self.copy:
            X = X.copy()

        n_samples, n_row, n_column = X.shape

        # Making sure the type of the data is float
        # X = as_float_array(X, copy=self.copy)
        # X = X.astype(np.float64)

        # Center data
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_

        # As we don't want to change the default n_row_components let's copy it

        self.n_row_components_ = np.copy(self.n_row_components)
        self.n_column_components_ = np.copy(self.n_column_components)

        # Computing Alternate 2DPCA
        if (self.n_row_components_ is not None and
                self.n_row_components_ > 0):
            U, S, V = linalg.svd(np.tensordot(X, X, axes=([0, 2], [0, 2]))
                                 / n_samples, full_matrices=False)

            if self.n_row_components_ < 1:
                self.n_row_components_ = np.sum((S / np.sum(S)).cumsum() <
                                                self.n_row_components_) + 1

            if self.row_whiten:
                # self.row_reg = np.diag(1/np.sqrt(S)+self.epsilon)

                # To save time we can multiply the whitening matrix
                # by U beforehand
                self.row_components_ = (U[:, :self.n_row_components_].
                                        dot(np.diag(1. / np.sqrt(S[:self.
                                            n_row_components_] + self.
                                                                 epsilon))))
            else:
                self.row_components_ = U[:, :self.n_row_components_]
        else:

            self.row_components_ = np.eye(n_row, n_row)

        # Computing 2DPCA
        if (self.n_column_components_ is not None and
                self.n_column_components_ > 0):
            U, S, V = linalg.svd(np.tensordot(X, X, axes=([0, 1], [0, 1]))
                                 / n_samples, full_matrices=False)

            if self.n_column_components_ < 1:
                self.n_column_components_ = np.sum((S / np.sum(S)).cumsum() <
                                                   self.n_column_components_
                                                   ) + 1

            if self.column_whiten:
                # self.column_reg = np.diag(1/np.sqrt(S)+self.epsilon)

                self.column_components_ = (U[:, :self.n_column_components_].
                                           dot(np.diag(1. / np.sqrt(S[:self.
                                               n_column_components_] +
                                               self.epsilon))))
            else:
                self.column_components_ = U[:, :self.n_column_components_]
        else:

            self.column_components_ = np.eye(n_column, n_column)
        return self

    def transform(self, X):
        """Apply the dimensionality reduction on X.
        Parameters
        ----------

        X: array-like, shape (n_samples, n_row, n_column)
            Training matrix where n_samples is the number of samples
            and n_row x n_column is the dimension of the features.

        Returns
        -------

        X_new : array-like,
                shape (n_samples, n_row_components, n_column_components)
                according to X.


        """

        # X -= self.mean_
        if issparse(X):
            raise TypeError("sparse matrices are not currently supported")

        X = np.asarray(X, np.float64)
        assert_all_finite(X)
        X = np.atleast_3d(X)

        # Disabling this features
        # if X.ndim == 2:
         #   return (self.row_components_.T.dot(X - self.mean_).
          #          dot(self.column_components_))

        # elif X.ndim == 3:
        if X.ndim == 3:
            return ((X - self.mean_).dot(self.column_components_).
                    transpose((0, 2, 1)).dot(self.row_components_).
                    transpose((0, 2, 1)))

    def inverse_transform(self, X):
        """ Transform data back to its original space, i.e.,
            return an input X_original whose transform would be X

            Parameters
            ----------

            X: array-like,
                shape (n_samples, n_row_components, n_column_components)
                New data  where n_samples is the number of samples
                and n_row_components, n_column_components are the number of
                components on the row and column direction respectively.

            Returns
            -------

            X_original : array-like, shape (n_samples, n_row, n_column)



            Notes
            -----

            If whitening is enabled, inverse_transform does not compute the
            exact inverse operation as transform.
        """

        if issparse(X):
            raise TypeError("sparse matrices are not currently supported")

        X = np.asarray(X, np.float64)
        assert_all_finite(X)
        X = np.atleast_3d(X)

        # Disabling this features
        # if X.ndim == 2:
         #   return (self.row_components_.dot(X).
          #          dot(self.column_components_.T) + self.mean_)

        # elif X.ndim == 3:
        if X.ndim == 3:

            return (X.dot(self.column_components_.T).transpose((0, 2, 1)).
                    dot(self.row_components_.T).transpose((0, 2, 1)) +
                    self.mean_)
