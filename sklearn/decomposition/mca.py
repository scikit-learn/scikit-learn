# coding=utf-8

""" Multiple Correspondance Analysis
not possible without Vivek Yadak's help
http://vxy10.github.io/2016/06/10/intro-MCA/
"""

# TODO: Handle Sparse Arrays?
# TODO: ~~handle benzecri correction~~
# TODO: ~~Handle greenacre correction~~
# TODO: ~~derive n_components if not provided?~~
# TODO: ~~Add documentation~~
# TODO: ~~Add lots of comments~~
# TODO: Add BaseClass
# TODO: add .score methods
# TODO: ~~Make Expl_var accessible~~
# TODO: SVD svd_solver parameter

from scipy.linalg import diagsvd
from numpy.linalg import svd
import numpy as np

import functools


def Matrix_mult(*args):
    """An internal method to multiply matrices."""
    return functools.reduce(np.dot, args)


class MCA():
    """ Multiple Correspondance Analysis

    Dimensionality reduction for categorical, one-hot encoded data.

    Singular Value Decomposition is applied to the categorical
    data to get linear rows or factors that describe the input
    data.

    Parameters
    ----------
    n_components: int
        Number of dimensions to reduce the input data to

    correction: string {'auto', 'benzecri', 'greenacre'}
        Since categorical data inputs are not uncorrelated, apply
        a correction algorithm to improve the variances between
        output data dimensions

    K: int
        The number of true independent measures (pre-encoding
        columns) for Benzécri or Greenacre corrections

    J: int
        The number of categorical variables for Greenacre
        correction


    Attributes
    ----------
    explained_variance_ : The amount of variance explained by each
    of the selected components.

    References
    ----------

    Vivek Yadak, "Multiple Correspondance Analysis: Principal
        Component Analysis
        for Catergorical variables"
        See http://vxy10.github.io/2016/06/10/intro-MCA/

    Examples
    --------
    >>> from sklearn.decomposition import MCA
    >>> from sklearn.preprocessing import OneHotEncoder
    >>> from numpy import array
    >>> from numpy import argmax
    >>> from sklearn.preprocessing import LabelEncoder
    >>> from sklearn.preprocessing import OneHotEncoder
    >>> data = ['bitter', 'bitter', 'spicy', 'bitter',
        'sweet', 'sweet', 'spicy', 'bitter', 'spicy', 'sweet']
    >>> values = array(data)
    >>> label_encoder = LabelEncoder()
    >>> integer_encoded = label_encoder.fit_transform(values)
    >>> onehot_encoder = OneHotEncoder(sparse=False)
    >>> integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    >>> onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
    >>> mca = MCA(n_components=3)
    >>> mca.fit_transform(X=onehot_encoded)
    array([[-0.12247449,  0.        ,  0.1       ],
           [-0.12247449,  0.        ,  0.1       ],
           [ 0.08164966,  0.12909944,  0.1       ],
           [-0.12247449,  0.        ,  0.1       ],
           [ 0.08164966, -0.12909944,  0.1       ],
           [ 0.08164966, -0.12909944,  0.1       ]])
    >>> r.explained_variance_
    array([0.5, 0.5, 0. ])

    """

    def __init__(self, n_components=None, correction='auto', K=10, J=None):
        self.n_components = n_components
        self.correction = correction
        self.explained_variance_ = None
        self.K = K
        self.J = J

    def fit(self, X, y=None):
        """
        X : array-like, shape (n_samples, n_features)
        """
        self._fit(X)
        return self

    def _fit(self, X, y=None):
        # determine if we can call full _fit
        # if smaller than 500 dimensions

        return self._fit_full(X, self.n_components)

    def _fit_full(self, X, n_components):
        n_samples, n_features = X.shape
        # i_sup = X.tail(1)

        if self.J is None:
            self.J = n_features

        N_all = np.sum(X)
        Z = X / N_all

        # Get 2 vectors corresponding to the sum of rows and colums.
        Sum_r = np.sum(Z, axis=1)
        Sum_c = np.sum(Z, axis=0)

        # Compute residual matrix by subtracting the expected indicator matrix
        # (outer product of rows and columns sums computed in step 2)
        Z_expected = np.outer(Sum_r, Sum_c)
        Z_residual = Z - Z_expected

        # Scale residual by the square root of column and row sums.
        D_r_sqrt_mi = np.sqrt(np.diag(Sum_r**-1))
        D_c_sqrt_mi = np.sqrt(np.diag(Sum_c**-1))

        MCA_mat = Matrix_mult(D_r_sqrt_mi, Z_residual, D_c_sqrt_mi)
        P, S, Q = svd(MCA_mat)  # most costly part, whole matrices set to false
        S_d = diagsvd(S, n_samples, n_features)

        G = Matrix_mult(D_c_sqrt_mi, Q.T, S_d.T)
        # Row space, contains linear combinations of rows

        Lam = S**2
        Expl_var = Lam/np.sum(Lam)
        # Explained variance before correction

        # Benzecri correction
        if self.correction == 'auto' or self.correction == 'benzecri':
            K = self.K
            E_ = [(K/(K-1.)*(lm - 1./K))**2 if lm > 1./K else 0 for lm in S**2]
            E = np.array(E_)
            Expl_var_bn = E/np.sum(E)
        elif self.correction == 'greenacre':
            K = self.K
            J = self.J
            green_norm = (K / (K - 1.) * (np.sum(S**4) - (J - K) / K**2.))
            Expl_var_bn = E/green_norm
        elif self.correction == 'none':
            Expl_var_bn = Expl_var
        else:
            raise ValueError("Unknown correction %s" % (self.correction))

        self.explained_variance_ = Expl_var_bn
        # explained variance after correction

        dim = self.n_components

        X_pjn = []
        for i in np.arange(0, 6):
            X_pjn.append(np.dot(X[i], G[:, :dim])/S[:dim]/10)

        X_pjn = np.asarray(X_pjn)

        # THE RESULT
        return X_pjn

    def fit_transform(self, X, y=None):
        X = self._fit(X)
        return X
