# -*- coding: utf-8 -*-

"""
The :mod:`sklearn.polynomial_classifier` module implements Polynomial
Classifier algorithms.
"""

# Author: Christoph Hermes <hermes(at)hausmilbe(dot)net>
#
# License: BSD Style.

import numpy as np
import scipy as sp
from base import BaseEstimator, ClassifierMixin


class PC(BaseEstimator, ClassifierMixin):
    """Polynomial Classifier.

    This type of classifier maps features to polynomial space and applies
    linear regression.

    Parameters
    ----------
    degree : int, optional (default=2)
        Polynomial degree, usual are values from 1 to 3. Higher values support
        over-training and dramatically increase memory consumption.

    Attributes
    ----------
    `A_` : matrix, shape = [n_polynomials, n_classes]
        Estimated parameter matrix.

    Notes
    -----
    Its not uncommon to apply a dimension reduction method (e.g. PCA) to the
    feature space before using the Polynomial Classifier. This reduces the
    amount of internal memory by the polynomials.

    To Be Done
    ----------
    * replace direct computation of parameter matrix A by Gauss-Jordan
      algorithm with pivot selection
    * add rejection of a predicted sample by rad criterion

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1, 1], [-1, -1], [-1, 1], [1, -1]]) # XOR problem
    >>> y = np.array([1, 1, 2, 2])
    >>> clf = PC(degree=2)
    >>> clf.fit(X, y)
    >>> print clf.predict([[1, 0.5], [1, -0.5]])
    [1, 2]

    References
    ----------
    Schuermann, J.: Pattern Classification: A Unified View of Statistical and
    Neural Approaches, John Wiley & Sons, Inc., 1996
    Niemann, H.: Klassifikation von Mustern, University Erlangen-Nuernberg,
    2003
    http://www5.informatik.uni-erlangen.de/fileadmin/Persons/
    NiemannHeinrich/klassifikation-von-mustern/m00links.html
    """

    def __init__(self, degree=2):
        self.degree = degree
        pass

    def _build_pc_features(self, X):
        """Build polynomial features from X up to self.degree

        Assume Xj = [c1, c2, c3, ..., cD] with cI as the Ith feature of the jth
        element.
        Then build polynomial feature list, e.g. for self.degree=2:
        Pj = [1, c1, c2, ..., cD, c1*c1, c1*c2, ..., cD*cD]

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : ndarray, shape = [n_samples, n_polynomials]
            Polynomial features.
            n_polynomials = (n_features+self.degree)! /
                            (n_features! * self.degree!)
        """
        N, D = X.shape

        # number of polynomials
        numP = lambda n, p: (sp.factorial(n + p) /
               (sp.factorial(n) * sp.factorial(p)))
        # build up index list to combine features, -1 indicates unused feature
        I = np.zeros((numP(D, self.degree), self.degree), dtype=int) - 1
        for i in range(1, I.shape[0]):
            I[i, :] = I[i - 1, :]
            for j in range(self.degree):
                if I[i - 1, j] + 1 < D:
                    I[i, j] = I[i - 1, j] + 1
                    break
            j -= 1
            while j >= 0:
                I[i, j] = I[i, j + 1]
                j -= 1

        # use index list to build combined polynomial features P
        P = np.ones((N, numP(D, self.degree)))
        for i in range(I.shape[0]):
            for d in range(self.degree):
                if I[i, d] > -1:
                    P[:, i] = P[:, i] * X[:, I[i, d]]

        return P

    def fit(self, X, y):
        """Fit polynomials according to X, y

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """

        X = np.asarray(X)
        y = np.asarray(y)

        N, D = X.shape

        self._classes = np.unique(y)

        # create discriminant vector Y from y
        Y = np.zeros((N, len(self._classes)))
        for i in range(N):
            Y[i, self._classes == y[i]] = 1
        Y = np.matrix(Y)

        PX = np.matrix(self._build_pc_features(X))

        # direct computation of A
        self.A = np.linalg.pinv(PX.T * PX / N) * (PX.T * Y / N)

        return self

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Predicted target values for X
        """

        X = np.asarray(X)

        # build up polynomial features
        PX = np.matrix(self._build_pc_features(X))

        # apply classifier
        D = self.A.T * PX.T
        d = np.array(np.argmax(D, axis=0))[0]

        # apply class labels
        return [self._classes[d[i]] for i in range(len(d))]
