# coding=utf8
"""
Transfer Component Analysis (TCA) is a famous domain adaptation method
that can transform the different inputs from two domains into the same
underlying subspace, i.e., the reproduing kernel Hilbert space (RKHS).
If the two domains or two datasets are having totally different probability
distributions, TCA can perform domain adaptation to mitigate the domain
mismatch in two domains. Eventually, better performance can be achieved
after using TCA for feature transformation.

Model Features
--------------
Transfer component analysis:
    Perform TCA on any two given datasets.

Examples
--------
>>> import numpy as np
>>> Xs, Xt = np.random.randn(100, 50), np.random.randn(90, 50)
>>> Ys, Yt = np.random.randint(0, 3, 100), np.random.randint(0, 3, 90)
>>> tca = TransferComponentAnalysis()
>>> acc, ypre = tca.fit_predict(Xs, Ys, Xt, Yt)
>>> TransferComponentAnalysis(...)

Notes
-----
References:
[1] Pan S J, Tsang I W, Kwok J T, et al. Domain adaptation via
transfer component analysis[J].
IEEE Transactions on Neural Networks, 2010, 22(2): 199-210.
"""

# Authors: Jindong Wang <jindongwang@outlook.com>
# License: BSD

import numpy as np
import scipy.io
import scipy.linalg
from ..metrics import pairwise, accuracy_score


def kernel(X1, X2, ker='primal', gamma=1.0):
    """Kernel function used for feature mapping.

    Parameters
    ----------
    X1 : array, (n1 x d), n1 is the number of rows and d is the dimension.
        Raw data from the first domain.

    X2 : array, (n2 x d), n2 is the number of rows and d is the dimension.
        Raw data from the second domain.

    ker : {'primal', 'linear', 'rbf'}, default='primal'
        Kernal type used for feature mapping.

    gamma : float, default=1.0
        Parameter for rbf kernel.

    Returns
    -------
    K : array, (n1 + n2) x (n1 + n2)
        Tranformed data matrix.
    """
    K = None
    if not ker or ker == 'primal':
        K = X1
    elif ker == 'linear':
        if X2 is not None:
            K = pairwise.linear_kernel(np.asarray(X1).T,
                                       np.asarray(X2).T)
        else:
            K = pairwise.linear_kernel(np.asarray(X1).T)
    elif ker == 'rbf':
        if X2 is not None:
            K = pairwise.rbf_kernel(np.asarray(X1).T,
                                    np.asarray(X2).T, gamma)
        else:
            K = pairwise.rbf_kernel(np.asarray(X1).T, None, gamma)
    return K


class TransferComponentAnalysis:
    """Class for transfer component analysis.

    Parameters
    ----------
    kernel_type : {'primal', 'linear', 'rbf'}, default='primal'
        Kernal type used for feature mapping.

    dim : int, default=30
        Dimension of the data after TCA transformation.

    lamb : float, default=1
        Regularization tradeoff parameter.

    gamma : float, default=1
        Parameter for rbf kernel.
    """

    def __init__(self, kernel_type='primal', dim=30, lamb=1, gamma=1):
        self.kernel_type = kernel_type
        self.dim = dim
        self.lamb = lamb
        self.gamma = gamma

    def fit(self, Xs, Xt):
        """Fit the data using TCA algorithm.

        Parameters
        ----------
        Xs : array, (n1 x d), n1 is number of row and d is dimension.
            Raw data from the first domain.

        Xt : array, (n2 x d), n2 is number of row and d is dimension.
            Raw data from the second domain.

        Returns
        -------
        Xs_new : array, (n1 x dim), n1 is number of row and dim is dimension.
            Transformed data from the first domain.

        Xt_new : array, (n2 x d), n2 is number of row and dim is dimension.
            Transformed data from the second domain.
        """
        X = np.hstack((Xs.T, Xt.T))
        X /= np.linalg.norm(X, axis=0)
        m, n = X.shape
        ns, nt = len(Xs), len(Xt)
        e = np.vstack((1 / ns * np.ones((ns, 1)), -1 / nt * np.ones((nt, 1))))
        M = e * e.T
        M = M / np.linalg.norm(M, 'fro')
        H = np.eye(n) - 1 / n * np.ones((n, n))
        K = kernel(self.kernel_type, X, None, gamma=self.gamma)
        n_eye = m if self.kernel_type == 'primal' else n
        a, b = np.linalg.multi_dot(
            [K, M, K.T]) + self.lamb * np.eye(n_eye),
        np.linalg.multi_dot([K, H, K.T])
        w, V = scipy.linalg.eig(a, b)
        ind = np.argsort(w)
        A = V[:, ind[:self.dim]]
        Z = np.dot(A.T, K)
        Z /= np.linalg.norm(Z, axis=0)
        Xs_new, Xt_new = Z[:, :ns].T, Z[:, ns:].T
        return Xs_new, Xt_new

    def fit_predict(self, Xs, Ys, Xt, Yt, clf):
        ''' Fit the data using TCA and use a classifier for classification.
        Xs : array, (n1 x d), n1 is the number of rows and d is the dimension.
            Raw data from the first domain.

        Ys : array, (n1 x 1), n1 is the number of rows.
            Raw labels from the first domain.

        Xt : array, (n2 x d), n2 is the number of rows and d is the dimension.
            Raw data from the second domain.

        Yt : array, (n2 x 1), n2 is the number of rows.
            Raw labels from the second domain.

        Clf: classifier
            Classifier used for classification.

        Returns
        -------
        acc : float
            Classification accuracy of TCA using a certain classifier.

        y_pred : array, (n2 x 1), where n2 is the number of rows.
            Predicted label for the second domain.
        '''
        Xs_new, Xt_new = self.fit(Xs, Xt)
        if not clf:
            print('You should use a certain classifier!')
            return 0
        clf.fit(Xs_new, Ys.ravel())
        y_pred = clf.predict(Xt_new)
        acc = accuracy_score(Yt, y_pred)
        return acc, y_pred
