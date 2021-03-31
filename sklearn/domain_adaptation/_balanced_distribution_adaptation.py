# coding=utf8
"""
Balanced Distribution Adaptation (BDA) is a popular domain adaptation method
that can match the marginal and conditional probability distributions between
two domains with tunable weight to dynamically tune their relative importance.
Eventually, better performance can be achieved after using TCA for feature
transformation.

Model Features
--------------
Balanced distribution adaptation:
    Perform BDA on any two given datasets.

Examples
--------
>>> import numpy as np
>>> Xs, Xt = np.random.randn(100, 50), np.random.randn(90, 50)
>>> Ys, Yt = np.random.randint(0, 3, 100), np.random.randint(0, 3, 90)
>>> bda = BalancedDistributionAdaptation()
>>> acc, ypre, _, _ = bda.fit_predict(Xs, Ys, Xt, Yt)
>>> BalancedDistributionAdaptation(...)

Notes
-----
References:
[1] Wang J, Chen Y, Hao S, et al. Balanced distribution adaptation for
transfer learning. 2017 IEEE international conference on data mining
(ICDM).IEEE, 2017: 1129-1134.
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
            K = pairwise.linear_kernel(np.asarray(X1).T, np.asarray(X2).T)
        else:
            K = pairwise.linear_kernel(np.asarray(X1).T)
    elif ker == 'rbf':
        if X2 is not None:
            K = pairwise.rbf_kernel(np.asarray(X1).T, np.asarray(X2).T, gamma)
        else:
            K = pairwise.rbf_kernel(np.asarray(X1).T, None, gamma)
    return K


class BalancedDistributionAdaptation:
    """Class for balanced distribution adaptation.

    Parameters
    ----------
    kernel_type : {'primal', 'linear', 'rbf'}, default='primal'
        Kernal type used for feature mapping.

    dim : int, default=30
        Dimension of the data after TCA transformation.

    lamb : float, default=1
        Regularization tradeoff parameter.

    mu : float, default=0.5
        Relative importance between marginal and conditional distributions.

    gamma : float, default=1
        Parameter for rbf kernel.

    T : int, default=10
        Number of iterations to solve the problem.

    mode: {'BDA', 'WBDA'}, default='BDA'
        Model for plain BDA or weighted BDA for class imbalanced learning.
    """

    def __init__(self, kernel_type='primal', dim=30, lamb=1, mu=0.5,
                 gamma=1, T=10,
                 mode='BDA'):
        self.kernel_type = kernel_type
        self.dim = dim
        self.lamb = lamb
        self.mu = mu
        self.gamma = gamma
        self.T = T
        self.mode = mode

    def fit_predict(self, Xs, Ys, Xt, Yt, clf):
        ''' Fit the data using BDA and use a classifier for classification.
        Xs : array, (n1 x d), n1 is the number of rows and d is the dimension.
            Raw data from the first domain.

        Ys : array, (n1 x 1), n1 is the number of rows.
            Raw labels from the first domain.

        Xt : array, (n2 x d), n2 is the number of rows and d is the dimension.
            Raw data from the second domain.

        Yt : array, (n2 x 1), n2 is the number of rows.
            Raw labels from the second domain.

        clf: classifier
            Classifier used for classification.

        Returns
        -------
        acc : float
            Classification accuracy of TCA using a certain classifier.

        Y_tar_pseudo : array, (n2 x 1), where n2 is the number of rows.
            Predicted label for the second domain.

        list_acc : list
            Predicted accuracy for the target domain in all iterations.

        A : array, (n1+n2) x (n1+n2)
            Transformation matrix for BDA.

        '''
        list_acc = []
        X = np.hstack((Xs.T, Xt.T))
        X /= np.linalg.norm(X, axis=0)
        m, n = X.shape
        ns, nt = len(Xs), len(Xt)
        e = np.vstack((1 / ns * np.ones((ns, 1)), -1 / nt * np.ones((nt, 1))))
        C = len(np.unique(Ys))
        H = np.eye(n) - 1 / n * np.ones((n, n))
        mu = self.mu
        M = 0
        Y_tar_pseudo = None
        Xs_new = None
        for t in range(self.T):
            N = 0
            M0 = e * e.T * C
            if Y_tar_pseudo is not None and len(Y_tar_pseudo) == nt:
                for c in range(1, C + 1):
                    e = np.zeros((n, 1))
                    Ns = len(Ys[np.where(Ys == c)])
                    Nt = len(Y_tar_pseudo[np.where(Y_tar_pseudo == c)])

                    if self.mode == 'WBDA':
                        Ps = Ns / len(Ys)
                        Pt = Nt / len(Y_tar_pseudo)
                        alpha = Pt / Ps
                        mu = 1
                    else:
                        alpha = 1

                    tt = Ys == c
                    e[tt] = 1 / Ns
                    yy = Y_tar_pseudo == c
                    ind = np.where(yy)
                    inds = [item + ns for item in ind]
                    e[tuple(inds)] = -alpha / Nt
                    e[np.isinf(e)] = 0
                    N = N + np.dot(e, e.T)
            M = (1 - mu) * M0 + mu * N
            M /= np.linalg.norm(M, 'fro')
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

            clf.fit(Xs_new, Ys.ravel())
            Y_tar_pseudo = clf.predict(Xt_new)
            acc = accuracy_score(Yt, Y_tar_pseudo)
            list_acc.append(acc)
        return acc, Y_tar_pseudo, list_acc, A
