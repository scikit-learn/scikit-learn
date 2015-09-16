# -*- coding: utf-8 -*-
"""
@author: thiolliere and Yuan Tang (terrytangyuan)
"""
import numpy as np
import scipy.optimize as opt
from sklearn.base import BaseEstimator

class NCAcost(object):

    @staticmethod
    def cost(A, X, y, threshold=None):
        """Compute the cost function and the gradient
        This is the objective function to be minimized

        Parameters:
        -----------
        A : array-like
            Projection matrix, shape = [dim, n_features] with dim <= n_features
        X : array-like
            Training data, shape = [n_features, n_samples]
        y : array-like
            Target values, shape = [n_samples]

        Returns:
        --------
        f : float
            The value of the objective function
        gradf : array-like
            The gradient of the objective function, shape = [dim * n_features]
        """

        (D, N) = np.shape(X)
        A = np.reshape(A, (np.size(A) / np.size(X, axis=0), np.size(X, axis=0)))
        (d, aux) = np.shape(A)
        assert D == aux

        AX = np.dot(A, X)
        normAX = np.linalg.norm(AX[:, :, None] - AX[:, None, :], axis=0)

        denomSum = np.sum(np.exp(-normAX[:, :]), axis=0)
        Pij = np.exp(- normAX) / denomSum[:, None]
        if threshold is not None:
            Pij[Pij < threshold] = 0
            Pij[Pij > 1-threshold] = 1

        mask = (y != y[:, None])
        Pijmask = np.ma.masked_array(Pij, mask)
        P = np.array(np.sum(Pijmask, axis=1))
        mask = np.negative(mask)

        f = np.sum(P)

        Xi = X[:, :, None] - X[:, None, :]
        Xi = np.swapaxes(Xi, 0, 2)

        Xi = Pij[:, :, None] * Xi

        Xij = Xi[:, :, :, None] * Xi[:, :, None, :]

        gradf = np.sum(P[:, None, None] * np.sum(Xij[:], axis=1), axis=0)

        # To optimize (use mask ?)
        for i in range(N):
            aux = np.sum(Xij[i, mask[i]], axis=0)
            gradf -= aux

        gradf = 2 * np.dot(A, gradf)
        gradf = -np.reshape(gradf, np.size(gradf))
        f = np.size(X, 1) - f

        return [f, gradf]

    @staticmethod
    def f(A, X, y):
        return cost(A, X, y)[0]

    @staticmethod
    def grad(A, X, y):
        return cost(A, X, y)[1]

    @staticmethod
    def cost_g(A, X, y, threshold=None):
        """Compute the cost function and the gradient for the K-L divergence

        Parameters:
        -----------
        A : array-like
            Projection matrix, shape = [dim, n_features] with dim <= n_features
        X : array-like
            Training data, shape = [n_features, n_samples]
        y : array-like
            Target values, shape = [n_samples]

        Returns:
        --------
        g : float
            The value of the objective function
        gradg : array-like
            The gradient of the objective function, shape = [dim * n_features]
        """

        (D, N) = np.shape(X)
        A = np.reshape(A, (np.size(A) / np.size(X, axis=0), np.size(X, axis=0)))
        (d, aux) = np.shape(A)
        assert D == aux

        AX = np.dot(A, X)
        normAX = np.linalg.norm(AX[:, :, None] - AX[:, None, :], axis=0)

        denomSum = np.sum(np.exp(-normAX[:, :]), axis=0)
        Pij = np.exp(- normAX) / denomSum[:, None]
        if threshold is not None:
            Pij[Pij < threshold] = 0
            Pij[Pij > 1-threshold] = 1

        mask = (y != y[:, None])
        Pijmask = np.ma.masked_array(Pij, mask)
        P = np.array(np.sum(Pijmask, axis=1))
        mask = np.negative(mask)

        g = np.sum(np.log(P))

        Xi = X[:, :, None] - X[:, None, :]
        Xi = np.swapaxes(Xi, 0, 2)

        Xi = Pij[:, :, None] * Xi

        Xij = Xi[:, :, :, None] * Xi[:, :, None, :]

        gradg = np.sum(np.sum(Xij[:], axis=1), axis=0)

        # To optimize (use mask ?)
        for i in range(N):
            aux = np.sum(Xij[i, mask[i]], axis=0) / P[i]
            gradg -= aux

        gradg = 2 * np.dot(A, gradg)
        gradg = -np.reshape(gradg, np.size(gradg))
        g = -g

        return [g, gradg]


class NCA(BaseEstimator):

    def __init__(self, metric=None, dim=None,
                 threshold=None, objective='mahalanobis', **kwargs):
        """Classification and/or dimensionality reduction with the neighborhood
        component analysis.

        The algorithm apply the softmax function on the transformed space and
        tries to maximise the leave-one-out classification.

        Parameters:
        -----------
        metric : array-like, optional
            The initial distance metric, if not precised, the algorithm will
            use a poor projection of the Mahalanobis distance.
            shape = [dim, n_features] with dim <= n_features being the
            dimension of the output space
        dim : int, optional
            The number of dimensions to keep for dimensionality reduction. If
            not precised, the algorithm wont perform dimensionality reduction.
        threshold : float, otpional
            Threshold for the softmax function, set it higher to discard
            further neighbors.
        objective : string, optional
            The objective function to optimize. The two implemented cost
            functions are for Mahalanobis distance and KL-divergence.
        **kwargs : keyword arguments, optional
            See scipy.optimise.minimize for the list of additional arguments.
            Those arguments include:

            method : string
                The algorithm to use for optimization.
            options : dict
                a dictionary of solver options
            hess, hessp : callable
                Hessian matrix
            bounds : sequence
                Bounds for variables
            constraints : dict or sequence of dict
                Constraints definition
            tol : float
                Tolerance for termination

        Attributes:
        -----------
        metric : array-like
            The trained disctance metric
        """
        self.metric = metric
        self.dim = dim
        self.threshold = threshold
        if objective == 'mahalanobis':
            self.objective = NCAcost.cost
        elif objective == 'kl-divergence':
            self.objective = NCAcost.cost_g
        self.kwargs = kwargs

    def fit(self, X, y):
        """Fit the model using X as training data and y as target values.

        Parameters:
        -----------
        X : array-like
            Training data, shape = [n_features, n_samples]
        y : array-like
            Target values, shape = [n_samples]
        """
        if self.metric is None:
            if self.dim is None:
                self.metric = np.eye(np.size(X, 1))
                self.dim = np.size(X, 1)
            else:
                self.metric = np.eye(self.dim, np.size(X, 1) - self.dim)

        res = opt.minimize(fun=self.objective,
                           x0=self.metric,
                           args=(X, y, self.threshold),
                           jac=True,
                           **self.kwargs
                           )

        self.metric = np.reshape(res.x,
                                 (np.size(res.x) / np.size(X, 0),
                                  np.size(X, 0)))

    def fit_transform(self, X, y):
        """Fit the model with X and apply the dimensionality reduction on X.

        Parameters:
        -----------
        X : array-like
            Training data, shape = [n_features, n_samples]
        y : array-like
            Target values, shape = [n_samples]

        Returns:
        --------
        X_new : array-like
            shape = [dim, n_samples]
        """
        self.fit(self, X, y)
        return np.dot(self.metric, X)

    def score(self, X, y):
        """Returns the proportion of X correctly classified by the leave-one-
        out classification

        Parameters:
        -----------
        X : array-like
            Training data, shape = [n_features, n_samples]
        y : array-like
            Target values, shape = [n_samples]

        Returns:
        --------
        score : float
            The proportion of X correctly classified
        """
        return 1 - NCAcost.cost(self.metric, X, y)[0]/np.size(X, 1)

    def getParameters(self):
        """Returns a dictionary of the parameters
        """
        return dict(metric=self.metric, dim=self.dim, objective=self.objective,
                    threshold=self.threshold, **self.kwargs)
