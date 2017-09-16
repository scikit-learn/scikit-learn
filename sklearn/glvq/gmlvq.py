# -*- coding: utf-8 -*-

# Author: Joris Jensen <jjensen@techfak.uni-bielefeld.de>
#
# License: BSD 3 clause

from __future__ import division

import math
from math import log

import numpy as np
from scipy.optimize import minimize

from .glvq import GlvqModel, _squared_euclidean
from ..utils import validation


class GmlvqModel(GlvqModel):
    """Generalized Matrix Learning Vector Quantization

    Parameters
    ----------

    prototypes_per_class : int or list of int, optional (default=1)
        Number of prototypes per class. Use list to specify different numbers
        per class.

    initial_prototypes : array-like,
     shape =  [n_prototypes, n_features + 1], optional
        Prototypes to start with. If not given initialization near the class
        means. Class label must be placed as last entry of each prototype

    initial_matrix : array-like, shape = [dim, n_features], optional
        Relevance matrix to start with.
        If not given random initialization for rectangular matrix and unity
        for squared matrix.

    regularization : float, optional (default=0.0)
        Value between 0 and 1. Regularization is done by the log determinant
        of the relevance matrix. Without regularization relevances may
        degenerate to zero.

    dim : int, optional (default=nb_features)
        Maximum rank or projection dimensions

    max_iter : int, optional (default=2500)
        The maximum number of iterations.

    gtol : float, optional (default=1e-5)
        Gradient norm must be less than gtol before successful
        termination of l-bfgs-b.

    display : boolean, optional (default=False)
        Print information about the bfgs steps.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------

    w_ : array-like, shape = [n_prototypes, n_features]
        Prototype vector, where n_prototypes in the number of prototypes and
        n_features is the number of features

    c_w_ : array-like, shape = [n_prototypes]
        Prototype classes

    classes_ : array-like, shape = [n_classes]
        Array containing labels.

    dim_ : int
        Maximum rank or projection dimensions

    omega_ : array-like, shape = [dim, n_features]
        Relevance matrix

    See also
    --------
    GlvqModel, GrlvqModel, LgmlvqModel
    """

    def __init__(self, prototypes_per_class=1, initial_prototypes=None,
                 initial_matrix=None, regularization=0.0,
                 dim=None, max_iter=2500, gtol=1e-5, display=False,
                 random_state=None):
        super(GmlvqModel, self).__init__(prototypes_per_class,
                                         initial_prototypes, max_iter, gtol,
                                         display, random_state)
        self.regularization = regularization
        self.initial_matrix = initial_matrix
        self.initialdim = dim

    def _optgrad(self, variables, training_data, label_equals_prototype,
                 random_state, lr_relevances=0, lr_prototypes=1):
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w_.shape[0]
        omegaT = variables[nb_prototypes:].conj().T
        dist = _squared_euclidean(training_data.dot(omegaT),
                                  variables[:nb_prototypes].dot(omegaT))
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)
        pidxwrong = d_wrong.argmin(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)
        pidxcorrect = d_correct.argmin(1)

        distcorrectpluswrong = distcorrect + distwrong

        G = np.zeros(variables.shape)
        distcorrectpluswrong = 4 / distcorrectpluswrong ** 2

        if lr_relevances > 0:
            Gw = np.zeros([omegaT.shape[0], n_dim])

        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong

            dcd = distcorrect[idxw] * distcorrectpluswrong[idxw]
            dwd = distwrong[idxc] * distcorrectpluswrong[idxc]
            if lr_relevances > 0:
                difc = training_data[idxc] - variables[i]
                difw = training_data[idxw] - variables[i]
                Gw = Gw - np.dot(difw * dcd[np.newaxis].T, omegaT).T \
                    .dot(difw) + np.dot(difc * dwd[np.newaxis].T,
                                        omegaT).T.dot(difc)
                if lr_prototypes > 0:
                    G[i] = dcd.dot(difw) - dwd.dot(difc)
            elif lr_prototypes > 0:
                G[i] = dcd.dot(training_data[idxw]) - \
                       dwd.dot(training_data[idxc]) + \
                       (dwd.sum(0) - dcd.sum(0)) * variables[i]
        f3 = 0
        if self.regularization:
            f3 = np.linalg.pinv(omegaT.conj().T).conj().T
        if lr_relevances > 0:
            G[nb_prototypes:] = 2 / n_data \
                                * lr_relevances * Gw - self.regularization * f3
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes \
                                * G[:nb_prototypes].dot(omegaT.dot(omegaT.T))
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return G.ravel()

    def _optfun(self, variables, training_data, label_equals_prototype):
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w_.shape[0]
        omegaT = variables[nb_prototypes:]  # .conj().T

        dist = _squared_euclidean(training_data.dot(omegaT),
                                  variables[:nb_prototypes].dot(omegaT))
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)

        distcorrectpluswrong = distcorrect + distwrong
        distcorectminuswrong = distcorrect - distwrong
        mu = distcorectminuswrong / distcorrectpluswrong

        if self.regularization > 0:
            regTerm = self.regularization * log(
                np.linalg.det(omegaT.conj().T.dot(omegaT)))
            return mu.sum(0) - regTerm  # f
        return mu.sum(0)

    def _optimize(self, X, y, random_state):
        if not isinstance(self.regularization,
                          float) or self.regularization < 0:
            raise ValueError("regularization must be a positive float ")
        nb_prototypes, nb_features = self.w_.shape
        if self.initialdim is None:
            self.dim_ = nb_features
        elif not isinstance(self.initialdim, int) or self.initialdim <= 0:
            raise ValueError("dim must be an positive int")
        else:
            self.dim_ = self.initialdim

        if self.initial_matrix is None:
            if self.dim_ == nb_features:
                self.omega_ = np.eye(nb_features)
            else:
                self.omega_ = random_state.rand(self.dim_, nb_features) * 2 - 1
        else:
            self.omega_ = validation.check_array(self.initial_matrix)
            if self.omega_.shape[1] != nb_features:  # TODO: check dim
                raise ValueError(
                    "initial matrix has wrong number of features\n"
                    "found=%d\n"
                    "expected=%d" % (self.omega_.shape[1], nb_features))

        variables = np.append(self.w_, self.omega_, axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        method = 'l-bfgs-b'
        res = minimize(
            fun=lambda x:
            self._optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x:
            self._optgrad(x, X, label_equals_prototype=label_equals_prototype,
                          random_state=random_state,
                          lr_prototypes=1, lr_relevances=0),
            method=method, x0=variables,
            options={'disp': self.display, 'gtol': self.gtol,
                     'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            fun=lambda x:
            self._optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x:
            self._optgrad(x, X, label_equals_prototype=label_equals_prototype,
                          random_state=random_state,
                          lr_prototypes=0, lr_relevances=1),
            method=method, x0=res.x,
            options={'disp': self.display, 'gtol': self.gtol,
                     'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        res = minimize(
            fun=lambda x:
            self._optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x:
            self._optgrad(x, X, label_equals_prototype=label_equals_prototype,
                          random_state=random_state,
                          lr_prototypes=1, lr_relevances=1),
            method=method, x0=res.x,
            options={'disp': self.display, 'gtol': self.gtol,
                     'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        self.omega_ = out[nb_prototypes:]
        self.omega_ /= math.sqrt(
            np.sum(np.diag(self.omega_.T.dot(self.omega_))))
        return n_iter

    def _compute_distance(self, X, w=None, omega=None):
        if w is None:
            w = self.w_
        if omega is None:
            omega = self.omega_
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        distance = np.zeros([nb_prototypes, nb_samples])
        for i in range(nb_prototypes):
            distance[i] = np.sum((X - w[i]).dot(omega.T) ** 2, 1)
        return distance.T

    def project(self, X, dims, print_variance_covered=False):
        """Projects the data input data X using the relevance matrix of trained
        model to dimension dim

        Parameters
        ----------
        X : array-like, shape = [n,n_features]
          input data for project
        dims : int
          dimension to project to
        print_variance_covered : boolean
          flag to print the covered variance of the projection

        Returns
        --------
        C : array, shape = [n,n_features]
            Returns predicted values.
        """
        v, u = np.linalg.eig(self.omega_.conj().T.dot(self.omega_))
        idx = v.argsort()[::-1]
        if print_variance_covered:
            print('variance coverd by projection:',
                  v[idx][:dims].sum() / v.sum() * 100)
        return X.dot(u[:, idx][:, :dims].dot(np.diag(np.sqrt(v[idx][:dims]))))
