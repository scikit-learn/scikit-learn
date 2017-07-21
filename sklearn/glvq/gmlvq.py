from math import log

import math
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from sklearn.glvq.glvq import GlvqModel
from sklearn.utils.multiclass import unique_labels, check_classification_targets

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import validation
from sklearn.utils.validation import check_is_fitted

import matplotlib.pyplot as plt


def _squared_euclidean(A, B=None):
    if B is None:
        d = np.sum(A ** 2, 1)[np.newaxis].T + np.sum(A ** 2, 1) - 2 * A.dot(A.T)
    else:
        d = np.sum(A ** 2, 1)[np.newaxis].T + np.sum(B ** 2, 1) - 2 * A.dot(B.T)
    return np.maximum(d, 0)


class GmlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0, initial_matrix=None, dim=None,
                 nb_reiterations=100):
        super().__init__(random_state, initial_prototypes, initial_rototype_labels, prototypes_per_class,
                         display, max_iter, gtol)
        if not isinstance(regularization, int):
            raise ValueError("nb_reiterations must be a int")
        self.regularization = regularization
        self.initial_matrix = initial_matrix
        self.initialdim = dim
        self.nb_reiterations = nb_reiterations

    def optgrad(self, variables, training_data, label_equals_prototype, random_state, lr_relevances=0, lr_prototypes=1):
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w_.shape[0]
        omegaT = variables[nb_prototypes:].conj().T
        dist = _squared_euclidean(training_data.dot(omegaT), variables[:nb_prototypes].dot(omegaT))
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
                Gw = Gw - np.dot(difw * dcd[np.newaxis].T, omegaT).T.dot(difw) + \
                     np.dot(difc * dwd[np.newaxis].T, omegaT).T.dot(difc)
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
            G[nb_prototypes:] = 2 / n_data * lr_relevances * Gw - self.regularization * f3
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes * G[:nb_prototypes].dot(omegaT.dot(omegaT.T))
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return G.ravel()

    def optfun(self, variables, training_data, label_equals_prototype):
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w_.shape[0]
        omegaT = variables[nb_prototypes:]  # .conj().T

        # dist = self._compute_distance(training_data, variables[:nb_prototypes],
        #                          np.diag(omegaT))  # change dist function ?
        dist = _squared_euclidean(training_data.dot(omegaT), variables[:nb_prototypes].dot(omegaT))
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
            regTerm = self.regularization * log(np.linalg.det(omegaT.conj().T.dot(omegaT)))
            return mu.sum(0) - regTerm  # f
        return mu.sum(0)

    def _optimize(self, X, y, random_state):
        nb_prototypes, nb_features = self.w_.shape
        if self.initialdim is None:
            self.dim_ = nb_features
        elif not isinstance(self.initialdim, int) or self.initialdim <= 0:
            raise ValueError("dim must be an int above 0")
        else:
            self.dim_ = self.initialdim

        if not isinstance(self.nb_reiterations, int):
            raise ValueError("nb_reiterations must be a int")
        elif self.nb_reiterations < 1:
            raise ValueError("nb_reiterations must be above 0")

        if self.initial_matrix is None:
            if self.dim_ == nb_features:
                self.omega_ = np.eye(nb_features)
            else:
                self.omega_ = random_state.rand(self.dim_, nb_features) * 2 - 1
        else:
            self.omega_ = validation.check_array(self.initial_matrix)
        variables = np.append(self.w_, self.omega_, axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, random_state=random_state,
                                       lr_prototypes=1, lr_relevances=0),
            x0=variables, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, random_state=random_state,
                                       lr_prototypes=0, lr_relevances=1),
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, random_state=random_state,
                                       lr_prototypes=1, lr_relevances=1),
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        self.omega_ = out[nb_prototypes:]
        self.omega_ /= math.sqrt(np.sum(np.diag(self.omega_.T.dot(self.omega_))))
        return n_iter

    def _compute_distance(self, X, w=None, omega=None):  # catch case where omega is not initialized
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
