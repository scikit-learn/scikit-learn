from math import log

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from lvq.glvq import GlvqModel
from sklearn.utils.multiclass import unique_labels, check_classification_targets

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import validation
from sklearn.utils.validation import check_is_fitted

import matplotlib.pyplot as plt


class GrlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0, initial_relevances=None):
        super().__init__(random_state, initial_prototypes, initial_rototype_labels, prototypes_per_class,
                         display, max_iter, gtol)
        self.regularization = regularization
        self.initial_relevances = initial_relevances

    def optfun(self, variables, training_data, label_equals_prototype, random_state, lr_relevances=0,
               lr_prototypes=1, calc_gradient=True):
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w_.shape[0]
        omegaT = variables[nb_prototypes:]

        dist = self._compute_distance(training_data, variables[:nb_prototypes],
                                      np.diag(omegaT))  # change dist function ?
        # dist = cdist(training_data, prototypes, 'sqeuclidean')
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)
        pidxwrong = d_wrong.argmin(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)
        pidxcorrect = d_correct.argmin(1)

        distcorrectpluswrong = distcorrect + distwrong
        distcorectminuswrong = distcorrect - distwrong
        mu = distcorectminuswrong / distcorrectpluswrong

        regTerm = self.regularization * log(np.linalg.det(omegaT.conj().T.dot(omegaT)))

        f = mu.sum(0) - regTerm

        if not calc_gradient:
            return f

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
                Gw -= np.dot(difw * dcd[np.newaxis].T, omegaT).T.dot(difw) + np.dot(difc * dwd[np.newaxis].T,
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
            G[nb_prototypes:] = 2 / n_data * lr_relevances * Gw - self.regularization * f3
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes * G[:nb_prototypes].dot(omegaT.dot(omegaT.T))
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return f, G.ravel()

    def _optimize(self, X, y, random_state):
        nb_prototypes, nb_features = self.w_.shape
        if self.initial_relevances is None:
            self.lambda_ = np.ones([nb_features])
        else:
            self.lambda_ = validation.column_or_1d(self.initial_relevances)
        variables = np.append(self.w_, np.diag(self.lambda_), axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype, random_state=random_state,
                                  lr_prototypes=1, lr_relevances=0),
            x0=variables, jac=True,
            options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype, random_state=random_state,
                                  lr_prototypes=0, lr_relevances=1),
            x0=res.x, jac=True,
            options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        self.lambda_ = np.diag(out[nb_prototypes:])
        return n_iter

    def _compute_distance(self, X, w=None, lambda_=None):
        if w is None:
            w = self.w_
        if lambda_ is None:
            lambda_ = self.lambda_
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        distance = np.zeros([nb_prototypes, nb_samples])
        for i in range(nb_prototypes):
            delta = X - w[i]
            distance[i] = np.sum(delta ** 2 * lambda_, 1)
        return np.transpose(distance)