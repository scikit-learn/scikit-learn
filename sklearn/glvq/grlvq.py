from math import log

import numpy as np
from scipy.optimize import minimize

from sklearn.glvq.glvq import GlvqModel, _squared_euclidean
from sklearn.utils import validation


# TODO: implement custom optfun for grlvq without omega

class GrlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0.0, initial_relevances=None):
        super().__init__(random_state, initial_prototypes, prototypes_per_class,
                         display, max_iter, gtol)
        if not isinstance(regularization, float) or regularization<0:
            raise ValueError("regularization must be a positive float")
        self.regularization = regularization
        self.initial_relevances = initial_relevances

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
        if self.initial_relevances is None:
            self.lambda_ = np.ones([nb_features])
        else:
            self.lambda_ = validation.column_or_1d(validation.check_array(self.initial_relevances, dtype='float', ensure_2d=False))
            if self.lambda_.size != nb_features:
                raise ValueError("length of initial relevances is wrong"
                                 "features=%d"
                                 "length=%d" % (nb_features, self.lambda_.size))
        self.lambda_ /= np.sum(self.lambda_)
        variables = np.append(self.w_, np.diag(np.sqrt(self.lambda_)), axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=1,
                                       lr_relevances=0, random_state=random_state),
            method='BFGS', x0=variables,
            options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=0,
                                       lr_relevances=1, random_state=random_state),
            method='BFGS', x0=variables,
            options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        res = minimize(
            fun=lambda x: self.optfun(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=1,
                                       lr_relevances=1, random_state=random_state),
            method='BFGS', x0=variables,
            options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        self.lambda_ = np.diag(out[nb_prototypes:].T.dot(out[nb_prototypes:]))
        self.lambda_ = self.lambda_ / self.lambda_.sum()
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
        return distance.T

    def project(self, X, dims):
        idx = self.lambda_.argsort()[::-1]
        print('projection procent:', self.lambda_[idx][:dims].sum() / self.lambda_.sum())
        return X.dot(np.diag(self.lambda_)[idx][:, :dims])
