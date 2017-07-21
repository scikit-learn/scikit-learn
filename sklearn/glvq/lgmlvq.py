from math import log

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


class LgmlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0, initial_matrices=None, classwise=False,
                 dim=None,
                 nb_reiterations=100):
        super().__init__(random_state, initial_prototypes, initial_rototype_labels, prototypes_per_class,
                         display, max_iter, gtol)
        if not isinstance(regularization, int):
            raise ValueError("nb_reiterations must be a int")
        self.regularization = regularization
        self.initial_matrices = initial_matrices
        self.classwise = classwise
        self.initialdim = dim
        self.nb_reiterations = nb_reiterations

    def g(self, variables, training_data, label_equals_prototype, random_state, lr_relevances=0,
          lr_prototypes=1):
        nb_samples, nb_features = training_data.shape
        nb_prototypes = self.c_w_.shape[0]
        variables = variables.reshape(variables.size // nb_features, nb_features)
        # dim to indices
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        psis = np.split(variables[nb_prototypes:], indices[:-1])#.conj().T

        dist = self._compute_distance(training_data, variables[:nb_prototypes], psis)  # change dist function ?
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

        G = np.zeros(variables.shape)
        normfactors = 4 / distcorrectpluswrong ** 2

        if lr_relevances > 0:
            Gw = []
            for i in range(len(psis)):
                Gw.append(np.zeros(psis[i].shape))
        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong
            if self.classwise:
                rightIdx = np.where((self.c_w_(i) == self.classes_) == 1)[0]  # test if works
            else:
                rightIdx = i
            dcd = distcorrect[idxw] * normfactors[idxw]
            dwd = distwrong[idxc] * normfactors[idxc]

            difc = training_data[idxc] - variables[i]
            difw = training_data[idxw] - variables[i]
            if lr_prototypes > 0:
                G[i] = (dcd.dot(difw) - dwd.dot(difc)).dot(psis[rightIdx].conj().T).dot(psis[rightIdx])
            if lr_relevances > 0:
                Gw[rightIdx] = Gw[rightIdx] - (difw * dcd[np.newaxis].T).dot(psis[rightIdx].conj().T).T.dot(difw) + \
                               (difc * dwd[np.newaxis].T).dot(psis[rightIdx].conj().T).T.dot(difc)
        if lr_relevances > 0:
            if sum(self.regularization) > 0:
                regmatrices = np.zeros(sum(self.dim_), nb_features)
                for i in range(len(psis)):
                    regmatrices[sum(self.dim_[:i]) - self.dim_[i] + 1:sum(self.dim_[:i])] = self.regularization[
                                                                                                i] * np.linalg.pinv(
                        psis[i])
                G[nb_prototypes:] = 2 / nb_samples * lr_relevances * np.concatenate(Gw) - regmatrices
            else:
                G[nb_prototypes:] = 2 / nb_samples * lr_relevances * np.concatenate(Gw)
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / nb_samples * lr_prototypes * G[:nb_prototypes]
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return G.ravel()

    def f(self, variables, training_data, label_equals_prototype):
        nb_samples, nb_features = training_data.shape
        nb_prototypes = self.c_w_.shape[0]
        variables = variables.reshape(variables.size // nb_features, nb_features)
        # dim to indices
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        psis = np.split(variables[nb_prototypes:], indices[:-1])#.conj().T

        dist = self._compute_distance(training_data, variables[:nb_prototypes], psis)  # change dist function ?
        # dist = cdist(training_data, prototypes, 'sqeuclidean')
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)

        distcorrectpluswrong = distcorrect + distwrong
        distcorectminuswrong = distcorrect - distwrong
        mu = distcorectminuswrong / distcorrectpluswrong

        if sum(self.regularization) > 0:
            func = np.vectorize(lambda x: np.log(np.linalg.det(x * x.conj().T)))
            regTerm = self.regularization * func(psis)
            return mu - 1 / nb_samples * regTerm(pidxcorrect) - 1 / nb_samples * regTerm(pidxwrong)
        return mu.sum(0)

    def _optimize(self, X, y, random_state):
        nb_prototypes, nb_features = self.w_.shape
        if self.initialdim is None:
            self.dim_ = [nb_features]
        elif not isinstance(self.initialdim, list):
            raise ValueError("dim must be a list of ints")
        else:
            for i in self.initialdim:
                if not isinstance(i, int) or i <= 0:
                    raise ValueError("dim must be a list of positive ints")
            self.dim_ = self.initialdim

        if not isinstance(self.nb_reiterations, int):
            raise ValueError("nb_reiterations must be a int")
        elif self.nb_reiterations < 1:
            raise ValueError("nb_reiterations must be above 0")

        # initialize psis (psis is list of arrays)
        nb_classes = len(self.classes_)
        if self.initial_matrices is None:
            if self.classwise:
                if len(self.dim_) == 1:
                    self.dim_ *= np.ones(nb_classes)
                else:
                    if len(self.dim_) > nb_classes:
                        self.dim_ = self.dim_[:nb_classes]
                    else:
                        self.dim_ = np.append(self.dim_, nb_features * np.ones(nb_classes - len(self.dim_)), axis=1)
            else:
                if len(self.dim_) == 1:
                    self.dim_ *= np.ones(self.w_.shape[0], dtype=np.int)
                elif len(self.dim_) > self.w_.shape[0]:
                    self.dim_ = self.dim_[self.w_.shape[0]]
                else:
                    self.dim_ = np.append(self.dim_, nb_features * np.ones(self.w_.shape[0] - len(self.dim_)), axis=1)
            self.psis_ = []
            for d in self.dim_:
                self.psis_.append(random_state.rand(d, nb_features) * 2 - 1)
        else:
            if self.classwise:
                if len(self.initial_matrices) != nb_classes:
                    raise ValueError("Length of matrices wrong\n"
                                     "found=%d\n"
                                     "expected=%d" % (len(self.initial_matrices), nb_classes))
                elif np.sum(map(lambda x: x.shape[1], self.initial_matrices)) != nb_features * len(
                        self.initial_matrices):
                    raise ValueError("Each matrix should have %d columns" % nb_features)
                self.psis_ = list(map(lambda x: validation.check_array(x), self.initial_matrices))
            elif len(self.initial_matrices) != self.w_.shape[0]:
                raise ValueError("Length of matrices wrong\n"
                                 "found=%d\n"
                                 "expected=%d" % (len(self.initial_matrices), nb_classes))
            elif np.sum([x.shape[1] for x in self.initial_matrices]) != nb_features * len(self.initial_matrices):
                raise ValueError("Each matrix should have %d columns" % nb_features)
            self.psis_ = list(map(lambda x: validation.check_array(x), self.initial_matrices))

        if isinstance(self.regularization, int) or isinstance(self.regularization, float):
            self.regularization = np.repeat(self.regularization, len(self.psis_))
        else:
            self.regularization = validation.column_or_1d(self.regularization)
            if self.classwise:
                if self.regularization.size != nb_classes:
                    raise ValueError("length of regularization must be number of classes")
            else:
                if self.regularization.size != self.w_.shape[0]:
                    raise ValueError("length of regularization must be number of prototypes")

        variables = np.append(self.w_, np.concatenate(self.psis_), axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            fun=lambda x: self.f(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.g(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=1,
                                       lr_relevances=0, random_state=random_state),
            x0=variables, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            fun=lambda x: self.f(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.g(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=0,
                                       lr_relevances=1, random_state=random_state),
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        res = minimize(
            fun=lambda x: self.f(x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.g(x, X, label_equals_prototype=label_equals_prototype, lr_prototypes=1,
                                       lr_relevances=1, random_state=random_state),
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        self.psis_ = np.split(variables[nb_prototypes:], indices[:-1])#.conj().T
        return n_iter

    def _compute_distance(self, X, w=None, psis=None):  # catch case where omega is not initialized
        if w is None:
            w = self.w_
        if psis is None:
            psis = self.psis_
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        if len(psis) == nb_prototypes:
            distance = np.zeros([nb_prototypes, nb_samples])
            for i in range(nb_prototypes):
                distance[i] = np.sum(np.dot(X - w[i], psis[i].conj().T) ** 2, 1)
            return np.transpose(distance)
        distance = np.zeros([nb_prototypes, nb_samples])
        for i in range(nb_prototypes):
            matrixIdx = self.classes_ == self.c_w_[i]
            distance[i] = np.sum(np.dot(X - w[i], psis[matrixIdx].conj().T) ** 2, 1)
        return np.transpose(distance)
