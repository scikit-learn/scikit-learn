# -*- coding: utf-8 -*-

# Author: Joris Jensen <jjensen@techfak.uni-bielefeld.de>
#
# License: BSD 3 clause

from __future__ import division

import numpy as np
from scipy.optimize import minimize

from .glvq import GlvqModel
from ..utils import validation


class LgmlvqModel(GlvqModel):
    """Localized Generalized Matrix Learning Vector Quantization

    Parameters
    ----------

    prototypes_per_class : int or list of int, optional (default=1)
        Number of prototypes per class. Use list to specify different numbers
        per class.

    initial_prototypes : array-like, shape =  [n_prototypes, n_features + 1],
     optional
        Prototypes to start with. If not given initialization near the class
        means. Class label must be placed as last entry of each prototype.

    initial_matrices : list of array-like, optional
        Matrices to start with. If not given random initialization

    regularization : float or array-like, shape = [n_classes/n_prototypes],
     optional (default=0.0)
        Values between 0 and 1. Regularization is done by the log determinant
        of the relevance matrix. Without regularization relevances may
        degenerate to zero.

    dim : int, optional
        Maximum rank or projection dimensions

    classwise : boolean, optional
        If true, each class has one relevance matrix.
        If false, each prototype has one relevance matrix.

    max_iter : int, optional (default=2500)
        The maximum number of iterations.

    gtol : float, optional (default=1e-5)
        Gradient norm must be less than gtol before successful termination
        of l-bfgs-b.

    display : boolean, optional (default=False)
        Print information about the bfgs steps.

    random_state : int, RandomState instance or None, optional
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

    omegas_ : list of array-like
        Relevance Matrices

    dim_ : list of int
        Maximum rank of projection

    regularization_ : array-like, shape = [n_classes/n_prototypes]
        Values between 0 and 1

    See also
    --------
    GlvqModel, GrlvqModel, GmlvqModel
    """

    def __init__(self, prototypes_per_class=1, initial_prototypes=None,
                 initial_matrices=None, regularization=0.0,
                 dim=None, classwise=False, max_iter=2500, gtol=1e-5,
                 display=False, random_state=None):
        super(LgmlvqModel, self).__init__(prototypes_per_class,
                                          initial_prototypes, max_iter, gtol,
                                          display, random_state)
        self.regularization = regularization
        self.initial_matrices = initial_matrices
        self.classwise = classwise
        self.initialdim = dim

    def _g(self, variables, training_data, label_equals_prototype,
           random_state, lr_relevances=0,
           lr_prototypes=1):
        # print("g")
        nb_samples, nb_features = training_data.shape
        nb_prototypes = self.c_w_.shape[0]
        variables = variables.reshape(variables.size // nb_features,
                                      nb_features)
        # dim to indices
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        psis = np.split(variables[nb_prototypes:], indices[:-1])  # .conj().T

        dist = self._compute_distance(training_data, variables[:nb_prototypes],
                                      psis)  # change dist function ?
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
                rightIdx = np.where((self.c_w_[i] == self.classes_) == 1)[0][
                    0]  # test if works
            else:
                rightIdx = i
            dcd = distcorrect[idxw] * normfactors[idxw]
            dwd = distwrong[idxc] * normfactors[idxc]

            difc = training_data[idxc] - variables[i]
            difw = training_data[idxw] - variables[i]
            if lr_prototypes > 0:
                G[i] = (dcd.dot(difw) - dwd.dot(difc)).dot(
                    psis[rightIdx].conj().T).dot(psis[rightIdx])
            if lr_relevances > 0:
                Gw[rightIdx] = Gw[rightIdx] - (difw * dcd[np.newaxis].T).dot(
                    psis[rightIdx].conj().T).T.dot(difw) + \
                               (difc * dwd[np.newaxis].T).dot(
                                   psis[rightIdx].conj().T).T.dot(difc)
        if lr_relevances > 0:
            if sum(self.regularization_) > 0:
                regmatrices = np.zeros([sum(self.dim_), nb_features])
                for i in range(len(psis)):
                    regmatrices[sum(self.dim_[:i + 1]) - self.dim_[i]:sum(
                        self.dim_[:i + 1])] = \
                        self.regularization_[i] * np.linalg.pinv(psis[i])
                G[nb_prototypes:] = 2 / nb_samples * lr_relevances * \
                    np.concatenate(Gw) - regmatrices
            else:
                G[nb_prototypes:] = 2 / nb_samples * lr_relevances * \
                                    np.concatenate(Gw)
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / nb_samples * \
                                lr_prototypes * G[:nb_prototypes]
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return G.ravel()

    def _f(self, variables, training_data, label_equals_prototype):
        # print("f")
        nb_samples, nb_features = training_data.shape
        nb_prototypes = self.c_w_.shape[0]
        variables = variables.reshape(variables.size // nb_features,
                                      nb_features)
        # dim to indices
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        psis = np.split(variables[nb_prototypes:], indices[:-1])  # .conj().T

        dist = self._compute_distance(training_data, variables[:nb_prototypes],
                                      psis)  # change dist function ?
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

        if sum(self.regularization_) > 0:
            def test(x):
                return np.log(np.linalg.det(x.dot(x.conj().T)))

            t = np.array([test(x) for x in psis])
            regTerm = self.regularization_ * t
            return mu - 1 / nb_samples * regTerm[
                pidxcorrect] - 1 / nb_samples * regTerm[pidxwrong]
        return mu.sum(0)

    def _optimize(self, X, y, random_state):
        nb_prototypes, nb_features = self.w_.shape
        nb_classes = len(self.classes_)
        if not isinstance(self.classwise, bool):
            raise ValueError("classwise must be a boolean")
        if self.initialdim is None:
            if self.classwise:
                self.dim_ = nb_features * np.ones(nb_classes, dtype=np.int)
            else:
                self.dim_ = nb_features * np.ones(nb_prototypes, dtype=np.int)
        else:
            self.dim_ = validation.column_or_1d(self.initialdim)
            if self.dim_.size == 1:
                if self.classwise:
                    self.dim_ = self.dim_[0] * np.ones(nb_classes,
                                                       dtype=np.int)
                else:
                    self.dim_ = self.dim_[0] * np.ones(nb_prototypes,
                                                       dtype=np.int)
            elif self.classwise and self.dim_.size != nb_classes:
                raise ValueError("dim length must be number of classes")
            elif self.dim_.size != nb_prototypes:
                raise ValueError("dim length must be number of prototypes")
            if self.dim_.min() <= 0:
                raise ValueError("dim must be a list of positive ints")

        # initialize psis (psis is list of arrays)
        if self.initial_matrices is None:
            self.omegas_ = []
            for d in self.dim_:
                self.omegas_.append(
                    random_state.rand(d, nb_features) * 2.0 - 1.0)
        else:
            if not isinstance(self.initial_matrices, list):
                raise ValueError("initial matrices must be a list")
            self.omegas_ = list(map(lambda x: validation.check_array(x),
                                    self.initial_matrices))
            if self.classwise:
                if len(self.omegas_) != nb_classes:
                    raise ValueError("length of matrices wrong\n"
                                     "found=%d\n"
                                     "expected=%d" % (
                                         len(self.omegas_), nb_classes))
                elif np.sum(map(lambda x: x.shape[1],
                                self.omegas_)) != nb_features * \
                        len(self.omegas_):
                    raise ValueError(
                        "each matrix should have %d columns" % nb_features)
            elif len(self.omegas_) != nb_prototypes:
                raise ValueError("length of matrices wrong\n"
                                 "found=%d\n"
                                 "expected=%d" % (
                                     len(self.omegas_), nb_classes))
            elif np.sum([x.shape[1] for x in self.omegas_]) != \
                    nb_features * len(self.omegas_):
                raise ValueError(
                    "each matrix should have %d columns" % nb_features)

        if isinstance(self.regularization, float):
            if self.regularization < 0:
                raise ValueError('regularization must be a positive float')
            self.regularization_ = np.repeat(self.regularization,
                                             len(self.omegas_))
        else:
            self.regularization_ = validation.column_or_1d(self.regularization)
            if self.classwise:
                if self.regularization_.size != nb_classes:
                    raise ValueError(
                        "length of regularization must be number of classes")
            else:
                if self.regularization_.size != self.w_.shape[0]:
                    raise ValueError(
                        "length of regularization "
                        "must be number of prototypes")

        variables = np.append(self.w_, np.concatenate(self.omegas_), axis=0)
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            fun=lambda x: self._f(
                x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self._g(
                x, X, label_equals_prototype=label_equals_prototype,
                lr_prototypes=1, lr_relevances=0, random_state=random_state),
            method='L-BFGS-B',
            x0=variables, options={'disp': self.display, 'gtol': self.gtol,
                                   'maxiter': self.max_iter})
        n_iter = res.nit
        res = minimize(
            fun=lambda x: self._f(
                x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self._g(
                x, X, label_equals_prototype=label_equals_prototype,
                lr_prototypes=0, lr_relevances=1, random_state=random_state),
            method='L-BFGS-B',
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol,
                               'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        res = minimize(
            fun=lambda x: self._f(
                x, X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self._g(
                x, X, label_equals_prototype=label_equals_prototype,
                lr_prototypes=1, lr_relevances=1, random_state=random_state),
            method='L-BFGS-B',
            x0=res.x, options={'disp': self.display, 'gtol': self.gtol,
                               'maxiter': self.max_iter})
        n_iter = max(n_iter, res.nit)
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w_ = out[:nb_prototypes]
        indices = []
        for i in range(len(self.dim_)):
            indices.append(sum(self.dim_[:i + 1]))
        self.omegas_ = np.split(out[nb_prototypes:], indices[:-1])  # .conj().T
        return n_iter

    def _compute_distance(self, X, w=None,
                          psis=None):
        if w is None:
            w = self.w_
        if psis is None:
            psis = self.omegas_
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        distance = np.zeros([nb_prototypes, nb_samples])
        if len(psis) == nb_prototypes:
            for i in range(nb_prototypes):
                distance[i] = np.sum(np.dot(X - w[i], psis[i].conj().T) ** 2,
                                     1)
            return np.transpose(distance)
        for i in range(nb_prototypes):
            matrixIdx = np.where(self.classes_ == self.c_w_[i])[0][0]
            distance[i] = np.sum(
                np.dot(X - w[i], psis[matrixIdx].conj().T) ** 2, 1)
        return np.transpose(distance)

    def project(self, X, prototype_idx, dims, print_variance_covered=False):
        """Projects the data input data X using the relevance matrix of the
        prototype specified by prototype_idx to dimension dim

        Parameters
        ----------
        X : array-like, shape = [n,n_features]
          input data for project
        prototype_idx : int
          index of the prototype
        dims : int
          dimension to project to
        print_variance_covered : boolean
          flag to print the covered variance of the projection

        Returns
        --------
        C : array, shape = [n,n_features]
            Returns predicted values.
        """
        nb_prototypes = self.w_.shape[0]
        if len(self.omegas_) != nb_prototypes \
                or self.prototypes_per_class != 1:
            print('project only possible with classwise relevance matrix')
        # y = self.predict(X)
        v, u = np.linalg.eig(
            self.omegas_[prototype_idx].T.dot(self.omegas_[prototype_idx]))
        idx = v.argsort()[::-1]
        if print_variance_covered:
            print('variance coverd by projection:',
                  v[idx][:dims].sum() / v.sum() * 100)
        return X.dot(u[:, idx][:, :dims].dot(np.diag(np.sqrt(v[idx][:dims]))))
