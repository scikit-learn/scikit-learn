from math import log

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import validation
from sklearn.utils.validation import check_is_fitted


def callback(xk):
    print(xk)


class GlvqModel(BaseEstimator, ClassifierMixin):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5):
        self.random_state = validation.check_random_state(random_state)
        self.w = initial_prototypes
        self.c_w = initial_rototype_labels
        if not isinstance(prototypes_per_class, int):
            raise ValueError("prototypes_per_class must be an integer")
        self.prototypes_per_class = prototypes_per_class
        if not isinstance(display, bool):
            raise ValueError("display must be a boolean")
        self.display = display
        if not isinstance(max_iter, int):
            raise ValueError("max_iter must be an integer")
        self.max_iter = max_iter
        if not isinstance(gtol, float):
            raise ValueError("gtol must be a float")
        self.gtol = gtol

    def optfun(self, variables, training_data, train_lab, lr_prototypes=1, calc_gradient=True):
        label_equals_prototype = np.equal(train_lab[np.newaxis].T, self.c_w)
        n_data, n_dim = training_data.shape
        nb_prototypes = self.c_w.size
        prototypes = variables.reshape(nb_prototypes, n_dim)

        dist = self._compute_distance(training_data, prototypes)
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

        f = mu.sum(0)

        if not calc_gradient:
            return f

        G = np.zeros(prototypes.shape)
        distcorrectpluswrong = 4 / distcorrectpluswrong ** 2

        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong

            dcd = distcorrect[idxw] * distcorrectpluswrong[idxw]
            dwd = distwrong[idxc] * distcorrectpluswrong[idxc]
            if lr_prototypes > 0:
                G[i] = dcd.dot(training_data[idxw]) - \
                       dwd.dot(training_data[idxc]) + \
                       (dwd.sum(0) - dcd.sum(0)) * prototypes[i]

        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes * G[:nb_prototypes]
        G = G * (1 + 0.0001 * self.random_state.rand(*G.shape) - 0.5)
        return f, G.ravel()

    def _validate_train_parms(self, train_set, train_lab):
        train_set = validation.as_float_array(validation.check_array(train_set))
        train_lab = validation.column_or_1d(validation.check_array(train_lab, ensure_2d=False))
        validation.check_X_y(train_set, train_lab)

        classes = np.unique(train_lab)
        nb_classes = len(classes)
        nb_samples = train_set.shape[0]
        nb_features = train_set.shape[1]

        # set prototypes per class
        if isinstance(self.prototypes_per_class, int):
            if self.prototypes_per_class < 0:
                raise ValueError("prototypes_per_class must be greater than 0")
            nb_ppc = np.ones([nb_classes], dtype='int') * self.prototypes_per_class
        else:
            nb_ppc = validation.column_or_1d(
                validation.check_array(self.prototypes_per_class, ensure_2d=False, dtype='int'))
            if nb_ppc.min() <= 0:
                raise ValueError("Values in prototypes_per_class must be greater than 0")
            if nb_ppc.shape[0] != nb_classes:
                raise ValueError("Length of prototypes per class does not fit the number of classes"
                                 "classes=%d"
                                 "length=%d" % (nb_classes, nb_ppc.shape[0]))
        # initialize prototypes
        if self.w is None or self.c_w is None:
            if self.w != self.c_w:
                raise ValueError("initialPrototypes and initialPrototypeLabels must both be defined")
            self.w = np.empty([np.sum(nb_ppc), nb_features], dtype='float')
            self.c_w = np.empty([np.sum(nb_ppc)])
            pos = 0
            for actClass in range(nb_classes):
                nb_prot = nb_ppc[actClass]
                mean = np.mean(train_set[train_lab == classes[actClass], :], 0)
                self.w[pos:pos + nb_prot] = mean + (self.random_state.rand(nb_prot, nb_features) * 2 - 1)
                self.c_w[pos:pos + nb_prot] = classes[actClass]
                pos += nb_prot
        else:
            self.w, self.c_w = validation.check_X_y(self.w, self.c_w)
            if self.w.shape != (np.sum(nb_ppc), nb_features):
                raise ValueError("The initial prototypes have wrong shape\n"
                                 "found=(%d,%d)\n"
                                 "expected=(%d,%d)" % (self.w.shape[0], self.w.shape[1], np.sum(nb_ppc), nb_features))
            if self.c_w.shape[0] != np.sum(nb_ppc):
                raise ValueError("Length of prototype labels wrong\n"
                                 "found=%d\n"
                                 "expected=%d" % (self.c_w.shape[0], np.sum(nb_ppc)))
            if set(self.c_w) != set(classes):
                raise ValueError("Prototype labels and test data classes dont match\n"
                                 "classes={}\n"
                                 "prototype labels={}\n".format(classes, self.c_w))
        return train_set, train_lab

    def fit(self, train_set, train_lab):
        train_set, train_lab = self._validate_train_parms(train_set, train_lab)
        res = minimize(lambda x: self.optfun(x, train_set, train_lab), #callback=callback,
                       x0=self.w, jac=True, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        self.w = res.x.reshape(self.w.shape)

    def _compute_distance(self, X, w=None):
        if w is None:
            w = self.w
        return cdist(X, w, 'euclidean')

    def predict(self, X):
        check_is_fitted(self, ['w', 'c_w'])
        dist = self._compute_distance(X)
        return (self.c_w[dist.argmin(1)])


class GrlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0, initial_relevances=None):
        super().__init__(random_state, initial_prototypes, initial_rototype_labels, prototypes_per_class,
                         display, max_iter, gtol)
        self.regularization = regularization
        self.lambda_ = initial_relevances

    def optfun(self, variables, training_data, train_lab, lr_relevances=0, lr_prototypes=1, regularization=0,
               calc_gradient=True):
        label_equals_prototype = np.equal(train_lab[np.newaxis].T, self.c_w)
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w.shape[0]
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

        regTerm = regularization * log(np.linalg.det(omegaT.conj().T.dot(omegaT)))

        f = mu.sum(0) - regTerm

        if not calc_gradient:
            return f

        G = np.zeros(variables.shape)
        distcorrectpluswrong = 4 / distcorrectpluswrong ** 2

        if lr_relevances > 0:
            Gw = np.zeros([omegaT.shape[0],n_dim])

        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong

            dcd = distcorrect[idxw] * distcorrectpluswrong[idxw]
            dwd = distwrong[idxc] * distcorrectpluswrong[idxc]
            if lr_relevances > 0:
                difc = training_data[idxc] - variables[i]
                difw = training_data[idxw] - variables[i]
                Gw -= np.dot(difw * dcd[np.newaxis].T, omegaT).T.dot(difw) + np.dot(difc * dwd[np.newaxis].T, omegaT).T.dot(difc)
                if lr_prototypes > 0:
                    G[i] = dcd.dot(difw) - dwd.dot(difc)
            elif lr_prototypes > 0:
                G[i] = dcd.dot(training_data[idxw]) - \
                       dwd.dot(training_data[idxc]) + \
                       (dwd.sum(0) - dcd.sum(0)) * variables[i]
        f3 = 0
        if regularization:
            f3 = np.linalg.pinv(omegaT.conj().T).conj().T
        if lr_relevances > 0:
            G[nb_prototypes:] = 2 / n_data * lr_relevances * Gw - regularization * f3
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes * G[:nb_prototypes].dot(omegaT.dot(omegaT.T))
        G = G * (1 + 0.0001 * self.random_state.rand(*G.shape) - 0.5)
        return f, G.ravel()

    def fit(self, train_set, train_lab):
        train_set, train_lab = self._validate_train_parms(train_set, train_lab)
        nb_prototypes, nb_features = self.w.shape
        if self.lambda_ is None:
            self.lambda_ = np.ones([nb_features])
        else:
            self.lambda_ = validation.column_or_1d(self.lambda_)
        variables = np.append(self.w, np.diag(self.lambda_), axis=0)
        res = minimize(lambda x: self.optfun(x, train_set, train_lab, lr_prototypes=1, lr_relevances=0),
                       x0=variables, jac=True,
                       options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        res = minimize(lambda x: self.optfun(x, train_set, train_lab, lr_prototypes=0, lr_relevances=1),
                       x0=res.x, jac=True,
                       options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w = out[:nb_prototypes]
        self.lambda_ = np.diag(out[nb_prototypes:])

    def _compute_distance(self, X, w=None, lambda_=None):
        if w is None:
            w = self.w
        if lambda_ is None:
            lambda_ = self.lambda_
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        distance = np.zeros([nb_prototypes, nb_samples])
        for i in range(nb_prototypes):
            delta = X - w[i]
            distance[i] = np.sum(delta ** 2 * lambda_, 1)
        return np.transpose(distance)


class GmlvqModel(GlvqModel):
    def __init__(self, random_state=None, initial_prototypes=None, initial_rototype_labels=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5, regularization=0, initial_matrix=None, dim=None,
                 nb_reiterations=100):
        super().__init__(random_state, initial_prototypes, initial_rototype_labels, prototypes_per_class,
                         display, max_iter, gtol)
        if not isinstance(regularization, int):
            raise ValueError("nb_reiterations must be a int")
        self.regularization = regularization
        self.omega = initial_matrix
        self.dim = dim
        if not isinstance(nb_reiterations, int):
            raise ValueError("nb_reiterations must be a int")
        elif nb_reiterations < 1:
            raise ValueError("nb_reiterations must be above 0")
        self.nb_reiterations = nb_reiterations

    def optfun(self, variables, training_data, train_lab, lr_relevances=0, lr_prototypes=1, regularization=0,
               calc_gradient=True):
        label_equals_prototype = np.equal(train_lab[np.newaxis].T, self.c_w)
        n_data, n_dim = training_data.shape
        variables = variables.reshape(variables.size // n_dim, n_dim)
        nb_prototypes = self.c_w.shape[0]
        omegaT = variables[nb_prototypes:]

        dist = self._compute_distance(training_data, variables[:nb_prototypes], omegaT)  # change dist function ?
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

        regTerm = regularization * log(np.linalg.det(omegaT.conj().T.dot(omegaT)))

        f = mu.sum(0) - regTerm

        if not calc_gradient:
            return f

        G = np.zeros(variables.shape)
        distcorrectpluswrong = 4 / distcorrectpluswrong ** 2

        if lr_relevances > 0:
            Gw = np.zeros([omegaT.shape[0],n_dim])

        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong

            dcd = distcorrect[idxw] * distcorrectpluswrong[idxw]
            dwd = distwrong[idxc] * distcorrectpluswrong[idxc]
            if lr_relevances > 0:
                difc = training_data[idxc] - variables[i]
                difw = training_data[idxw] - variables[i]
                Gw -= np.dot(difw * dcd[np.newaxis].T, omegaT).T.dot(difw) + np.dot(difc * dwd[np.newaxis].T, omegaT).T.dot(difc)
                if lr_prototypes > 0:
                    G[i] = dcd.dot(difw) - dwd.dot(difc)
            elif lr_prototypes > 0:
                G[i] = dcd.dot(training_data[idxw]) - \
                       dwd.dot(training_data[idxc]) + \
                       (dwd.sum(0) - dcd.sum(0)) * variables[i]
        f3 = 0
        if regularization:
            f3 = np.linalg.pinv(omegaT.conj().T).conj().T
        if lr_relevances > 0:
            G[nb_prototypes:] = 2 / n_data * lr_relevances * Gw - regularization * f3
        if lr_prototypes > 0:
            G[:nb_prototypes] = 1 / n_data * lr_prototypes * G[:nb_prototypes].dot(omegaT.dot(omegaT.T))
        G = G * (1 + 0.0001 * self.random_state.rand(*G.shape) - 0.5)
        return f, G.ravel()

    def fit(self, train_set, train_lab):
        train_set, train_lab = self._validate_train_parms(train_set, train_lab)
        nb_prototypes, nb_features = self.w.shape
        if self.dim is None:
            self.dim = nb_features
        elif not isinstance(self.dim, int) and self.dim <= 0:
            raise ValueError("dim must be an int above 0")
        if self.omega is None:
            if self.dim == nb_features:
                self.omega = np.eye(nb_features)
            else:
                self.omega = self.random_state.rand(self.dim, nb_features) * 2 - 1
        else:
            self.omega = validation.check_array(self.omega)
        variables = np.append(self.w, self.omega, axis=0)
        res = minimize(lambda x: self.optfun(x, train_set, train_lab, lr_prototypes=1, lr_relevances=0),
                       x0=variables, jac=True,
                       options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        res = minimize(lambda x: self.optfun(x, train_set, train_lab, lr_prototypes=0, lr_relevances=1),
                       x0=res.x, jac=True,
                       options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        out = res.x.reshape(res.x.size // nb_features, nb_features)
        self.w = out[:nb_prototypes]
        self.lambda_ = np.diag(out[nb_prototypes:])

    def _compute_distance(self, X, w=None, omega=None):  # catch case where omega is not initialized
        if w is None:
            w = self.w
        if omega is None:
            omega = self.omega
        nb_samples = X.shape[0]
        nb_prototypes = w.shape[0]
        distance = np.zeros([nb_prototypes, nb_samples])
        for i in range(nb_prototypes):
            distance[i] = np.sum(np.dot(X - w[i], omega.conj().T) ** 2, 1)
        return np.transpose(distance)
