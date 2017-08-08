import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

from ..base import BaseEstimator, ClassifierMixin
from ..utils import validation
from ..utils.multiclass import unique_labels
from ..utils.validation import check_is_fitted


def _squared_euclidean(A, B=None):
    if B is None:
        d = np.sum(A ** 2, 1)[np.newaxis].T + np.sum(A ** 2, 1) - 2 * A.dot(A.T)
    else:
        d = np.sum(A ** 2, 1)[np.newaxis].T + np.sum(B ** 2, 1) - 2 * A.dot(B.T)
    return np.maximum(d, 0)


class GlvqModel(BaseEstimator, ClassifierMixin):
    def __init__(self, random_state=None, initial_prototypes=None, prototypes_per_class=1,
                 display=False, max_iter=2500, gtol=1e-5):
        self.random_state = random_state
        self.initial_prototypes = initial_prototypes
        self.prototypes_per_class = prototypes_per_class
        self.display = display
        self.max_iter = max_iter
        self.gtol = gtol

    def optgrad(self, variables, training_data, label_equals_prototype, random_state):
        n_data, n_dim = training_data.shape
        nb_prototypes = self.c_w_.size
        prototypes = variables.reshape(nb_prototypes, n_dim)

        dist = _squared_euclidean(training_data, prototypes)
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)
        pidxwrong = d_wrong.argmin(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)
        pidxcorrect = d_correct.argmin(1)

        distcorrectpluswrong = distcorrect + distwrong

        G = np.zeros(prototypes.shape)
        distcorrectpluswrong = 4 / distcorrectpluswrong ** 2

        for i in range(nb_prototypes):
            idxc = i == pidxcorrect
            idxw = i == pidxwrong

            dcd = distcorrect[idxw] * distcorrectpluswrong[idxw]
            dwd = distwrong[idxc] * distcorrectpluswrong[idxc]
            G[i] = dcd.dot(training_data[idxw]) - \
                   dwd.dot(training_data[idxc]) + \
                   (dwd.sum(0) - dcd.sum(0)) * prototypes[i]
        G[:nb_prototypes] = 1 / n_data * G[:nb_prototypes]
        G = G * (1 + 0.0001 * random_state.rand(*G.shape) - 0.5)
        return G.ravel()

    def optfun(self, variables, training_data, label_equals_prototype):
        n_data, n_dim = training_data.shape
        nb_prototypes = self.c_w_.size
        prototypes = variables.reshape(nb_prototypes, n_dim)

        # dist = self._compute_distance(training_data, variables[:nb_prototypes],
        #                          np.diag(omegaT))  # change dist function ?
        dist = _squared_euclidean(training_data, prototypes)
        d_wrong = dist.copy()
        d_wrong[label_equals_prototype] = np.inf
        distwrong = d_wrong.min(1)

        d_correct = dist
        d_correct[np.invert(label_equals_prototype)] = np.inf
        distcorrect = d_correct.min(1)

        distcorrectpluswrong = distcorrect + distwrong
        distcorectminuswrong = distcorrect - distwrong
        mu = distcorectminuswrong / distcorrectpluswrong

        return mu.sum(0)

    def _validate_train_parms(self, train_set, train_lab):
        random_state = validation.check_random_state(self.random_state)
        if not isinstance(self.display, bool):
            raise ValueError("display must be a boolean")
        if not isinstance(self.max_iter, int) or self.max_iter<1:
            raise ValueError("max_iter must be an positive integer")
        if not isinstance(self.gtol, float) or self.gtol<=0:
            raise ValueError("gtol must be a positive float")
        train_set, train_lab = validation.check_X_y(train_set, train_lab)

        self.classes_ = unique_labels(train_lab)
        nb_classes = len(self.classes_)
        nb_samples, nb_features = train_set.shape  # nb_samples unused

        # set prototypes per class
        if isinstance(self.prototypes_per_class, int):
            if self.prototypes_per_class < 0:
                raise ValueError("prototypes_per_class must be positive")
            nb_ppc = np.ones([nb_classes], dtype='int') * self.prototypes_per_class
        else:
            nb_ppc = validation.column_or_1d(
                validation.check_array(self.prototypes_per_class, ensure_2d=False, dtype='int'))
            if nb_ppc.min() <= 0:
                raise ValueError("values in prototypes_per_class must be positive")
            if nb_ppc.size != nb_classes:
                raise ValueError("length of prototypes per class does not fit the number of classes"
                                 "classes=%d"
                                 "length=%d" % (nb_classes, nb_ppc.size))
        # initialize prototypes
        if self.initial_prototypes is None:
            self.w_ = np.empty([np.sum(nb_ppc), nb_features], dtype=np.double)
            self.c_w_ = np.empty([nb_ppc.sum()], dtype=self.classes_.dtype)
            pos = 0
            for actClass in range(nb_classes):
                nb_prot = nb_ppc[actClass]
                mean = np.mean(train_set[train_lab == self.classes_[actClass], :], 0)
                self.w_[pos:pos + nb_prot] = mean + (random_state.rand(nb_prot, nb_features) * 2 - 1)
                self.c_w_[pos:pos + nb_prot] = self.classes_[actClass]
                pos += nb_prot
        else:
            X = validation.check_array(self.initial_prototypes)
            self.w_ = X[:,:-1]
            self.c_w_ = X[:,-1]
            if self.w_.shape != (np.sum(nb_ppc), nb_features):
                raise ValueError("the initial prototypes have wrong shape\n"
                                 "found=(%d,%d)\n"
                                 "expected=(%d,%d)" % (self.w_.shape[0], self.w_.shape[1], nb_ppc.sum(), nb_features))
            if set(self.c_w_) != set(self.classes_):
                raise ValueError("prototype labels and test data classes do not match\n"
                                 "classes={}\n"
                                 "prototype labels={}\n".format(self.classes_, self.c_w_))
        return train_set, train_lab, random_state

    def _optimize(self, X, y, random_state):
        label_equals_prototype = y[np.newaxis].T == self.c_w_
        res = minimize(
            fun=lambda x: self.optfun(variables=x, training_data=X, label_equals_prototype=label_equals_prototype),
            jac=lambda x: self.optgrad(variables=x, training_data=X, label_equals_prototype=label_equals_prototype,
                                       random_state=random_state),
            x0=self.w_, options={'disp': self.display, 'gtol': self.gtol, 'maxiter': self.max_iter})
        self.w_ = res.x.reshape(self.w_.shape)
        return res.nit

    def fit(self, X, y):
        X, y, random_state = self._validate_train_parms(X, y)
        self.n_iter_ = self._optimize(X, y, random_state)
        return self

    def _compute_distance(self, X, w=None):
        if w is None:
            w = self.w_
        return cdist(X, w, 'euclidean')

    def predict(self, X):
        check_is_fitted(self, ['w_', 'c_w_'])
        X = validation.check_array(X)
        if X.shape[1] != self.w_.shape[1]:
            raise ValueError("X has wrong number of features\n"
                             "found=%d\n"
                             "expected=%d" % (self.w_.shape[1], X.shape[1]))
        dist = self._compute_distance(X)
        return (self.c_w_[dist.argmin(1)])
