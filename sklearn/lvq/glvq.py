import numpy as np
from scipy.spatial.distance import cdist
from sklearn.utils import validation

from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_is_fitted


class GlvqModel(BaseEstimator):
    def __init__(self, random_state=None):
        self.random_state = validation.check_random_state(random_state)
        self.prototypes_per_class = 1
        self.w = None
        self.c_w = None

    def initial_prototypes(self, w, c_w):
        self.w = validation.check_array(w, dtype='float')
        self.c_w = validation.column_or_1d(validation.check_array(c_w, ensure_2d=False))

    def fit(self, train_set, train_lab, nb_epochs=100, learning_rate_prototypes=None):
        # validate arguments
        train_set = validation.as_float_array(validation.check_array(train_set))
        train_lab = validation.column_or_1d(validation.check_array(train_lab, ensure_2d=False))
        validation.check_X_y(train_set, train_lab)
        if not isinstance(nb_epochs, int) or nb_epochs <= 0:
            raise ValueError("nb_epochs must be an int greater than 0")

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
            self.w = np.empty([np.sum(nb_ppc), nb_features], dtype='float')
            self.c_w = np.empty([np.sum(nb_ppc)])
            pos = 0
            for actClass in range(nb_classes):
                nb_prot = nb_ppc[actClass]
                mean = np.mean(train_set[train_lab == classes[actClass], :], 0)
                self.w[pos:pos+nb_prot] = mean + (self.random_state.rand(nb_prot,nb_features) * 2 - 1)
                self.c_w[pos:pos+nb_prot] = classes[actClass]
                pos += nb_prot
        else:
            self.w = validation.check_array(self.w, dtype='float64')
            self.c_w = validation.column_or_1d(validation.check_array(self.c_w, ensure_2d=False))
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

        # compute the vector of nb_epochs learning rates alpha for the prototype learning
        if learning_rate_prototypes is not None and len(learning_rate_prototypes) > 2:
            learning_rate_prototypes = validation.column_or_1d(
                validation.check_array(learning_rate_prototypes, ensure_2d=False, dtype='float'))
            if learning_rate_prototypes.shape[0] != nb_epochs:
                raise ValueError("The learning rate vector for the prototypes does not fit the nb_epochs\n"
                                 "nb_epochs=%d\n"
                                 "length of vector=%d" % (nb_epochs, learning_rate_prototypes.shape[0]))
            alphas = learning_rate_prototypes
        else:
            if learning_rate_prototypes is None or len(learning_rate_prototypes) is not 2:
                learning_rate_prototypes = [nb_features / 100, nb_features / 10000]
            alpha_start = learning_rate_prototypes[0]
            alpha_end = learning_rate_prototypes[1]
            alphas = np.array([alpha_start * pow(alpha_end / alpha_start, x / nb_epochs) for x in range(nb_epochs)])
        for epoch in range(nb_epochs):
            order = self.random_state.permutation(nb_samples)
            for j in range(nb_samples):
                index = order[j]
                xi = train_set[index]
                c_xi = train_lab[index]
                dist = cdist(xi.reshape(1, -1), self.w, 'euclidean').ravel()
                sortIdx = np.argsort(dist)
                count = 0
                while self.c_w[sortIdx[count]] != c_xi:
                    count += 1
                J = sortIdx[count]
                dJ = dist[sortIdx[count]]
                wJ = self.w[J]

                count = 0
                while self.c_w[sortIdx[count]] == c_xi:
                    count += 1
                K = sortIdx[count]
                dK = dist[sortIdx[count]]
                wK = self.w[K]
                norm_factor = pow((dJ + dK), 2)
                DJ = (xi - wJ)
                DK = (xi - wK)

                dwJ = (2 * dK / norm_factor) * 2 * DJ
                dwK = (2 * dJ / norm_factor) * 2 * DK

                self.w[J] += alphas[epoch] * dwJ
                self.w[K] -= alphas[epoch] * dwK

    def predict(self, X):
        check_is_fitted(self, ['w', 'c_w'])
        dist = cdist(X, self.w, 'euclidean')
        return (self.c_w[dist.argmin(1)])

    def __str__(self):
        return 'w:' + str(self.w) + '\n' \
               + 'c_w:' + str(self.c_w)
