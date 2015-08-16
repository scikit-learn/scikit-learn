import numpy as np
from scipy.sparse import csr_matrix

from sklearn.naive_bayes import MultinomialNB
from profile_support import profile

__author__ = 'bas'


# TODO: is there a better place/format for a test like this in the codebase?
def run():
    clf = MultinomialNB(fit_prior=False)
    clf_sp = MultinomialNB(sparse=True, fit_prior=False)
    X, y = dataset()

    y_pred = fit_set(X, clf, y)
    y_pred_sp = fit_set_sp(X, clf_sp, y)
    #
    print np.allclose(y_pred, y_pred_sp)


@profile
def fit_set(X, clf, y):
    clf.fit(X, y)
    y_pred = clf.predict(X)
    return y_pred


@profile
def fit_set_sp(X, clf, y):
    clf.fit(X, y)
    y_pred = clf.predict(X)
    return y_pred


def dataset():
    n_samples = 50
    n_features = 2 ** 20
    n_classes = 10
    n_informative = 7
    X = csr_matrix((n_samples, n_features))
    for i in xrange(n_samples):
        indexes = np.random.randint(0, n_features, n_informative)
        for j in indexes:
            X[i, j] = 1
    y = np.random.randint(0, n_classes, n_samples)
    return X, y


if __name__ == '__main__':
    run()
