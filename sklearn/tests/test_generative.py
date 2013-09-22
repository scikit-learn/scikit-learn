"""Tests for Generative Classification"""
import numpy as np
from sklearn.generative import GenerativeBayes
from sklearn.naive_bayes import GaussianNB
from sklearn.utils.testing import assert_array_almost_equal, assert_greater
from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_digits

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
y = np.array([1, 1, 1, 2, 2, 2])

# A bit more random tests
rng = np.random.RandomState(0)
X1 = rng.normal(size=(10, 3))
y1 = (rng.normal(size=(10)) > 0).astype(np.int)

# Data is 6 random integer points in a 100 dimensional space classified to
# three classes.
X2 = rng.randint(5, size=(6, 100))
y2 = np.array([1, 1, 2, 2, 3, 3])


MODEL_KWARGS = {'norm_approx': {},
                'gmm':{'n_components': 3,
                       'covariance_type': 'diag'},
                'kde':{'bandwidth': 4.0,
                       'kernel': 'gaussian',
                       'metric': 'euclidean'}}


def test_compare_generative_gnb():
    """Compare GenerativeBayes to GaussianNB"""
    # using norm_approx, the two should yield an identical model.
    clf1 = GenerativeBayes('norm_approx')
    clf2 = GaussianNB()

    p1 = clf1.fit(X, y).predict_proba(X)
    p2 = clf2.fit(X, y).predict_proba(X)

    assert_array_almost_equal(p1, p2)


def test_generative_model_prior():
    """Test whether class priors are properly set."""
    # class priors should sum to 1.
    for model, kwargs in MODEL_KWARGS.iteritems():
        clf = GenerativeBayes(model, **kwargs)

        clf.fit(X, y)
        assert_array_almost_equal(np.array([3, 3]) / 6.0,
                                  clf.class_prior_, 8)

        clf.fit(X1, y1)
        assert_array_almost_equal(clf.class_prior_.sum(), 1)

def test_check_accuracy_on_digits():
    digits = load_digits()
    X, y = digits.data, digits.target

    binary_3v8 = np.logical_or(digits.target == 3, digits.target == 8)
    X_3v8, y_3v8 = X[binary_3v8], y[binary_3v8]

    scores_cmp = {'kde': (0.98, 0.99),
                  'norm_approx': (0.89, 0.93),
                  'gmm': (0.92, 0.93)}

    for model, kwargs in MODEL_KWARGS.iteritems():
        scores = cross_val_score(GenerativeBayes(model, **kwargs),
                                 X, y, cv=4)
        assert_greater(scores.mean(), scores_cmp[model][0])

        scores = cross_val_score(GenerativeBayes(model, **kwargs),
                                 X_3v8, y_3v8, cv=4)
        assert_greater(scores.mean(), scores_cmp[model][1])


if __name__ == '__main__':
    import nose
    nose.runmodule()
