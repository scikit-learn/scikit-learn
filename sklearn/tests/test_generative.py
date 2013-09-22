"""Tests for Generative Classification"""
import numpy as np
from sklearn.generative import GenerativeBayes
from sklearn.naive_bayes import GaussianNB
from sklearn.utils.testing import\
    assert_array_almost_equal, assert_greater, assert_allclose
from sklearn.cross_validation import cross_val_score
from sklearn.datasets import load_digits, make_blobs

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


def test_generate_samples():
    # create a simple unbalanced dataset with 4 classes
    np.random.seed(0)
    X1, y1 = make_blobs(50, 2, centers=2)
    X2, y2 = make_blobs(100, 2, centers=2)
    X = np.vstack([X1, X2])
    y = np.concatenate([y1, y2 + 2])

    # test with norm_approx; other models have their sample() method
    # tested independently
    clf = GenerativeBayes('norm_approx')
    clf.fit(X, y)

    X_new, y_new = clf.sample(2000)

    for i in range(4):
        Xnew_i = X_new[y_new == i]
        X_i = X[y == i]

        # check the means
        assert_array_almost_equal(Xnew_i.mean(0),
                                  X_i.mean(0), decimal=1)

        # check the standard deviations
        assert_array_almost_equal(Xnew_i.std(0),
                                  X_i.std(0), decimal=1)

        # check the number of points
        assert_allclose(X_i.shape[0] * 1. / X.shape[0],
                        Xnew_i.shape[0] * 1. / X_new.shape[0],
                        rtol=0.1)


if __name__ == '__main__':
    import nose
    nose.runmodule()
