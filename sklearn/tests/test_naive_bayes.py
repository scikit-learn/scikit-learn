
import pickle
from io import BytesIO
import numpy as np
import scipy.sparse
import pytest

from sklearn.datasets import load_digits, load_iris

from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

from sklearn.utils._testing import assert_almost_equal
from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils._testing import assert_raises
from sklearn.utils._testing import assert_raise_message
from sklearn.utils._testing import assert_warns
from sklearn.utils._testing import assert_no_warnings

from sklearn.naive_bayes import GaussianNB, BernoulliNB
from sklearn.naive_bayes import MultinomialNB, ComplementNB
from sklearn.naive_bayes import CategoricalNB
from sklearn.naive_bayes import BaseNB, BaseDiscreteNB


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


def test_gnb():
    # Gaussian Naive Bayes classification.
    # This checks that GaussianNB implements fit and predict and returns
    # correct values for a simple toy dataset.

    clf = GaussianNB()
    y_pred = clf.fit(X, y).predict(X)
    assert_array_equal(y_pred, y)

    y_pred_proba = clf.predict_proba(X)
    y_pred_log_proba = clf.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)

    # Test whether label mismatch between target y and classes raises
    # an Error
    # FIXME Remove this test once the more general partial_fit tests are merged
    assert_raises(ValueError, GaussianNB().partial_fit, X, y, classes=[0, 1])


def test_gnb_prior():
    # Test whether class priors are properly set.
    clf = GaussianNB().fit(X, y)
    assert_array_almost_equal(np.array([3, 3]) / 6.0,
                              clf.class_prior_, 8)
    clf.fit(X1, y1)
    # Check that the class priors sum to 1
    assert_array_almost_equal(clf.class_prior_.sum(), 1)


def test_gnb_sample_weight():
    """Test whether sample weights are properly used in GNB. """
    # Sample weights all being 1 should not change results
    sw = np.ones(6)
    clf = GaussianNB().fit(X, y)
    clf_sw = GaussianNB().fit(X, y, sw)

    assert_array_almost_equal(clf.theta_, clf_sw.theta_)
    assert_array_almost_equal(clf.sigma_, clf_sw.sigma_)

    # Fitting twice with half sample-weights should result
    # in same result as fitting once with full weights
    sw = rng.rand(y.shape[0])
    clf1 = GaussianNB().fit(X, y, sample_weight=sw)
    clf2 = GaussianNB().partial_fit(X, y, classes=[1, 2], sample_weight=sw / 2)
    clf2.partial_fit(X, y, sample_weight=sw / 2)

    assert_array_almost_equal(clf1.theta_, clf2.theta_)
    assert_array_almost_equal(clf1.sigma_, clf2.sigma_)

    # Check that duplicate entries and correspondingly increased sample
    # weights yield the same result
    ind = rng.randint(0, X.shape[0], 20)
    sample_weight = np.bincount(ind, minlength=X.shape[0])

    clf_dupl = GaussianNB().fit(X[ind], y[ind])
    clf_sw = GaussianNB().fit(X, y, sample_weight)

    assert_array_almost_equal(clf_dupl.theta_, clf_sw.theta_)
    assert_array_almost_equal(clf_dupl.sigma_, clf_sw.sigma_)


def test_gnb_neg_priors():
    """Test whether an error is raised in case of negative priors"""
    clf = GaussianNB(priors=np.array([-1., 2.]))
    assert_raises(ValueError, clf.fit, X, y)


def test_gnb_priors():
    """Test whether the class prior override is properly used"""
    clf = GaussianNB(priors=np.array([0.3, 0.7])).fit(X, y)
    assert_array_almost_equal(clf.predict_proba([[-0.1, -0.1]]),
                              np.array([[0.825303662161683,
                                         0.174696337838317]]), 8)
    assert_array_almost_equal(clf.class_prior_, np.array([0.3, 0.7]))


def test_gnb_priors_sum_isclose():
    # test whether the class prior sum is properly tested"""
    X = np.array([[-1, -1], [-2, -1], [-3, -2], [-4, -5], [-5, -4],
                  [1, 1], [2, 1], [3, 2], [4, 4], [5, 5]])
    priors = np.array([0.08, 0.14, 0.03, 0.16, 0.11, 0.16, 0.07, 0.14,
                       0.11, 0.0])
    Y = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    clf = GaussianNB(priors)
    # smoke test for issue #9633
    clf.fit(X, Y)


def test_gnb_wrong_nb_priors():
    """ Test whether an error is raised if the number of prior is different
    from the number of class"""
    clf = GaussianNB(priors=np.array([.25, .25, .25, .25]))
    assert_raises(ValueError, clf.fit, X, y)


def test_gnb_prior_greater_one():
    """Test if an error is raised if the sum of prior greater than one"""
    clf = GaussianNB(priors=np.array([2., 1.]))
    assert_raises(ValueError, clf.fit, X, y)


def test_gnb_prior_large_bias():
    """Test if good prediction when class prior favor largely one class"""
    clf = GaussianNB(priors=np.array([0.01, 0.99]))
    clf.fit(X, y)
    assert clf.predict([[-0.1, -0.1]]) == np.array([2])


def test_gnb_check_update_with_no_data():
    """ Test when the partial fit is called without any data"""
    # Create an empty array
    prev_points = 100
    mean = 0.
    var = 1.
    x_empty = np.empty((0, X.shape[1]))
    tmean, tvar = GaussianNB._update_mean_variance(prev_points, mean,
                                                   var, x_empty)
    assert tmean == mean
    assert tvar == var


def test_gnb_pfit_wrong_nb_features():
    """Test whether an error is raised when the number of feature changes
    between two partial fit"""
    clf = GaussianNB()
    # Fit for the first time the GNB
    clf.fit(X, y)
    # Partial fit a second time with an incoherent X
    assert_raises(ValueError, clf.partial_fit, np.hstack((X, X)), y)


def test_gnb_partial_fit():
    clf = GaussianNB().fit(X, y)
    clf_pf = GaussianNB().partial_fit(X, y, np.unique(y))
    assert_array_almost_equal(clf.theta_, clf_pf.theta_)
    assert_array_almost_equal(clf.sigma_, clf_pf.sigma_)
    assert_array_almost_equal(clf.class_prior_, clf_pf.class_prior_)

    clf_pf2 = GaussianNB().partial_fit(X[0::2, :], y[0::2], np.unique(y))
    clf_pf2.partial_fit(X[1::2], y[1::2])
    assert_array_almost_equal(clf.theta_, clf_pf2.theta_)
    assert_array_almost_equal(clf.sigma_, clf_pf2.sigma_)
    assert_array_almost_equal(clf.class_prior_, clf_pf2.class_prior_)


def test_gnb_naive_bayes_scale_invariance():
    # Scaling the data should not change the prediction results
    iris = load_iris()
    X, y = iris.data, iris.target
    labels = [GaussianNB().fit(f * X, y).predict(f * X)
              for f in [1E-10, 1, 1E10]]
    assert_array_equal(labels[0], labels[1])
    assert_array_equal(labels[1], labels[2])


@pytest.mark.parametrize("cls", [MultinomialNB, BernoulliNB, CategoricalNB])
def test_discretenb_prior(cls):
    # Test whether class priors are properly set.
    clf = cls().fit(X2, y2)
    assert_array_almost_equal(np.log(np.array([2, 2, 2]) / 6.0),
                              clf.class_log_prior_, 8)


@pytest.mark.parametrize("cls", [MultinomialNB, BernoulliNB, CategoricalNB])
def test_discretenb_partial_fit(cls):
    clf1 = cls()
    clf1.fit([[0, 1], [1, 0], [1, 1]], [0, 1, 1])

    clf2 = cls()
    clf2.partial_fit([[0, 1], [1, 0], [1, 1]], [0, 1, 1], classes=[0, 1])
    assert_array_equal(clf1.class_count_, clf2.class_count_)
    if cls is CategoricalNB:
        for i in range(len(clf1.category_count_)):
            assert_array_equal(clf1.category_count_[i],
                               clf2.category_count_[i])
    else:
        assert_array_equal(clf1.feature_count_, clf2.feature_count_)

    clf3 = cls()
    # all categories have to appear in the first partial fit
    clf3.partial_fit([[0, 1]], [0], classes=[0, 1])
    clf3.partial_fit([[1, 0]], [1])
    clf3.partial_fit([[1, 1]], [1])
    assert_array_equal(clf1.class_count_, clf3.class_count_)
    if cls is CategoricalNB:
        # the categories for each feature of CategoricalNB are mapped to an
        # index chronologically with each call of partial fit and therefore
        # the category_count matrices cannot be compared for equality
        for i in range(len(clf1.category_count_)):
            assert_array_equal(clf1.category_count_[i].shape,
                               clf3.category_count_[i].shape)
            assert_array_equal(np.sum(clf1.category_count_[i], axis=1),
                               np.sum(clf3.category_count_[i], axis=1))

        # assert category 0 occurs 1x in the first class and 0x in the 2nd
        # class
        assert_array_equal(clf1.category_count_[0][0], np.array([1, 0]))
        # assert category 1 occurs 0x in the first class and 2x in the 2nd
        # class
        assert_array_equal(clf1.category_count_[0][1], np.array([0, 2]))

        # assert category 0 occurs 0x in the first class and 1x in the 2nd
        # class
        assert_array_equal(clf1.category_count_[1][0], np.array([0, 1]))
        # assert category 1 occurs 1x in the first class and 1x in the 2nd
        # class
        assert_array_equal(clf1.category_count_[1][1], np.array([1, 1]))
    else:
        assert_array_equal(clf1.feature_count_, clf3.feature_count_)


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, GaussianNB,
                                 CategoricalNB])
def test_discretenb_pickle(cls):
    # Test picklability of discrete naive Bayes classifiers

    clf = cls().fit(X2, y2)
    y_pred = clf.predict(X2)

    store = BytesIO()
    pickle.dump(clf, store)
    clf = pickle.load(BytesIO(store.getvalue()))

    assert_array_equal(y_pred, clf.predict(X2))

    # Test pickling of estimator trained with partial_fit
    clf2 = cls().partial_fit(X2[:3], y2[:3], classes=np.unique(y2))
    clf2.partial_fit(X2[3:], y2[3:])
    store = BytesIO()
    pickle.dump(clf2, store)
    clf2 = pickle.load(BytesIO(store.getvalue()))
    assert_array_equal(y_pred, clf2.predict(X2))


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, GaussianNB,
                                 CategoricalNB])
def test_discretenb_input_check_fit(cls):
    # Test input checks for the fit method

    # check shape consistency for number of samples at fit time
    assert_raises(ValueError, cls().fit, X2, y2[:-1])

    # check shape consistency for number of input features at predict time
    clf = cls().fit(X2, y2)
    assert_raises(ValueError, clf.predict, X2[:, :-1])


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, CategoricalNB])
def test_discretenb_input_check_partial_fit(cls):
    # check shape consistency
    assert_raises(ValueError, cls().partial_fit, X2, y2[:-1],
                  classes=np.unique(y2))

    # classes is required for first call to partial fit
    assert_raises(ValueError, cls().partial_fit, X2, y2)

    # check consistency of consecutive classes values
    clf = cls()
    clf.partial_fit(X2, y2, classes=np.unique(y2))
    assert_raises(ValueError, clf.partial_fit, X2, y2,
                  classes=np.arange(42))

    # check consistency of input shape for partial_fit
    assert_raises(ValueError, clf.partial_fit, X2[:, :-1], y2)

    # check consistency of input shape for predict
    assert_raises(ValueError, clf.predict, X2[:, :-1])


def test_discretenb_predict_proba():
    # Test discrete NB classes' probability scores

    # The 100s below distinguish Bernoulli from multinomial.
    # FIXME: write a test to show this.
    X_bernoulli = [[1, 100, 0], [0, 1, 0], [0, 100, 1]]
    X_multinomial = [[0, 1], [1, 3], [4, 0]]

    # test binary case (1-d output)
    y = [0, 0, 2]  # 2 is regression test for binary case, 02e673
    for cls, X in zip([BernoulliNB, MultinomialNB],
                      [X_bernoulli, X_multinomial]):
        clf = cls().fit(X, y)
        assert clf.predict(X[-1:]) == 2
        assert clf.predict_proba([X[0]]).shape == (1, 2)
        assert_array_almost_equal(clf.predict_proba(X[:2]).sum(axis=1),
                                  np.array([1., 1.]), 6)

    # test multiclass case (2-d output, must sum to one)
    y = [0, 1, 2]
    for cls, X in zip([BernoulliNB, MultinomialNB],
                      [X_bernoulli, X_multinomial]):
        clf = cls().fit(X, y)
        assert clf.predict_proba(X[0:1]).shape == (1, 3)
        assert clf.predict_proba(X[:2]).shape == (2, 3)
        assert_almost_equal(np.sum(clf.predict_proba([X[1]])), 1)
        assert_almost_equal(np.sum(clf.predict_proba([X[-1]])), 1)
        assert_almost_equal(np.sum(np.exp(clf.class_log_prior_)), 1)
        assert_almost_equal(np.sum(np.exp(clf.intercept_)), 1)


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, CategoricalNB])
def test_discretenb_uniform_prior(cls):
    # Test whether discrete NB classes fit a uniform prior
    # when fit_prior=False and class_prior=None

    clf = cls()
    clf.set_params(fit_prior=False)
    clf.fit([[0], [0], [1]], [0, 0, 1])
    prior = np.exp(clf.class_log_prior_)
    assert_array_almost_equal(prior, np.array([.5, .5]))


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, CategoricalNB])
def test_discretenb_provide_prior(cls):
    # Test whether discrete NB classes use provided prior

    clf = cls(class_prior=[0.5, 0.5])
    clf.fit([[0], [0], [1]], [0, 0, 1])
    prior = np.exp(clf.class_log_prior_)
    assert_array_almost_equal(prior, np.array([.5, .5]))

    # Inconsistent number of classes with prior
    assert_raises(ValueError, clf.fit, [[0], [1], [2]], [0, 1, 2])
    assert_raises(ValueError, clf.partial_fit, [[0], [1]], [0, 1],
                  classes=[0, 1, 1])


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, CategoricalNB])
def test_discretenb_provide_prior_with_partial_fit(cls):
    # Test whether discrete NB classes use provided prior
    # when using partial_fit

    iris = load_iris()
    iris_data1, iris_data2, iris_target1, iris_target2 = train_test_split(
        iris.data, iris.target, test_size=0.4, random_state=415)

    for prior in [None, [0.3, 0.3, 0.4]]:
        clf_full = cls(class_prior=prior)
        clf_full.fit(iris.data, iris.target)
        clf_partial = cls(class_prior=prior)
        clf_partial.partial_fit(iris_data1, iris_target1,
                                classes=[0, 1, 2])
        clf_partial.partial_fit(iris_data2, iris_target2)
        assert_array_almost_equal(clf_full.class_log_prior_,
                                  clf_partial.class_log_prior_)


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB, CategoricalNB])
def test_discretenb_sample_weight_multiclass(cls):
    # check shape consistency for number of samples at fit time
    X = [
        [0, 0, 1],
        [0, 1, 1],
        [0, 1, 1],
        [1, 0, 0],
    ]
    y = [0, 0, 1, 2]
    sample_weight = np.array([1, 1, 2, 2], dtype=np.float64)
    sample_weight /= sample_weight.sum()
    clf = cls().fit(X, y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X), [0, 1, 1, 2])

    # Check sample weight using the partial_fit method
    clf = cls()
    clf.partial_fit(X[:2], y[:2], classes=[0, 1, 2],
                    sample_weight=sample_weight[:2])
    clf.partial_fit(X[2:3], y[2:3], sample_weight=sample_weight[2:3])
    clf.partial_fit(X[3:], y[3:], sample_weight=sample_weight[3:])
    assert_array_equal(clf.predict(X), [0, 1, 1, 2])


@pytest.mark.parametrize('cls', [BernoulliNB, MultinomialNB])
def test_discretenb_coef_intercept_shape(cls):
    # coef_ and intercept_ should have shapes as in other linear models.
    # Non-regression test for issue #2127.
    X = [[1, 0, 0], [1, 1, 1]]
    y = [1, 2]  # binary classification
    clf = cls()

    clf.fit(X, y)
    assert clf.coef_.shape == (1, 3)
    assert clf.intercept_.shape == (1,)


@pytest.mark.parametrize('kind', ('dense', 'sparse'))
def test_mnnb(kind):
    # Test Multinomial Naive Bayes classification.
    # This checks that MultinomialNB implements fit and predict and returns
    # correct values for a simple toy dataset.

    if kind == 'dense':
        X = X2
    elif kind == 'sparse':
        X = scipy.sparse.csr_matrix(X2)

    # Check the ability to predict the learning set.
    clf = MultinomialNB()
    assert_raises(ValueError, clf.fit, -X, y2)
    y_pred = clf.fit(X, y2).predict(X)

    assert_array_equal(y_pred, y2)

    # Verify that np.log(clf.predict_proba(X)) gives the same results as
    # clf.predict_log_proba(X)
    y_pred_proba = clf.predict_proba(X)
    y_pred_log_proba = clf.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)

    # Check that incremental fitting yields the same results
    clf2 = MultinomialNB()
    clf2.partial_fit(X[:2], y2[:2], classes=np.unique(y2))
    clf2.partial_fit(X[2:5], y2[2:5])
    clf2.partial_fit(X[5:], y2[5:])

    y_pred2 = clf2.predict(X)
    assert_array_equal(y_pred2, y2)

    y_pred_proba2 = clf2.predict_proba(X)
    y_pred_log_proba2 = clf2.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba2), y_pred_log_proba2, 8)
    assert_array_almost_equal(y_pred_proba2, y_pred_proba)
    assert_array_almost_equal(y_pred_log_proba2, y_pred_log_proba)

    # Partial fit on the whole data at once should be the same as fit too
    clf3 = MultinomialNB()
    clf3.partial_fit(X, y2, classes=np.unique(y2))

    y_pred3 = clf3.predict(X)
    assert_array_equal(y_pred3, y2)
    y_pred_proba3 = clf3.predict_proba(X)
    y_pred_log_proba3 = clf3.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba3), y_pred_log_proba3, 8)
    assert_array_almost_equal(y_pred_proba3, y_pred_proba)
    assert_array_almost_equal(y_pred_log_proba3, y_pred_log_proba)


def test_mnb_prior_unobserved_targets():
    # test smoothing of prior for yet unobserved targets

    # Create toy training data
    X = np.array([[0, 1], [1, 0]])
    y = np.array([0, 1])

    clf = MultinomialNB()

    assert_no_warnings(
        clf.partial_fit, X, y, classes=[0, 1, 2]
    )

    assert clf.predict([[0, 1]]) == 0
    assert clf.predict([[1, 0]]) == 1
    assert clf.predict([[1, 1]]) == 0

    # add a training example with previously unobserved class
    assert_no_warnings(
        clf.partial_fit, [[1, 1]], [2]
    )

    assert clf.predict([[0, 1]]) == 0
    assert clf.predict([[1, 0]]) == 1
    assert clf.predict([[1, 1]]) == 2


def test_mnb_sample_weight():
    clf = MultinomialNB()
    clf.fit([[1, 2], [1, 2], [1, 0]],
            [0, 0, 1],
            sample_weight=[1, 1, 4])
    assert_array_equal(clf.predict([[1, 0]]), [1])
    positive_prior = np.exp(clf.intercept_[0])
    assert_array_almost_equal([1 - positive_prior, positive_prior],
                              [1 / 3., 2 / 3.])


def test_bnb():
    # Tests that BernoulliNB when alpha=1.0 gives the same values as
    # those given for the toy example in Manning, Raghavan, and
    # Schuetze's "Introduction to Information Retrieval" book:
    # https://nlp.stanford.edu/IR-book/html/htmledition/the-bernoulli-model-1.html

    # Training data points are:
    # Chinese Beijing Chinese (class: China)
    # Chinese Chinese Shanghai (class: China)
    # Chinese Macao (class: China)
    # Tokyo Japan Chinese (class: Japan)

    # Features are Beijing, Chinese, Japan, Macao, Shanghai, and Tokyo
    X = np.array([[1, 1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 1, 0],
                  [0, 1, 0, 1, 0, 0],
                  [0, 1, 1, 0, 0, 1]])

    # Classes are China (0), Japan (1)
    Y = np.array([0, 0, 0, 1])

    # Fit BernoulliBN w/ alpha = 1.0
    clf = BernoulliNB(alpha=1.0)
    clf.fit(X, Y)

    # Check the class prior is correct
    class_prior = np.array([0.75, 0.25])
    assert_array_almost_equal(np.exp(clf.class_log_prior_), class_prior)

    # Check the feature probabilities are correct
    feature_prob = np.array([[0.4, 0.8, 0.2, 0.4, 0.4, 0.2],
                             [1 / 3.0, 2 / 3.0, 2 / 3.0, 1 / 3.0, 1 / 3.0,
                              2 / 3.0]])
    assert_array_almost_equal(np.exp(clf.feature_log_prob_), feature_prob)

    # Testing data point is:
    # Chinese Chinese Chinese Tokyo Japan
    X_test = np.array([[0, 1, 1, 0, 0, 1]])

    # Check the predictive probabilities are correct
    unnorm_predict_proba = np.array([[0.005183999999999999,
                                      0.02194787379972565]])
    predict_proba = unnorm_predict_proba / np.sum(unnorm_predict_proba)
    assert_array_almost_equal(clf.predict_proba(X_test), predict_proba)


def test_bnb_feature_log_prob():
    # Test for issue #4268.
    # Tests that the feature log prob value computed by BernoulliNB when
    # alpha=1.0 is equal to the expression given in Manning, Raghavan,
    # and Schuetze's "Introduction to Information Retrieval" book:
    # http://nlp.stanford.edu/IR-book/html/htmledition/the-bernoulli-model-1.html

    X = np.array([[0, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 1], [0, 1, 0]])
    Y = np.array([0, 0, 1, 2, 2])

    # Fit Bernoulli NB w/ alpha = 1.0
    clf = BernoulliNB(alpha=1.0)
    clf.fit(X, Y)

    # Manually form the (log) numerator and denominator that
    # constitute P(feature presence | class)
    num = np.log(clf.feature_count_ + 1.0)
    denom = np.tile(np.log(clf.class_count_ + 2.0), (X.shape[1], 1)).T

    # Check manual estimate matches
    assert_array_almost_equal(clf.feature_log_prob_, (num - denom))


def test_cnb():
    # Tests ComplementNB when alpha=1.0 for the toy example in Manning,
    # Raghavan, and Schuetze's "Introduction to Information Retrieval" book:
    # https://nlp.stanford.edu/IR-book/html/htmledition/the-bernoulli-model-1.html

    # Training data points are:
    # Chinese Beijing Chinese (class: China)
    # Chinese Chinese Shanghai (class: China)
    # Chinese Macao (class: China)
    # Tokyo Japan Chinese (class: Japan)

    # Features are Beijing, Chinese, Japan, Macao, Shanghai, and Tokyo.
    X = np.array([[1, 1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 1, 0],
                  [0, 1, 0, 1, 0, 0],
                  [0, 1, 1, 0, 0, 1]])

    # Classes are China (0), Japan (1).
    Y = np.array([0, 0, 0, 1])

    # Check that weights are correct. See steps 4-6 in Table 4 of
    # Rennie et al. (2003).
    theta = np.array([
        [
            (0 + 1) / (3 + 6),
            (1 + 1) / (3 + 6),
            (1 + 1) / (3 + 6),
            (0 + 1) / (3 + 6),
            (0 + 1) / (3 + 6),
            (1 + 1) / (3 + 6)
        ],
        [
            (1 + 1) / (6 + 6),
            (3 + 1) / (6 + 6),
            (0 + 1) / (6 + 6),
            (1 + 1) / (6 + 6),
            (1 + 1) / (6 + 6),
            (0 + 1) / (6 + 6)
        ]])

    weights = np.zeros(theta.shape)
    normed_weights = np.zeros(theta.shape)
    for i in range(2):
        weights[i] = -np.log(theta[i])
        normed_weights[i] = weights[i] / weights[i].sum()

    # Verify inputs are nonnegative.
    clf = ComplementNB(alpha=1.0)
    assert_raises(ValueError, clf.fit, -X, Y)

    clf.fit(X, Y)

    # Check that counts/weights are correct.
    feature_count = np.array([[1, 3, 0, 1, 1, 0], [0, 1, 1, 0, 0, 1]])
    assert_array_equal(clf.feature_count_, feature_count)
    class_count = np.array([3, 1])
    assert_array_equal(clf.class_count_, class_count)
    feature_all = np.array([1, 4, 1, 1, 1, 1])
    assert_array_equal(clf.feature_all_, feature_all)
    assert_array_almost_equal(clf.feature_log_prob_, weights)

    clf = ComplementNB(alpha=1.0, norm=True)
    clf.fit(X, Y)
    assert_array_almost_equal(clf.feature_log_prob_, normed_weights)


def test_categoricalnb():
    # Check the ability to predict the training set.
    clf = CategoricalNB()
    y_pred = clf.fit(X2, y2).predict(X2)
    assert_array_equal(y_pred, y2)

    X3 = np.array([[1, 4], [2, 5]])
    y3 = np.array([1, 2])
    clf = CategoricalNB(alpha=1, fit_prior=False)

    clf.fit(X3, y3)

    # Check error is raised for X with negative entries
    X = np.array([[0, -1]])
    y = np.array([1])
    error_msg = "X must not contain negative values."
    assert_raise_message(ValueError, error_msg, clf.predict, X)
    assert_raise_message(ValueError, error_msg, clf.fit, X, y)

    # Test alpha
    X3_test = np.array([[2, 5]])
    # alpha=1 increases the count of all categories by one so the final
    # probability for each category is not 50/50 but 1/3 to 2/3
    bayes_numerator = np.array([[1/3*1/3, 2/3*2/3]])
    bayes_denominator = bayes_numerator.sum()
    assert_array_almost_equal(clf.predict_proba(X3_test),
                              bayes_numerator / bayes_denominator)

    # Assert category_count has counted all features
    assert len(clf.category_count_) == X3.shape[1]

    # Check sample_weight
    X = np.array([[0, 0], [0, 1], [0, 0], [1, 1]])
    y = np.array([1, 1, 2, 2])
    clf = CategoricalNB(alpha=1, fit_prior=False)
    clf.fit(X, y)
    assert_array_equal(clf.predict(np.array([[0, 0]])), np.array([1]))

    for factor in [1., 0.3, 5, 0.0001]:
        X = np.array([[0, 0], [0, 1], [0, 0], [1, 1]])
        y = np.array([1, 1, 2, 2])
        sample_weight = np.array([1, 1, 10, 0.1]) * factor
        clf = CategoricalNB(alpha=1, fit_prior=False)
        clf.fit(X, y, sample_weight=sample_weight)
        assert_array_equal(clf.predict(np.array([[0, 0]])), np.array([2]))


def test_alpha():
    # Setting alpha=0 should not output nan results when p(x_i|y_j)=0 is a case
    X = np.array([[1, 0], [1, 1]])
    y = np.array([0, 1])
    nb = BernoulliNB(alpha=0.)
    assert_warns(UserWarning, nb.partial_fit, X, y, classes=[0, 1])
    assert_warns(UserWarning, nb.fit, X, y)
    prob = np.array([[1, 0], [0, 1]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    nb = MultinomialNB(alpha=0.)
    assert_warns(UserWarning, nb.partial_fit, X, y, classes=[0, 1])
    assert_warns(UserWarning, nb.fit, X, y)
    prob = np.array([[2. / 3, 1. / 3], [0, 1]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    nb = CategoricalNB(alpha=0.)
    assert_warns(UserWarning, nb.fit, X, y)
    prob = np.array([[1., 0.], [0., 1.]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    # Test sparse X
    X = scipy.sparse.csr_matrix(X)
    nb = BernoulliNB(alpha=0.)
    assert_warns(UserWarning, nb.fit, X, y)
    prob = np.array([[1, 0], [0, 1]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    nb = MultinomialNB(alpha=0.)
    assert_warns(UserWarning, nb.fit, X, y)
    prob = np.array([[2. / 3, 1. / 3], [0, 1]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    # Test for alpha < 0
    X = np.array([[1, 0], [1, 1]])
    y = np.array([0, 1])
    expected_msg = ('Smoothing parameter alpha = -1.0e-01. '
                    'alpha should be > 0.')
    b_nb = BernoulliNB(alpha=-0.1)
    m_nb = MultinomialNB(alpha=-0.1)
    c_nb = CategoricalNB(alpha=-0.1)
    assert_raise_message(ValueError, expected_msg, b_nb.fit, X, y)
    assert_raise_message(ValueError, expected_msg, m_nb.fit, X, y)
    assert_raise_message(ValueError, expected_msg, c_nb.fit, X, y)

    b_nb = BernoulliNB(alpha=-0.1)
    m_nb = MultinomialNB(alpha=-0.1)
    assert_raise_message(ValueError, expected_msg, b_nb.partial_fit,
                         X, y, classes=[0, 1])
    assert_raise_message(ValueError, expected_msg, m_nb.partial_fit,
                         X, y, classes=[0, 1])


def test_alpha_vector():
    X = np.array([[1, 0], [1, 1]])
    y = np.array([0, 1])

    # Setting alpha=np.array with same length
    # as number of features should be fine
    alpha = np.array([1, 2])
    nb = MultinomialNB(alpha=alpha)
    nb.partial_fit(X, y, classes=[0, 1])

    # Test feature probabilities uses pseudo-counts (alpha)
    feature_prob = np.array([[1 / 2, 1 / 2], [2 / 5, 3 / 5]])
    assert_array_almost_equal(nb.feature_log_prob_, np.log(feature_prob))

    # Test predictions
    prob = np.array([[5 / 9, 4 / 9], [25 / 49, 24 / 49]])
    assert_array_almost_equal(nb.predict_proba(X), prob)

    # Test alpha non-negative
    alpha = np.array([1., -0.1])
    expected_msg = ('Smoothing parameter alpha = -1.0e-01. '
                    'alpha should be > 0.')
    m_nb = MultinomialNB(alpha=alpha)
    assert_raise_message(ValueError, expected_msg, m_nb.fit, X, y)

    # Test that too small pseudo-counts are replaced
    ALPHA_MIN = 1e-10
    alpha = np.array([ALPHA_MIN / 2, 0.5])
    m_nb = MultinomialNB(alpha=alpha)
    m_nb.partial_fit(X, y, classes=[0, 1])
    assert_array_almost_equal(m_nb._check_alpha(),
                              [ALPHA_MIN, 0.5],
                              decimal=12)

    # Test correct dimensions
    alpha = np.array([1., 2., 3.])
    m_nb = MultinomialNB(alpha=alpha)
    expected_msg = ('alpha should be a scalar or a numpy array '
                    'with shape [n_features]')
    assert_raise_message(ValueError, expected_msg, m_nb.fit, X, y)


def test_check_accuracy_on_digits():
    # Non regression test to make sure that any further refactoring / optim
    # of the NB models do not harm the performance on a slightly non-linearly
    # separable dataset
    X, y = load_digits(return_X_y=True)
    binary_3v8 = np.logical_or(y == 3, y == 8)
    X_3v8, y_3v8 = X[binary_3v8], y[binary_3v8]

    # Multinomial NB
    scores = cross_val_score(MultinomialNB(alpha=10), X, y, cv=10)
    assert scores.mean() > 0.86

    scores = cross_val_score(MultinomialNB(alpha=10), X_3v8, y_3v8, cv=10)
    assert scores.mean() > 0.94

    # Bernoulli NB
    scores = cross_val_score(BernoulliNB(alpha=10), X > 4, y, cv=10)
    assert scores.mean() > 0.83

    scores = cross_val_score(BernoulliNB(alpha=10), X_3v8 > 4, y_3v8, cv=10)
    assert scores.mean() > 0.92

    # Gaussian NB
    scores = cross_val_score(GaussianNB(), X, y, cv=10)
    assert scores.mean() > 0.77

    scores = cross_val_score(GaussianNB(var_smoothing=0.1), X, y, cv=10)
    assert scores.mean() > 0.89

    scores = cross_val_score(GaussianNB(), X_3v8, y_3v8, cv=10)
    assert scores.mean() > 0.86


# TODO: remove in 0.24
def test_deprecations():

    class A(BaseNB, GaussianNB):
        pass

    class B(BaseDiscreteNB, CategoricalNB):
        pass

    with pytest.warns(FutureWarning, match="is deprecated in version 0.22"):
        A()

    with pytest.warns(FutureWarning, match="is deprecated in version 0.22"):
        B()
