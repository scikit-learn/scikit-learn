"""
Testing for Isolation Forest algorithm (sklearn.ensemble.iforest).
"""

# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import pytest

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_allclose

from sklearn.model_selection import ParameterGrid
from sklearn.ensemble import IsolationForest
from sklearn.ensemble.iforest import _average_path_length
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_boston, load_iris
from sklearn.utils import check_random_state
from sklearn.metrics import roc_auc_score

from scipy.sparse import csc_matrix, csr_matrix
from unittest.mock import Mock, patch

rng = check_random_state(0)

# load the iris dataset
# and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


@pytest.mark.filterwarnings('ignore:threshold_ attribute')
def test_iforest():
    """Check Isolation Forest for various parameter settings."""
    X_train = np.array([[0, 1], [1, 2]])
    X_test = np.array([[2, 1], [1, 1]])

    grid = ParameterGrid({"n_estimators": [3],
                          "max_samples": [0.5, 1.0, 3],
                          "bootstrap": [True, False]})

    with ignore_warnings():
        for params in grid:
            IsolationForest(random_state=rng,
                            **params).fit(X_train).predict(X_test)


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:threshold_ attribute')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_iforest_sparse():
    """Check IForest for various parameter settings on sparse input."""
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "bootstrap": [True, False]})

    for sparse_format in [csc_matrix, csr_matrix]:
        X_train_sparse = sparse_format(X_train)
        X_test_sparse = sparse_format(X_test)

        for params in grid:
            # Trained on sparse format
            sparse_classifier = IsolationForest(
                n_estimators=10, random_state=1, **params).fit(X_train_sparse)
            sparse_results = sparse_classifier.predict(X_test_sparse)

            # Trained on dense format
            dense_classifier = IsolationForest(
                n_estimators=10, random_state=1, **params).fit(X_train)
            dense_results = dense_classifier.predict(X_test)

            assert_array_equal(sparse_results, dense_results)


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:threshold_ attribute')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_iforest_error():
    """Test that it gives proper exception on deficient input."""
    X = iris.data

    # Test max_samples
    assert_raises(ValueError,
                  IsolationForest(max_samples=-1).fit, X)
    assert_raises(ValueError,
                  IsolationForest(max_samples=0.0).fit, X)
    assert_raises(ValueError,
                  IsolationForest(max_samples=2.0).fit, X)
    # The dataset has less than 256 samples, explicitly setting
    # max_samples > n_samples should result in a warning. If not set
    # explicitly there should be no warning
    assert_warns_message(UserWarning,
                         "max_samples will be set to n_samples for estimation",
                         IsolationForest(max_samples=1000).fit, X)
    # note that assert_no_warnings does not apply since it enables a
    # PendingDeprecationWarning triggered by scipy.sparse's use of
    # np.matrix. See issue #11251.
    with pytest.warns(None) as record:
        IsolationForest(max_samples='auto').fit(X)
    user_warnings = [each for each in record
                     if issubclass(each.category, UserWarning)]
    assert len(user_warnings) == 0
    with pytest.warns(None) as record:
        IsolationForest(max_samples=np.int64(2)).fit(X)
    user_warnings = [each for each in record
                     if issubclass(each.category, UserWarning)]
    assert len(user_warnings) == 0

    assert_raises(ValueError, IsolationForest(max_samples='foobar').fit, X)
    assert_raises(ValueError, IsolationForest(max_samples=1.5).fit, X)

    # test X_test n_features match X_train one:
    assert_raises(ValueError, IsolationForest().fit(X).predict, X[:, 1:])

    # test threshold_ attribute error when behaviour is not old:
    msg = "threshold_ attribute does not exist when behaviour != 'old'"
    assert_raises_regex(AttributeError, msg, getattr,
                        IsolationForest(behaviour='new'), 'threshold_')


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_recalculate_max_depth():
    """Check max_depth recalculation when max_samples is reset to n_samples"""
    X = iris.data
    clf = IsolationForest().fit(X)
    for est in clf.estimators_:
        assert_equal(est.max_depth, int(np.ceil(np.log2(X.shape[0]))))


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_max_samples_attribute():
    X = iris.data
    clf = IsolationForest().fit(X)
    assert_equal(clf.max_samples_, X.shape[0])

    clf = IsolationForest(max_samples=500)
    assert_warns_message(UserWarning,
                         "max_samples will be set to n_samples for estimation",
                         clf.fit, X)
    assert_equal(clf.max_samples_, X.shape[0])

    clf = IsolationForest(max_samples=0.4).fit(X)
    assert_equal(clf.max_samples_, 0.4*X.shape[0])


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:threshold_ attribute')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_iforest_parallel_regression():
    """Check parallel regression."""
    rng = check_random_state(0)

    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    ensemble = IsolationForest(n_jobs=3,
                               random_state=0).fit(X_train)

    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = IsolationForest(n_jobs=1,
                               random_state=0).fit(X_train)

    y3 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y3)


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_iforest_performance():
    """Test Isolation Forest performs well"""

    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2)
    X_train = np.r_[X + 2, X - 2]
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2))
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model
    clf = IsolationForest(max_samples=100, random_state=rng).fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = - clf.decision_function(X_test)

    # check that there is at most 6 errors (false positive or false negative)
    assert_greater(roc_auc_score(y_test, y_pred), 0.98)


@pytest.mark.parametrize("contamination", [0.25, "auto"])
@pytest.mark.filterwarnings("ignore:threshold_ attribute")
def test_iforest_works(contamination):
    # toy sample (the last two samples are outliers)
    X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [6, 3], [-4, 7]]

    # Test IsolationForest
    clf = IsolationForest(
        behaviour="new", random_state=rng, contamination=contamination
    )
    clf.fit(X)
    decision_func = -clf.decision_function(X)
    pred = clf.predict(X)
    # assert detect outliers:
    assert_greater(np.min(decision_func[-2:]), np.max(decision_func[:-2]))
    assert_array_equal(pred, 6 * [1] + 2 * [-1])


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_max_samples_consistency():
    # Make sure validated max_samples in iforest and BaseBagging are identical
    X = iris.data
    clf = IsolationForest().fit(X)
    assert_equal(clf.max_samples_, clf._max_samples)


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:threshold_ attribute')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_iforest_subsampled_features():
    # It tests non-regression for #5732 which failed at predict.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)
    clf = IsolationForest(max_features=0.8)
    clf.fit(X_train, y_train)
    clf.predict(X_test)


def test_iforest_average_path_length():
    # It tests non-regression for #8549 which used the wrong formula
    # for average path length, strictly for the integer case
    # Updated to check average path length when input is <= 2 (issue #11839)
    result_one = 2.0 * (np.log(4.0) + np.euler_gamma) - 2.0 * 4.0 / 5.0
    result_two = 2.0 * (np.log(998.0) + np.euler_gamma) - 2.0 * 998.0 / 999.0
    assert_allclose(_average_path_length([0]), [0.0])
    assert_allclose(_average_path_length([1]), [0.0])
    assert_allclose(_average_path_length([2]), [1.0])
    assert_allclose(_average_path_length([5]), [result_one])
    assert_allclose(_average_path_length([999]), [result_two])
    assert_allclose(
        _average_path_length(np.array([1, 2, 5, 999])),
        [0.0, 1.0, result_one, result_two],
    )
    # _average_path_length is increasing
    avg_path_length = _average_path_length(np.arange(5))
    assert_array_equal(avg_path_length, np.sort(avg_path_length))


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_score_samples():
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = IsolationForest(contamination=0.1).fit(X_train)
    clf2 = IsolationForest().fit(X_train)
    assert_array_equal(clf1.score_samples([[2., 2.]]),
                       clf1.decision_function([[2., 2.]]) + clf1.offset_)
    assert_array_equal(clf2.score_samples([[2., 2.]]),
                       clf2.decision_function([[2., 2.]]) + clf2.offset_)
    assert_array_equal(clf1.score_samples([[2., 2.]]),
                       clf2.score_samples([[2., 2.]]))


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_deprecation():
    X = [[0.0], [1.0]]
    clf = IsolationForest()

    assert_warns_message(FutureWarning,
                         'default contamination parameter 0.1 will change '
                         'in version 0.22 to "auto"',
                         clf.fit, X)

    assert_warns_message(FutureWarning,
                         'behaviour="old" is deprecated and will be removed '
                         'in version 0.22',
                         clf.fit, X)

    clf = IsolationForest().fit(X)
    assert_warns_message(DeprecationWarning,
                         "threshold_ attribute is deprecated in 0.20 and will"
                         " be removed in 0.22.",
                         getattr, clf, "threshold_")


@pytest.mark.filterwarnings('ignore:default contamination')
@pytest.mark.filterwarnings('ignore:behaviour="old"')
def test_behaviour_param():
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = IsolationForest(behaviour='old').fit(X_train)
    clf2 = IsolationForest(behaviour='new', contamination='auto').fit(X_train)
    assert_array_equal(clf1.decision_function([[2., 2.]]),
                       clf2.decision_function([[2., 2.]]))


# mock get_chunk_n_rows to actually test more than one chunk (here one
# chunk = 3 rows:
@patch(
    "sklearn.ensemble.iforest.get_chunk_n_rows",
    side_effect=Mock(**{"return_value": 3}),
)
@pytest.mark.parametrize(
    "contamination, n_predict_calls", [(0.25, 3), ("auto", 2)]
)
@pytest.mark.filterwarnings("ignore:threshold_ attribute")
def test_iforest_chunks_works1(
    mocked_get_chunk, contamination, n_predict_calls
):
    test_iforest_works(contamination)
    assert mocked_get_chunk.call_count == n_predict_calls


# idem with chunk_size = 5 rows
@patch(
    "sklearn.ensemble.iforest.get_chunk_n_rows",
    side_effect=Mock(**{"return_value": 10}),
)
@pytest.mark.parametrize(
    "contamination, n_predict_calls", [(0.25, 3), ("auto", 2)]
)
@pytest.mark.filterwarnings("ignore:threshold_ attribute")
def test_iforest_chunks_works2(
    mocked_get_chunk, contamination, n_predict_calls
):
    test_iforest_works(contamination)
    assert mocked_get_chunk.call_count == n_predict_calls
