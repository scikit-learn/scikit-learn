# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import sparse

from sklearn.base import BaseEstimator
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import LeaveOneOut, train_test_split

from sklearn.utils._testing import (assert_array_almost_equal,
                                    assert_almost_equal,
                                    assert_array_equal,
                                    assert_raises, ignore_warnings)
from sklearn.utils.extmath import softmax
from sklearn.exceptions import NotFittedError
from sklearn.datasets import make_classification, make_blobs
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.naive_bayes import MultinomialNB
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.svm import LinearSVC
from sklearn.isotonic import IsotonicRegression
from sklearn.feature_extraction import DictVectorizer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import brier_score_loss
from sklearn.calibration import CalibratedClassifierCV, _CalibratedClassifier
from sklearn.calibration import _sigmoid_calibration, _SigmoidCalibration
from sklearn.calibration import calibration_curve


@pytest.fixture(scope="module")
def data():
    X, y = make_classification(
        n_samples=200, n_features=6, random_state=42
    )
    return X, y


@pytest.mark.parametrize('method', ['sigmoid', 'isotonic'])
@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration(data, method, ensemble):
    # Test calibration objects with isotonic and sigmoid
    n_samples = 100
    X, y = data
    sample_weight = np.random.RandomState(seed=42).uniform(size=y.size)

    X -= X.min()  # MultinomialNB only allows positive X

    # split train and test
    X_train, y_train, sw_train = \
        X[:n_samples], y[:n_samples], sample_weight[:n_samples]
    X_test, y_test = X[n_samples:], y[n_samples:]

    # Naive-Bayes
    clf = MultinomialNB().fit(X_train, y_train, sample_weight=sw_train)
    prob_pos_clf = clf.predict_proba(X_test)[:, 1]

    cal_clf = CalibratedClassifierCV(clf, cv=y.size + 1, ensemble=ensemble)
    assert_raises(ValueError, cal_clf.fit, X, y)

    # Naive Bayes with calibration
    for this_X_train, this_X_test in [(X_train, X_test),
                                      (sparse.csr_matrix(X_train),
                                       sparse.csr_matrix(X_test))]:
        cal_clf = CalibratedClassifierCV(
            clf, method=method, cv=5, ensemble=ensemble
        )
        # Note that this fit overwrites the fit on the entire training
        # set
        cal_clf.fit(this_X_train, y_train, sample_weight=sw_train)
        prob_pos_cal_clf = cal_clf.predict_proba(this_X_test)[:, 1]

        # Check that brier score has improved after calibration
        assert (brier_score_loss(y_test, prob_pos_clf) >
                brier_score_loss(y_test, prob_pos_cal_clf))

        # Check invariance against relabeling [0, 1] -> [1, 2]
        cal_clf.fit(this_X_train, y_train + 1, sample_weight=sw_train)
        prob_pos_cal_clf_relabeled = cal_clf.predict_proba(this_X_test)[:, 1]
        assert_array_almost_equal(prob_pos_cal_clf,
                                  prob_pos_cal_clf_relabeled)

        # Check invariance against relabeling [0, 1] -> [-1, 1]
        cal_clf.fit(this_X_train, 2 * y_train - 1, sample_weight=sw_train)
        prob_pos_cal_clf_relabeled = cal_clf.predict_proba(this_X_test)[:, 1]
        assert_array_almost_equal(prob_pos_cal_clf, prob_pos_cal_clf_relabeled)

        # Check invariance against relabeling [0, 1] -> [1, 0]
        cal_clf.fit(this_X_train, (y_train + 1) % 2, sample_weight=sw_train)
        prob_pos_cal_clf_relabeled = cal_clf.predict_proba(this_X_test)[:, 1]
        if method == "sigmoid":
            assert_array_almost_equal(prob_pos_cal_clf,
                                      1 - prob_pos_cal_clf_relabeled)
        else:
            # Isotonic calibration is not invariant against relabeling
            # but should improve in both cases
            assert (brier_score_loss(y_test, prob_pos_clf) >
                    brier_score_loss((y_test + 1) % 2,
                                     prob_pos_cal_clf_relabeled))


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_bad_method(data, ensemble):
    # Check only "isotonic" and "sigmoid" are accepted as methods
    X, y = data
    clf = LinearSVC()
    clf_invalid_method = CalibratedClassifierCV(
        clf, method="foo", ensemble=ensemble
    )
    with pytest.raises(ValueError):
        clf_invalid_method.fit(X, y)


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_regressor(data, ensemble):
    # `base-estimator` should provide either decision_function or
    # predict_proba (most regressors, for instance, should fail)
    X, y = data
    clf_base_regressor = \
        CalibratedClassifierCV(RandomForestRegressor(), ensemble=ensemble)
    with pytest.raises(RuntimeError):
        clf_base_regressor.fit(X, y)


def test_calibration_default_estimator(data):
    # Check base_estimator default is LinearSVC
    X, y = data
    calib_clf = CalibratedClassifierCV(cv=2)
    calib_clf.fit(X, y)

    base_est = calib_clf.calibrated_classifiers_[0].base_estimator
    assert isinstance(base_est, LinearSVC)


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_cv_splitter(data, ensemble):
    # Check when `cv` is a CV splitter
    X, y = data

    splits = 5
    kfold = KFold(n_splits=splits)
    calib_clf = CalibratedClassifierCV(cv=kfold, ensemble=ensemble)
    assert isinstance(calib_clf.cv, KFold)
    assert calib_clf.cv.n_splits == splits

    calib_clf.fit(X, y)
    expected_n_clf = splits if ensemble else 1
    assert len(calib_clf.calibrated_classifiers_) == expected_n_clf


@pytest.mark.parametrize('method', ['sigmoid', 'isotonic'])
@pytest.mark.parametrize('ensemble', [True, False])
def test_sample_weight(data, method, ensemble):
    n_samples = 100
    X, y = data

    sample_weight = np.random.RandomState(seed=42).uniform(size=len(y))
    X_train, y_train, sw_train = \
        X[:n_samples], y[:n_samples], sample_weight[:n_samples]
    X_test = X[n_samples:]

    base_estimator = LinearSVC(random_state=42)
    calibrated_clf = CalibratedClassifierCV(
        base_estimator, method=method, ensemble=ensemble
    )
    calibrated_clf.fit(X_train, y_train, sample_weight=sw_train)
    probs_with_sw = calibrated_clf.predict_proba(X_test)

    # As the weights are used for the calibration, they should still yield
    # different predictions
    calibrated_clf.fit(X_train, y_train)
    probs_without_sw = calibrated_clf.predict_proba(X_test)

    diff = np.linalg.norm(probs_with_sw - probs_without_sw)
    assert diff > 0.1


@pytest.mark.parametrize('method', ['sigmoid', 'isotonic'])
@pytest.mark.parametrize('ensemble', [True, False])
def test_parallel_execution(data, method, ensemble):
    """Test parallel calibration"""
    X, y = data
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

    base_estimator = LinearSVC(random_state=42)

    cal_clf_parallel = CalibratedClassifierCV(
        base_estimator, method=method, n_jobs=2, ensemble=ensemble
    )
    cal_clf_parallel.fit(X_train, y_train)
    probs_parallel = cal_clf_parallel.predict_proba(X_test)

    cal_clf_sequential = CalibratedClassifierCV(
        base_estimator, method=method, n_jobs=1, ensemble=ensemble
    )
    cal_clf_sequential.fit(X_train, y_train)
    probs_sequential = cal_clf_sequential.predict_proba(X_test)

    assert_allclose(probs_parallel, probs_sequential)


@pytest.mark.parametrize('method', ['sigmoid', 'isotonic'])
@pytest.mark.parametrize('ensemble', [True, False])
# increase the number of RNG seeds to assess the statistical stability of this
# test:
@pytest.mark.parametrize('seed', range(2))
def test_calibration_multiclass(method, ensemble, seed):

    def multiclass_brier(y_true, proba_pred, n_classes):
        Y_onehot = np.eye(n_classes)[y_true]
        return np.sum((Y_onehot - proba_pred) ** 2) / Y_onehot.shape[0]

    # Test calibration for multiclass with classifier that implements
    # only decision function.
    clf = LinearSVC(random_state=7)
    X, y = make_blobs(n_samples=500, n_features=100, random_state=seed,
                      centers=10, cluster_std=15.0)

    # Use an unbalanced dataset by collapsing 8 clusters into one class
    # to make the naive calibration based on a softmax more unlikely
    # to work.
    y[y > 2] = 2
    n_classes = np.unique(y).shape[0]
    X_train, y_train = X[::2], y[::2]
    X_test, y_test = X[1::2], y[1::2]

    clf.fit(X_train, y_train)

    cal_clf = CalibratedClassifierCV(
        clf, method=method, cv=5, ensemble=ensemble
    )
    cal_clf.fit(X_train, y_train)
    probas = cal_clf.predict_proba(X_test)
    # Check probabilities sum to 1
    assert_allclose(np.sum(probas, axis=1), np.ones(len(X_test)))

    # Check that the dataset is not too trivial, otherwise it's hard
    # to get interesting calibration data during the internal
    # cross-validation loop.
    assert 0.65 < clf.score(X_test, y_test) < 0.95

    # Check that the accuracy of the calibrated model is never degraded
    # too much compared to the original classifier.
    assert cal_clf.score(X_test, y_test) > 0.95 * clf.score(X_test, y_test)

    # Check that Brier loss of calibrated classifier is smaller than
    # loss obtained by naively turning OvR decision function to
    # probabilities via a softmax
    uncalibrated_brier = \
        multiclass_brier(y_test, softmax(clf.decision_function(X_test)),
                         n_classes=n_classes)
    calibrated_brier = multiclass_brier(y_test, probas,
                                        n_classes=n_classes)

    assert calibrated_brier < 1.1 * uncalibrated_brier

    # Test that calibration of a multiclass classifier decreases log-loss
    # for RandomForestClassifier
    clf = RandomForestClassifier(n_estimators=30, random_state=42)
    clf.fit(X_train, y_train)
    clf_probs = clf.predict_proba(X_test)
    uncalibrated_brier = multiclass_brier(y_test, clf_probs,
                                          n_classes=n_classes)

    cal_clf = CalibratedClassifierCV(
        clf, method=method, cv=5, ensemble=ensemble
    )
    cal_clf.fit(X_train, y_train)
    cal_clf_probs = cal_clf.predict_proba(X_test)
    calibrated_brier = multiclass_brier(y_test, cal_clf_probs,
                                        n_classes=n_classes)
    assert calibrated_brier < 1.1 * uncalibrated_brier


def test_calibration_zero_probability():
    # Test an edge case where _CalibratedClassifier avoids numerical errors
    # in the multiclass normalization step if all the calibrators output
    # are zero all at once for a given sample and instead fallback to uniform
    # probabilities.
    class ZeroCalibrator():
        # This function is called from _CalibratedClassifier.predict_proba.
        def predict(self, X):
            return np.zeros(X.shape[0])

    X, y = make_blobs(n_samples=50, n_features=10, random_state=7,
                      centers=10, cluster_std=15.0)
    clf = DummyClassifier().fit(X, y)
    calibrator = ZeroCalibrator()
    cal_clf = _CalibratedClassifier(
        base_estimator=clf, calibrators=[calibrator], classes=clf.classes_)

    probas = cal_clf.predict_proba(X)

    # Check that all probabilities are uniformly 1. / clf.n_classes_
    assert_allclose(probas, 1. / clf.n_classes_)


def test_calibration_prefit():
    """Test calibration for prefitted classifiers"""
    n_samples = 50
    X, y = make_classification(n_samples=3 * n_samples, n_features=6,
                               random_state=42)
    sample_weight = np.random.RandomState(seed=42).uniform(size=y.size)

    X -= X.min()  # MultinomialNB only allows positive X

    # split train and test
    X_train, y_train, sw_train = \
        X[:n_samples], y[:n_samples], sample_weight[:n_samples]
    X_calib, y_calib, sw_calib = \
        X[n_samples:2 * n_samples], y[n_samples:2 * n_samples], \
        sample_weight[n_samples:2 * n_samples]
    X_test, y_test = X[2 * n_samples:], y[2 * n_samples:]

    # Naive-Bayes
    clf = MultinomialNB()
    # Check error if clf not prefit
    unfit_clf = CalibratedClassifierCV(clf, cv="prefit")
    with pytest.raises(NotFittedError):
        unfit_clf.fit(X_calib, y_calib)

    clf.fit(X_train, y_train, sw_train)
    prob_pos_clf = clf.predict_proba(X_test)[:, 1]

    # Naive Bayes with calibration
    for this_X_calib, this_X_test in [(X_calib, X_test),
                                      (sparse.csr_matrix(X_calib),
                                       sparse.csr_matrix(X_test))]:
        for method in ['isotonic', 'sigmoid']:
            cal_clf = CalibratedClassifierCV(clf, method=method, cv="prefit")

            for sw in [sw_calib, None]:
                cal_clf.fit(this_X_calib, y_calib, sample_weight=sw)
                y_prob = cal_clf.predict_proba(this_X_test)
                y_pred = cal_clf.predict(this_X_test)
                prob_pos_cal_clf = y_prob[:, 1]
                assert_array_equal(y_pred,
                                   np.array([0, 1])[np.argmax(y_prob, axis=1)])

                assert (brier_score_loss(y_test, prob_pos_clf) >
                        brier_score_loss(y_test, prob_pos_cal_clf))


@pytest.mark.parametrize('method', ['sigmoid', 'isotonic'])
def test_calibration_ensemble_false(data, method):
    # Test that `ensemble=False` is the same as using predictions from
    # `cross_val_predict` to train calibrator.
    X, y = data
    clf = LinearSVC(random_state=7)

    cal_clf = CalibratedClassifierCV(clf, method=method, cv=3, ensemble=False)
    cal_clf.fit(X, y)
    cal_probas = cal_clf.predict_proba(X)

    # Get probas manually
    unbiased_preds = cross_val_predict(
        clf, X, y, cv=3, method='decision_function'
    )
    if method == 'isotonic':
        calibrator = IsotonicRegression(out_of_bounds='clip')
    else:
        calibrator = _SigmoidCalibration()
    calibrator.fit(unbiased_preds, y)
    # Use `clf` fit on all data
    clf.fit(X, y)
    clf_df = clf.decision_function(X)
    manual_probas = calibrator.predict(clf_df)
    assert_allclose(cal_probas[:, 1], manual_probas)


def test_sigmoid_calibration():
    """Test calibration values with Platt sigmoid model"""
    exF = np.array([5, -4, 1.0])
    exY = np.array([1, -1, -1])
    # computed from my python port of the C++ code in LibSVM
    AB_lin_libsvm = np.array([-0.20261354391187855, 0.65236314980010512])
    assert_array_almost_equal(AB_lin_libsvm,
                              _sigmoid_calibration(exF, exY), 3)
    lin_prob = 1. / (1. + np.exp(AB_lin_libsvm[0] * exF + AB_lin_libsvm[1]))
    sk_prob = _SigmoidCalibration().fit(exF, exY).predict(exF)
    assert_array_almost_equal(lin_prob, sk_prob, 6)

    # check that _SigmoidCalibration().fit only accepts 1d array or 2d column
    # arrays
    assert_raises(ValueError, _SigmoidCalibration().fit,
                  np.vstack((exF, exF)), exY)


def test_calibration_curve():
    """Check calibration_curve function"""
    y_true = np.array([0, 0, 0, 1, 1, 1])
    y_pred = np.array([0., 0.1, 0.2, 0.8, 0.9, 1.])
    prob_true, prob_pred = calibration_curve(y_true, y_pred, n_bins=2)
    prob_true_unnormalized, prob_pred_unnormalized = \
        calibration_curve(y_true, y_pred * 2, n_bins=2, normalize=True)
    assert len(prob_true) == len(prob_pred)
    assert len(prob_true) == 2
    assert_almost_equal(prob_true, [0, 1])
    assert_almost_equal(prob_pred, [0.1, 0.9])
    assert_almost_equal(prob_true, prob_true_unnormalized)
    assert_almost_equal(prob_pred, prob_pred_unnormalized)

    # probabilities outside [0, 1] should not be accepted when normalize
    # is set to False
    assert_raises(ValueError, calibration_curve, [1.1], [-0.1],
                  normalize=False)

    # test that quantiles work as expected
    y_true2 = np.array([0, 0, 0, 0, 1, 1])
    y_pred2 = np.array([0., 0.1, 0.2, 0.5, 0.9, 1.])
    prob_true_quantile, prob_pred_quantile = calibration_curve(
        y_true2, y_pred2, n_bins=2, strategy='quantile')

    assert len(prob_true_quantile) == len(prob_pred_quantile)
    assert len(prob_true_quantile) == 2
    assert_almost_equal(prob_true_quantile, [0, 2 / 3])
    assert_almost_equal(prob_pred_quantile, [0.1, 0.8])

    # Check that error is raised when invalid strategy is selected
    assert_raises(ValueError, calibration_curve, y_true2, y_pred2,
                  strategy='percentile')


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_nan_imputer(ensemble):
    """Test that calibration can accept nan"""
    X, y = make_classification(n_samples=10, n_features=2,
                               n_informative=2, n_redundant=0,
                               random_state=42)
    X[0, 0] = np.nan
    clf = Pipeline(
        [('imputer', SimpleImputer()),
         ('rf', RandomForestClassifier(n_estimators=1))])
    clf_c = CalibratedClassifierCV(
        clf, cv=2, method='isotonic', ensemble=ensemble
    )
    clf_c.fit(X, y)
    clf_c.predict(X)


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_prob_sum(ensemble):
    # Test that sum of probabilities is 1. A non-regression test for
    # issue #7796
    num_classes = 2
    X, y = make_classification(n_samples=10, n_features=5,
                               n_classes=num_classes)
    clf = LinearSVC(C=1.0, random_state=7)
    clf_prob = CalibratedClassifierCV(
        clf, method="sigmoid", cv=LeaveOneOut(), ensemble=ensemble
    )
    clf_prob.fit(X, y)

    probs = clf_prob.predict_proba(X)
    assert_array_almost_equal(probs.sum(axis=1), np.ones(probs.shape[0]))


@pytest.mark.parametrize('ensemble', [True, False])
def test_calibration_less_classes(ensemble):
    # Test to check calibration works fine when train set in a test-train
    # split does not contain all classes
    # Since this test uses LOO, at each iteration train set will not contain a
    # class label
    X = np.random.randn(10, 5)
    y = np.arange(10)
    clf = LinearSVC(C=1.0, random_state=7)
    cal_clf = CalibratedClassifierCV(
        clf, method="sigmoid", cv=LeaveOneOut(), ensemble=ensemble
    )
    cal_clf.fit(X, y)

    for i, calibrated_classifier in \
            enumerate(cal_clf.calibrated_classifiers_):
        proba = calibrated_classifier.predict_proba(X)
        if ensemble:
            # Check that the unobserved class has proba=0
            assert_array_equal(proba[:, i], np.zeros(len(y)))
            # Check for all other classes proba>0
            assert np.all(proba[:, :i] > 0)
            assert np.all(proba[:, i + 1:] > 0)
        else:
            # Check `proba` are all 1/n_classes
            assert np.allclose(proba, 1 / proba.shape[0])


@ignore_warnings(category=FutureWarning)
@pytest.mark.parametrize('X', [np.random.RandomState(42).randn(15, 5, 2),
                               np.random.RandomState(42).randn(15, 5, 2, 6)])
def test_calibration_accepts_ndarray(X):
    """Test that calibration accepts n-dimensional arrays as input"""
    y = [1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0]

    class MockTensorClassifier(BaseEstimator):
        """A toy estimator that accepts tensor inputs"""

        def fit(self, X, y):
            self.classes_ = np.unique(y)
            return self

        def decision_function(self, X):
            # toy decision function that just needs to have the right shape:
            return X.reshape(X.shape[0], -1).sum(axis=1)

    calibrated_clf = CalibratedClassifierCV(MockTensorClassifier())
    # we should be able to fit this classifier with no error
    calibrated_clf.fit(X, y)


@pytest.fixture
def text_data():
    text_data = [
        {'state': 'NY', 'age': 'adult'},
        {'state': 'TX', 'age': 'adult'},
        {'state': 'VT', 'age': 'child'},
    ]
    text_labels = [1, 0, 1]
    return text_data, text_labels


@pytest.fixture
def text_data_pipeline(text_data):
    X, y = text_data
    pipeline_prefit = Pipeline([
        ('vectorizer', DictVectorizer()),
        ('clf', RandomForestClassifier())
    ])
    return pipeline_prefit.fit(X, y)


def test_calibration_pipeline(text_data, text_data_pipeline):
    # Test that calibration works in prefit pipeline with transformer,
    # where `X` is not array-like, sparse matrix or dataframe at the start.
    # See https://github.com/scikit-learn/scikit-learn/issues/8710
    X, y = text_data
    clf = text_data_pipeline
    calib_clf = CalibratedClassifierCV(clf, cv='prefit')
    calib_clf.fit(X, y)
    # Check attributes are obtained from fitted estimator
    assert_array_equal(calib_clf.classes_, clf.classes_)
    msg = "'CalibratedClassifierCV' object has no attribute"
    with pytest.raises(AttributeError, match=msg):
        calib_clf.n_features_in_


@pytest.mark.parametrize('clf, cv', [
    pytest.param(LinearSVC(C=1), 2),
    pytest.param(LinearSVC(C=1), 'prefit'),
])
def test_calibration_attributes(clf, cv):
    # Check that `n_features_in_` and `classes_` attributes created properly
    X, y = make_classification(n_samples=10, n_features=5,
                               n_classes=2, random_state=7)
    if cv == 'prefit':
        clf = clf.fit(X, y)
    calib_clf = CalibratedClassifierCV(clf, cv=cv)
    calib_clf.fit(X, y)

    if cv == 'prefit':
        assert_array_equal(calib_clf.classes_, clf.classes_)
        assert calib_clf.n_features_in_ == clf.n_features_in_
    else:
        classes = LabelEncoder().fit(y).classes_
        assert_array_equal(calib_clf.classes_, classes)
        assert calib_clf.n_features_in_ == X.shape[1]


# FIXME: remove in 1.1
def test_calibrated_classifier_cv_deprecation(data):
    # Check that we raise the proper deprecation warning if accessing
    # `calibrators_` from the `_CalibratedClassifier`.
    X, y = data
    calib_clf = CalibratedClassifierCV(cv=2).fit(X, y)

    with pytest.warns(FutureWarning):
        calibrators = calib_clf.calibrated_classifiers_[0].calibrators_

    for clf1, clf2 in zip(
        calibrators, calib_clf.calibrated_classifiers_[0].calibrators
    ):
        assert clf1 is clf2
