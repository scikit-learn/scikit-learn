# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import sparse

from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.model_selection import LeaveOneOut

from sklearn.utils._testing import (assert_array_almost_equal,
                                   assert_almost_equal,
                                   assert_array_equal,
                                   assert_raises, ignore_warnings)
from sklearn.exceptions import NotFittedError
from sklearn.datasets import make_classification, make_blobs, load_iris
from sklearn.naive_bayes import MultinomialNB
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import make_column_transformer
from sklearn.impute import SimpleImputer
from sklearn.metrics import brier_score_loss, log_loss
from sklearn.calibration import CalibratedClassifierCV
from sklearn.calibration import _sigmoid_calibration, _SigmoidCalibration
from sklearn.calibration import (calibration_curve, plot_calibration_curve,
                                 CalibrationDisplay)


def test_calibration():
    """Test calibration objects with isotonic and sigmoid"""
    n_samples = 100
    X, y = make_classification(n_samples=2 * n_samples, n_features=6,
                               random_state=42)
    sample_weight = np.random.RandomState(seed=42).uniform(size=y.size)

    X -= X.min()  # MultinomialNB only allows positive X

    # split train and test
    X_train, y_train, sw_train = \
        X[:n_samples], y[:n_samples], sample_weight[:n_samples]
    X_test, y_test = X[n_samples:], y[n_samples:]

    # Naive-Bayes
    clf = MultinomialNB().fit(X_train, y_train, sample_weight=sw_train)
    prob_pos_clf = clf.predict_proba(X_test)[:, 1]

    pc_clf = CalibratedClassifierCV(clf, cv=y.size + 1)
    assert_raises(ValueError, pc_clf.fit, X, y)

    # Naive Bayes with calibration
    for this_X_train, this_X_test in [(X_train, X_test),
                                      (sparse.csr_matrix(X_train),
                                       sparse.csr_matrix(X_test))]:
        for method in ['isotonic', 'sigmoid']:
            pc_clf = CalibratedClassifierCV(clf, method=method, cv=2)
            # Note that this fit overwrites the fit on the entire training
            # set
            pc_clf.fit(this_X_train, y_train, sample_weight=sw_train)
            prob_pos_pc_clf = pc_clf.predict_proba(this_X_test)[:, 1]

            # Check that brier score has improved after calibration
            assert (brier_score_loss(y_test, prob_pos_clf) >
                           brier_score_loss(y_test, prob_pos_pc_clf))

            # Check invariance against relabeling [0, 1] -> [1, 2]
            pc_clf.fit(this_X_train, y_train + 1, sample_weight=sw_train)
            prob_pos_pc_clf_relabeled = pc_clf.predict_proba(this_X_test)[:, 1]
            assert_array_almost_equal(prob_pos_pc_clf,
                                      prob_pos_pc_clf_relabeled)

            # Check invariance against relabeling [0, 1] -> [-1, 1]
            pc_clf.fit(this_X_train, 2 * y_train - 1, sample_weight=sw_train)
            prob_pos_pc_clf_relabeled = pc_clf.predict_proba(this_X_test)[:, 1]
            assert_array_almost_equal(prob_pos_pc_clf,
                                      prob_pos_pc_clf_relabeled)

            # Check invariance against relabeling [0, 1] -> [1, 0]
            pc_clf.fit(this_X_train, (y_train + 1) % 2,
                       sample_weight=sw_train)
            prob_pos_pc_clf_relabeled = \
                pc_clf.predict_proba(this_X_test)[:, 1]
            if method == "sigmoid":
                assert_array_almost_equal(prob_pos_pc_clf,
                                          1 - prob_pos_pc_clf_relabeled)
            else:
                # Isotonic calibration is not invariant against relabeling
                # but should improve in both cases
                assert (brier_score_loss(y_test, prob_pos_clf) >
                               brier_score_loss((y_test + 1) % 2,
                                                prob_pos_pc_clf_relabeled))

        # Check failure cases:
        # only "isotonic" and "sigmoid" should be accepted as methods
        clf_invalid_method = CalibratedClassifierCV(clf, method="foo")
        assert_raises(ValueError, clf_invalid_method.fit, X_train, y_train)

        # base-estimators should provide either decision_function or
        # predict_proba (most regressors, for instance, should fail)
        clf_base_regressor = \
            CalibratedClassifierCV(RandomForestRegressor(), method="sigmoid")
        assert_raises(RuntimeError, clf_base_regressor.fit, X_train, y_train)


def test_sample_weight():
    n_samples = 100
    X, y = make_classification(n_samples=2 * n_samples, n_features=6,
                               random_state=42)

    sample_weight = np.random.RandomState(seed=42).uniform(size=len(y))
    X_train, y_train, sw_train = \
        X[:n_samples], y[:n_samples], sample_weight[:n_samples]
    X_test = X[n_samples:]

    for method in ['sigmoid', 'isotonic']:
        base_estimator = LinearSVC(random_state=42)
        calibrated_clf = CalibratedClassifierCV(base_estimator, method=method)
        calibrated_clf.fit(X_train, y_train, sample_weight=sw_train)
        probs_with_sw = calibrated_clf.predict_proba(X_test)

        # As the weights are used for the calibration, they should still yield
        # a different predictions
        calibrated_clf.fit(X_train, y_train)
        probs_without_sw = calibrated_clf.predict_proba(X_test)

        diff = np.linalg.norm(probs_with_sw - probs_without_sw)
        assert diff > 0.1


def test_calibration_multiclass():
    """Test calibration for multiclass """
    # test multi-class setting with classifier that implements
    # only decision function
    clf = LinearSVC()
    X, y_idx = make_blobs(n_samples=100, n_features=2, random_state=42,
                          centers=3, cluster_std=3.0)

    # Use categorical labels to check that CalibratedClassifierCV supports
    # them correctly
    target_names = np.array(['a', 'b', 'c'])
    y = target_names[y_idx]

    X_train, y_train = X[::2], y[::2]
    X_test, y_test = X[1::2], y[1::2]

    clf.fit(X_train, y_train)
    for method in ['isotonic', 'sigmoid']:
        cal_clf = CalibratedClassifierCV(clf, method=method, cv=2)
        cal_clf.fit(X_train, y_train)
        probas = cal_clf.predict_proba(X_test)
        assert_array_almost_equal(np.sum(probas, axis=1), np.ones(len(X_test)))

        # Check that log-loss of calibrated classifier is smaller than
        # log-loss of naively turned OvR decision function to probabilities
        # via softmax
        def softmax(y_pred):
            e = np.exp(-y_pred)
            return e / e.sum(axis=1).reshape(-1, 1)

        uncalibrated_log_loss = \
            log_loss(y_test, softmax(clf.decision_function(X_test)))
        calibrated_log_loss = log_loss(y_test, probas)
        assert uncalibrated_log_loss >= calibrated_log_loss

    # Test that calibration of a multiclass classifier decreases log-loss
    # for RandomForestClassifier
    X, y = make_blobs(n_samples=100, n_features=2, random_state=42,
                      cluster_std=3.0)
    X_train, y_train = X[::2], y[::2]
    X_test, y_test = X[1::2], y[1::2]

    clf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf.fit(X_train, y_train)
    clf_probs = clf.predict_proba(X_test)
    loss = log_loss(y_test, clf_probs)

    for method in ['isotonic', 'sigmoid']:
        cal_clf = CalibratedClassifierCV(clf, method=method, cv=3)
        cal_clf.fit(X_train, y_train)
        cal_clf_probs = cal_clf.predict_proba(X_test)
        cal_loss = log_loss(y_test, cal_clf_probs)
        assert loss > cal_loss


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
    clf.fit(X_train, y_train, sw_train)
    prob_pos_clf = clf.predict_proba(X_test)[:, 1]

    # Naive Bayes with calibration
    for this_X_calib, this_X_test in [(X_calib, X_test),
                                      (sparse.csr_matrix(X_calib),
                                       sparse.csr_matrix(X_test))]:
        for method in ['isotonic', 'sigmoid']:
            pc_clf = CalibratedClassifierCV(clf, method=method, cv="prefit")

            for sw in [sw_calib, None]:
                pc_clf.fit(this_X_calib, y_calib, sample_weight=sw)
                y_prob = pc_clf.predict_proba(this_X_test)
                y_pred = pc_clf.predict(this_X_test)
                prob_pos_pc_clf = y_prob[:, 1]
                assert_array_equal(y_pred,
                                   np.array([0, 1])[np.argmax(y_prob, axis=1)])

                assert (brier_score_loss(y_test, prob_pos_clf) >
                               brier_score_loss(y_test, prob_pos_pc_clf))


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


def test_calibration_nan_imputer():
    """Test that calibration can accept nan"""
    X, y = make_classification(n_samples=10, n_features=2,
                               n_informative=2, n_redundant=0,
                               random_state=42)
    X[0, 0] = np.nan
    clf = Pipeline(
        [('imputer', SimpleImputer()),
         ('rf', RandomForestClassifier(n_estimators=1))])
    clf_c = CalibratedClassifierCV(clf, cv=2, method='isotonic')
    clf_c.fit(X, y)
    clf_c.predict(X)


def test_calibration_prob_sum():
    # Test that sum of probabilities is 1. A non-regression test for
    # issue #7796
    num_classes = 2
    X, y = make_classification(n_samples=10, n_features=5,
                               n_classes=num_classes)
    clf = LinearSVC(C=1.0)
    clf_prob = CalibratedClassifierCV(clf, method="sigmoid", cv=LeaveOneOut())
    clf_prob.fit(X, y)

    probs = clf_prob.predict_proba(X)
    assert_array_almost_equal(probs.sum(axis=1), np.ones(probs.shape[0]))


def test_calibration_less_classes():
    # Test to check calibration works fine when train set in a test-train
    # split does not contain all classes
    # Since this test uses LOO, at each iteration train set will not contain a
    # class label
    X = np.random.randn(10, 5)
    y = np.arange(10)
    clf = LinearSVC(C=1.0)
    cal_clf = CalibratedClassifierCV(clf, method="sigmoid", cv=LeaveOneOut())
    cal_clf.fit(X, y)

    for i, calibrated_classifier in \
            enumerate(cal_clf.calibrated_classifiers_):
        proba = calibrated_classifier.predict_proba(X)
        assert_array_equal(proba[:, i], np.zeros(len(y)))
        assert np.all(np.hstack([proba[:, :i],
                                 proba[:, i + 1:]]))


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


@pytest.fixture(scope="module")
def data():
    return load_iris(return_X_y=True)


@pytest.fixture(scope="module")
def data_binary(data):
    X, y = data
    return X[y < 2], y[y < 2]


def test_plot_calibration_curve_error_non_binary(pyplot, data):
    X, y = data
    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=msg):
        plot_calibration_curve(clf, X, y)


def test_plot_calibration_curve_no_method(pyplot, data_binary):
    X, y = data_binary

    class MyClassifier(ClassifierMixin):
        def fit(self, X, y):
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    msg = ("Neither response method 'predict_proba' nor 'decision_function' "
           "are defined in")
    with pytest.raises(ValueError, match=msg):
        plot_calibration_curve(clf, X, y)


def test_plot_calibration_curve_not_fitted(pyplot, data_binary):
    X, y = data_binary
    clf = LogisticRegression()

    with pytest.raises(NotFittedError):
        plot_calibration_curve(clf, X, y)


@pytest.mark.parametrize("n_bins", [5, 10])
@pytest.mark.parametrize("strategy", ["uniform", "quantile"])
@pytest.mark.parametrize("brier_score", [True, False])
@pytest.mark.parametrize("with_strings", [True, False])
def test_plot_calibration_curve(pyplot, data_binary, n_bins, strategy,
                                brier_score, with_strings):
    X, y = data_binary

    pos_label = None
    if with_strings:
        y = np.array(["c", "b"])[y]
        pos_label = "c"

    lr = LogisticRegression().fit(X, y)

    viz = plot_calibration_curve(
        lr, X, y, n_bins=n_bins, strategy=strategy, brier_score=brier_score,
        alpha=0.8
    )

    y_prob = lr.predict_proba(X)[:, 1]
    prob_true, prob_pred = calibration_curve(
        y, y_prob, n_bins=n_bins, strategy=strategy
    )

    assert_allclose(viz.prob_true, prob_true)
    assert_allclose(viz.prob_pred, prob_pred)
    assert_allclose(viz.y_prob, y_prob)

    assert viz.estimator_name == "LogisticRegression"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqa
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert viz.line_.get_alpha() == 0.8
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    assert viz.ax_.get_xlabel() == "Mean predicted probability"
    assert viz.ax_.get_ylabel() == "Fraction of positives"

    if brier_score:
        brier_value = brier_score_loss(
            y, y_prob, pos_label=pos_label
        )
        assert_allclose(brier_value, viz.brier_value)
        expected_label = \
            f"LogisticRegression (Brier: {viz.brier_value:.3f})"
        assert viz.line_.get_label() == expected_label
    else:
        assert viz.line_.get_label() == "LogisticRegression"


def test_plot_calibration_curve_decision_function(pyplot, data_binary):
    # Check when `estimator` has no `predict_proba` and `decision_function`
    # used
    X, y = data_binary
    svc = LinearSVC().fit(X, y)

    viz = plot_calibration_curve(svc, X, y)

    y_prob = svc.decision_function(X)
    y_prob_pos = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())
    prob_true, prob_pred = calibration_curve(y, y_prob_pos)

    assert_allclose(viz.prob_true, prob_true)
    assert_allclose(viz.prob_pred, prob_pred)
    assert_allclose(viz.y_prob, y_prob_pos)

    assert viz.estimator_name == "LinearSVC"

    # cannot fail thanks to pyplot fixture
    import matplotlib as mpl  # noqa
    assert isinstance(viz.line_, mpl.lines.Line2D)
    assert isinstance(viz.ax_, mpl.axes.Axes)
    assert isinstance(viz.figure_, mpl.figure.Figure)

    assert viz.ax_.get_xlabel() == "Mean predicted probability"
    assert viz.ax_.get_ylabel() == "Fraction of positives"


@pytest.mark.parametrize(
    "clf", [make_pipeline(StandardScaler(), LogisticRegression()),
            make_pipeline(make_column_transformer((StandardScaler(), [0, 1])),
                          LogisticRegression())])
def test_plot_calibration_curve_pipeline(pyplot, data_binary, clf):
    X, y = data_binary
    clf.fit(X, y)
    viz = plot_calibration_curve(clf, X, y)
    assert clf.__class__.__name__ in viz.line_.get_label()
    assert viz.estimator_name == clf.__class__.__name__


def test_plot_roc_curve_estimator_name_multiple_calls(pyplot, data_binary):
    # non-regression test checking that the `name` used when calling
    # `plot_calibration_curve` is used as well when calling `viz.plot()`
    X, y = data_binary
    clf_name = "my hand-crafted name"
    clf = LogisticRegression().fit(X, y)
    viz = plot_calibration_curve(clf, X, y, name=clf_name)
    assert viz.estimator_name == clf_name
    pyplot.close("all")
    viz.plot()
    assert clf_name in viz.line_.get_label()
    pyplot.close("all")
    clf_name = "another_name"
    viz.plot(name=clf_name)
    assert clf_name in viz.line_.get_label()


def test_plot_calibration_curve_ref_line(pyplot, data_binary):
    # Check that `ref_line` only appears once
    X, y = data_binary
    lr = LogisticRegression().fit(X, y)
    dt = DecisionTreeClassifier().fit(X, y)

    viz = plot_calibration_curve(lr, X, y)
    viz2 = plot_calibration_curve(dt, X, y, ax=viz.ax_)

    labels = viz2.ax_.get_legend_handles_labels()[1]
    assert labels.count('Perfectly calibrated') == 1


@pytest.mark.parametrize(
    "brier_value, estimator_name, expected_label",
    [(0.07, None, "Brier: 0.070"),
     (None, "my_est", "my_est"),
     (0.07, "my_est2", "my_est2 (Brier: 0.070)")]
)
def test_calibration_display_default_labels(pyplot, brier_value,
                                            estimator_name, expected_label):
    prob_true = np.array([0, 1, 1, 0])
    prob_pred = np.array([0.2, 0.8, 0.8, 0.4])
    y_prob = np.array([])

    viz = CalibrationDisplay(prob_true, prob_pred, y_prob,
                             brier_value=brier_value,
                             estimator_name=estimator_name)
    viz.plot()
    assert viz.line_.get_label() == expected_label
