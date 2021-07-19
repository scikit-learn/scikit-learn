import pytest

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.datasets import make_classification
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, SVR

from sklearn.metrics import PrecisionRecallDisplay, plot_precision_recall_curve

# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*"
)


def test_confusion_matrix_display_validation(pyplot):
    """Check that we raise the proper error when validating parameters."""
    X, y = make_classification(
        n_samples=100, n_informative=5, n_classes=5, random_state=0
    )

    with pytest.raises(NotFittedError):
        PrecisionRecallDisplay.from_estimator(SVC(), X, y)

    regressor = SVR().fit(X, y)
    y_pred_regressor = regressor.predict(X)
    classifier = SVC(probability=True).fit(X, y)
    y_pred_classifier = classifier.predict_proba(X)[:, -1]

    err_msg = "PrecisionRecallDisplay.from_estimator only supports classifiers"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_estimator(regressor, X, y)

    err_msg = "SVC should be a binary classifier"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_estimator(classifier, X, y)

    err_msg = "{} format is not supported"
    with pytest.raises(ValueError, match=err_msg.format("continuous")):
        # Force `y_true` to be seen as a regression problem
        PrecisionRecallDisplay.from_predictions(y + 0.5, y_pred_classifier, pos_label=1)
    with pytest.raises(ValueError, match=err_msg.format("multiclass")):
        PrecisionRecallDisplay.from_predictions(y, y_pred_regressor, pos_label=1)

    err_msg = "Found input variables with inconsistent numbers of samples"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_predictions(y, y_pred_classifier[::2])

    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    y += 10
    classifier.fit(X, y)
    y_pred_classifier = classifier.predict_proba(X)[:, -1]
    err_msg = r"y_true takes value in {10, 11} and pos_label is not specified"
    with pytest.raises(ValueError, match=err_msg):
        PrecisionRecallDisplay.from_predictions(y, y_pred_classifier)


# FIXME: Remove in 1.2
def test_plot_precision_recall_curve_deprecation(pyplot):
    """Check that we raise a FutureWarning when calling
    `plot_precision_recall_curve`."""

    X, y = make_classification(random_state=0)
    clf = LogisticRegression().fit(X, y)
    deprecation_warning = "Function plot_precision_recall_curve is deprecated"
    with pytest.warns(FutureWarning, match=deprecation_warning):
        plot_precision_recall_curve(clf, X, y)


@pytest.mark.parametrize(
    "response_method, msg",
    [
        (
            "predict_proba",
            "response method predict_proba is not defined in MyClassifier",
        ),
        (
            "decision_function",
            "response method decision_function is not defined in MyClassifier",
        ),
        (
            "auto",
            "response method decision_function or predict_proba is not "
            "defined in MyClassifier",
        ),
        (
            "bad_method",
            "response_method must be 'predict_proba', 'decision_function' or 'auto'",
        ),
    ],
)
def test_precision_recall_display_bad_response(pyplot, response_method, msg):
    """Check that the proper error is raised when passing a `response_method`
    not compatible with the estimator."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    class MyClassifier(ClassifierMixin, BaseEstimator):
        def fit(self, X, y):
            self.fitted_ = True
            self.classes_ = [0, 1]
            return self

    clf = MyClassifier().fit(X, y)

    with pytest.raises(ValueError, match=msg):
        PrecisionRecallDisplay.from_estimator(
            clf, X, y, response_method=response_method
        )


@pytest.mark.parametrize("constructor_name", ["from_estimator", "from_predictions"])
@pytest.mark.parametrize("response_method", ["predict_proba", "decision_function"])
def test_precision_recall_display_plotting(pyplot, constructor_name, response_method):
    """Check the overall plotting rendering."""
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)
    classifier = LogisticRegression().fit(X, y)
    classifier.fit(X, y)

    y_pred = getattr(classifier, response_method)(X)

    # safe guard for the binary if/else construction
    assert constructor_name in ("from_estimator", "from_predictions")

    if constructor_name == "from_estimator":
        PrecisionRecallDisplay.from_estimator(classifier, X, y)
    else:
        PrecisionRecallDisplay.from_predictions(y, y_pred)
