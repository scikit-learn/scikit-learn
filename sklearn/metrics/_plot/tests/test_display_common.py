import pytest

from sklearn.datasets import make_classification
from sklearn.svm import SVC, SVR

from sklearn.metrics import ConfusionMatrixDisplay, PrecisionRecallDisplay


@pytest.mark.parametrize("Display", [ConfusionMatrixDisplay, PrecisionRecallDisplay])
def test_confusion_matrix_display_validation(pyplot, Display):
    """Check that we raise the proper error when validating parameters for
    display handling binary classification."""
    X, y = make_classification(
        n_samples=100, n_informative=5, n_classes=5, random_state=0
    )

    regressor = SVR().fit(X, y)
    y_pred_regressor = regressor.predict(X)
    y_pred_classifier = SVC().fit(X, y).predict(X)

    err_msg = f"{Display.__name__}.from_estimator only supports classifiers"
    with pytest.raises(ValueError, match=err_msg):
        Display.from_estimator(regressor, X, y)

    err_msg = "Mix type of y not allowed, got types"
    with pytest.raises(ValueError, match=err_msg):
        # Force `y_true` to be seen as a regression problem
        Display.from_predictions(y + 0.5, y_pred_classifier)
    with pytest.raises(ValueError, match=err_msg):
        Display.from_predictions(y, y_pred_regressor)

    err_msg = "Found input variables with inconsistent numbers of samples"
    with pytest.raises(ValueError, match=err_msg):
        Display.from_predictions(y, y_pred_classifier[::2])
