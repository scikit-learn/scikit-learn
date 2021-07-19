import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

from sklearn.metrics import plot_precision_recall_curve

# TODO: Remove when https://github.com/numpy/numpy/issues/14397 is resolved
pytestmark = pytest.mark.filterwarnings(
    "ignore:In future, it will be an error for 'np.bool_':DeprecationWarning:"
    "matplotlib.*"
)


# FIXME: Remove in 1.2
def test_plot_precision_recall_curve_deprecation(pyplot):
    """Check that we raise a FutureWarning when calling
    `plot_precision_recall_curve`."""

    X, y = make_classification(random_state=0)
    clf = LogisticRegression().fit(X, y)
    deprecation_warning = "Function plot_precision_recall_curve is deprecated"
    with pytest.warns(FutureWarning, match=deprecation_warning):
        plot_precision_recall_curve(clf, X, y)
