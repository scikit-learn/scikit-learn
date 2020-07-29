import pytest

from numpy.testing import assert_allclose

from sklearn.datasets import load_diabetes
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

from sklearn.metrics import plot_prediction_error

X, y = load_diabetes(return_X_y=True)


@pytest.fixture
def regressor_fitted():
    return Ridge().fit(X, y)


@pytest.mark.parametrize(
    "regressor, params, err_type, err_msg",
    [(Ridge(), {}, NotFittedError, "is not fitted yet"),
     (Ridge().fit(X, y), {"subsample": -1}, ValueError,
      "When an integer, subsample=-1 should be"),
     (Ridge().fit(X, y), {"subsample": 20.0}, ValueError,
      "When a floating-point, subsample=20.0 should be"),
     (Ridge().fit(X, y), {"subsample": -20.0}, ValueError,
      "When a floating-point, subsample=-20.0 should be")]
)
def test_plot_prediction_error_raise_error(
    pyplot, regressor, params, err_type, err_msg
):
    # check that we raise the proper error when making the parameters
    # validation
    with pytest.raises(err_type, match=err_msg):
        plot_prediction_error(regressor, X, y, **params)


def test_plot_prediction_error(pyplot, regressor_fitted):
    # general check regarding the plot_prediction_error function
    disp = plot_prediction_error(regressor_fitted, X, y)

    assert_allclose(disp.line_.get_xdata(), disp.line_.get_ydata())
    assert disp.ax_.get_xlabel() == "Actual values"
    assert disp.ax_.get_ylabel() == "Predicted values"
    assert disp.residual_lines_ is None

    legend_text = disp.ax_.get_legend().get_texts()
    assert len(legend_text) == 1
    assert legend_text[0].get_text().startswith("r2 =")


def test_plot_prediction_error_scoring(pyplot, regressor_fitted):
    # check that the legend is properly display when scoring is passed

    # default behaviour -> r2 score will be computed
    disp = plot_prediction_error(regressor_fitted, X, y, scoring=None)
    legend_text = disp.ax_.get_legend().get_texts()[0].get_text()
    assert legend_text.startswith("r2 =")
    pyplot.close("all")

    # passing a single callable
    def my_scorer(estimator, X, y):
        y_pred = estimator.predict(X)
        return r2_score(y, y_pred)

    disp = plot_prediction_error(regressor_fitted, X, y, scoring=my_scorer)
    legend_text = disp.ax_.get_legend().get_texts()[0].get_text()
    assert legend_text.startswith("my scorer =")
    pyplot.close("all")

    # passing a list of scorer
    disp = plot_prediction_error(
        regressor_fitted, X, y, scoring=["r2", "neg_median_absolute_error"]
    )
    legend_text = disp.ax_.get_legend().get_texts()[0].get_text()
    for metric_name in ["r2 =", "neg median absolute error ="]:
        assert metric_name in legend_text


def test_plot_prediction_error_subsample(pyplot, regressor_fitted):
    # check the behaviour of subsample
    disp = plot_prediction_error(regressor_fitted, X, y, subsample=5)
