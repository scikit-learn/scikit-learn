import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LinearRegression

from sklearn.inspection import plot_individual_conditional_expectation


dummy_classification_data = make_classification(random_state=0)


@pytest.mark.parametrize(
    "data, params, err_msg",
    [(dummy_classification_data, {'features': [(1, 2)]},
      'Each entry in features must be either an int '),
     (dummy_classification_data, {'features': [1, {}]},
      'Each entry in features must be either an int ')]
)
def test_plot_ice_error(pyplot, data, params, err_msg):
    X, y = data
    estimator = LinearRegression().fit(X, y)

    with pytest.raises(ValueError, match=err_msg):
        plot_individual_conditional_expectation(estimator, X, **params)
