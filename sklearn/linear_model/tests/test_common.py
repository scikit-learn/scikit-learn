# Author: Maria Telenczuk <maikia>
#
# License: BSD 3 clause

import pytest

import numpy as np

from scipy import sparse

from sklearn.datasets import make_regression

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model._ridge import _RidgeGCV
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import _convert_container


# FIXME: 'normalize' to be removed in 1.2
@pytest.mark.filterwarnings("ignore:'normalize' was deprecated")
@pytest.mark.parametrize(
    "estimator",
    [LinearRegression]
)
@pytest.mark.parametrize(
    "is_sparse, with_mean",
    [(True, False),
     (False, True),
     (False, False)]
)
def test_linear_model_sample_weights_normalize_in_pipeline(
        estimator, is_sparse, with_mean
):
    # Test that the results for running linear regression LinearRegression with
    # sample_weight set and with normalize set to True gives similar results as
    # LinearRegression with no normalize in a pipeline with a StandardScaler
    # and set sample_weight.
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=20, n_features=5, noise=1e-2,
                           random_state=rng)
    # make sure the data is not centered to make the problem more
    # difficult
    X += 10
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5,
                                                        random_state=rng)
    if is_sparse:
        X_train = sparse.csr_matrix(X_train)
        X_test = _convert_container(X_train, 'sparse')

    sample_weight = rng.rand(X_train.shape[0])

    # linear estimator with explicit sample_weight
    reg_with_normalize = estimator(normalize=True)
    reg_with_normalize.fit(X_train, y_train, sample_weight=sample_weight)

    # linear estimator in a pipeline
    reg_with_scaler = make_pipeline(
        StandardScaler(with_mean=with_mean),
        estimator(normalize=False)
    )
    kwargs = {reg_with_scaler.steps[-1][0] + '__sample_weight':
              sample_weight}
    reg_with_scaler.fit(X_train, y_train, **kwargs)

    y_pred_norm = reg_with_normalize.predict(X_test)
    y_pred_pip = reg_with_scaler.predict(X_test)

    assert_allclose(
        reg_with_normalize.coef_ * reg_with_scaler[0].scale_,
        reg_with_scaler[1].coef_
    )
    assert_allclose(y_pred_norm, y_pred_pip)