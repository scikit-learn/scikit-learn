from sklearn.linear_model import PCR
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn import datasets
import numpy as np


def test_pcrRegressor_fit():
    diabetes = datasets.load_diabetes()
    X = diabetes.data
    Y = diabetes.target

    pcr = PCR(n_components=7)
    pcr.fit(X, Y)

    assert_array_almost_equal(
        [round(elem, 7) for elem in pcr.explained_variance_[0:2]],
        [0.009125, 0.003384])
    assert_almost_equal(pcr.coef_[0], 448.1948567)


def test_pcrRegressor_predict():
    diabetes = datasets.load_diabetes()
    X = diabetes.data
    Y = diabetes.target

    pcr = PCR(n_components=7)
    pcr.fit(X, Y)

    predictions = pcr.predict(X)
    mae = np.mean(abs(predictions - Y))
    assert_almost_equal(mae, 43.423082022528)


test_pcrRegressor_fit()
test_pcrRegressor_predict()
