import pytest
import numpy as np


from sklearn.datasets import load_iris

from sklearn.model_selection import train_test_split

from sklearn.preprocessing import QuantileTransformer

from sklearn.utils.testing import assert_array_equal

iris = load_iris()

TRANSFORMER_HANDLING_NAN = ['QuantileTransformer']


def _generate_tuple_transformer_missing_value():
    trans_handling_nan = [globals()[trans_name]()
                          for trans_name in TRANSFORMER_HANDLING_NAN]
    return [(trans, iris.data.copy(), 15)
            for trans in trans_handling_nan]


@pytest.mark.parametrize(
    "est, X, n_missing",
    _generate_tuple_transformer_missing_value()
)
def test_missing_value_handling(est, X, n_missing):
    # check that the preprocessing method let pass nan
    rng = np.random.RandomState(42)
    X[rng.randint(X.shape[0], size=n_missing),
      rng.randint(X.shape[1], size=n_missing)] = np.nan
    X_train, X_test = train_test_split(X)
    # sanity check
    assert not np.all(np.isnan(X_train), axis=0).any()
    X_test[:, 0] = np.nan  # make sure this boundary case is tested

    Xt = est.fit(X_train).transform(X_test)
    # missing values should still be missing, and only them
    assert_array_equal(np.isnan(Xt), np.isnan(X_test))

    for i in range(X.shape[1]):
        # train only on non-NaN
        est.fit(X_train[:, [i]][~np.isnan(X_train[:, i])])
        # check transforming with NaN works even when training without NaN
        Xt_col = est.transform(X_test[:, [i]])
        assert_array_equal(Xt_col, Xt[:, [i]])
        # check non-NaN is handled as before - the 1st column is all nan
        if not np.isnan(X_test[:, i]).all():
            Xt_col_nonan = est.transform(
                X_test[:, [i]][~np.isnan(X_test[:, i])])
            assert_array_equal(Xt_col_nonan,
                               Xt_col[~np.isnan(Xt_col.squeeze())])
