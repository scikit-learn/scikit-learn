import numpy as np
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from .._squared_residue import all_msr


def squared_residue_python(X):
    arr = X - X.mean(axis=1)[None].T - X.mean(axis=0) + X.mean()
    return np.power(arr, 2)


def test_squared_residue():
    generator = np.random.RandomState(0)

    def make_data(n_rows, n_cols):
        X = generator.uniform(0, 100, (n_rows, n_cols))
        rows = np.arange(n_rows)
        cols = np.arange(n_cols)
        return X, rows, cols

    for shape in ((20, 20), (10, 20), (20, 10)):
        X, rows, cols = make_data(*shape)
        rows_t = rows[None].T
        sr = squared_residue_python(X[rows_t, cols])

        msr, row_msr, col_msr = all_msr(rows, cols, X)
        assert_array_equal(row_msr, sr.mean(axis=1))
        assert_array_equal(col_msr, sr.mean(axis=0))
        assert_equal(msr, sr.mean())
