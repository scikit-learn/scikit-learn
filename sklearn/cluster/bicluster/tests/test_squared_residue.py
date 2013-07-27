import numpy as np
from sklearn.utils.testing import assert_array_equal
from .._squared_residue import squared_residue


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
        assert_array_equal(squared_residue(rows, cols, X),
                           squared_residue_python(X[rows[None].T, cols]))
