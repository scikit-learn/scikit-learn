import numpy as np
from sklearn.utils.testing import assert_array_almost_equal
from .._square_residue import square_residue


def square_residue_python(X):
    arr = X - X.mean(axis=1)[None].T - X.mean(axis=0) + X.mean()
    return np.power(arr, 2)


def test_square_residue():
    generator = np.random.RandomState(0)

    def make_data(n_rows, n_cols):
        X = generator.uniform(0, 100, (n_rows, n_cols))
        rows = np.arange(n_rows)
        cols = np.arange(n_cols)
        return X, rows, cols

    for shape in ((20, 20), (10, 20), (20, 10)):
        X, rows, cols = make_data(*shape)
        assert_array_almost_equal(square_residue(rows, cols, X),
                                  square_residue_python(X[rows[None].T, cols]),
                                  decimal=10)
