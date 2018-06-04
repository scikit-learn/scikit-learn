import pytest
pytest.importorskip('sklearn')

from sklearn.linear_model import SGDClassifier
import dask.array as da
import numpy as np
import dask


x = np.array([[1, 0],
              [2, 0],
              [3, 0],
              [4, 0],
              [0, 1],
              [0, 2],
              [3, 3],
              [4, 4]])

y = np.array([1, 1, 1, 1, -1, -1, 0, 0])

z = np.array([[1, -1],
              [-1, 1],
              [10, -10],
              [-10, 10]])

X = da.from_array(x, chunks=(3, 2))
Y = da.from_array(y, chunks=(3,))
Z = da.from_array(z, chunks=(2, 2))


@pytest.mark.skipif(reason="Hangs")
def test_fit():
    sgd = SGDClassifier()

    sgd = da.learn.fit(sgd, X, Y, get=dask.get, classes=np.array([-1, 0, 1]))

    sol = sgd.predict(z)
    result = da.learn.predict(sgd, Z)
    assert result.chunks == ((2, 2),)
    assert result.compute(get=dask.get).tolist() == sol.tolist()
