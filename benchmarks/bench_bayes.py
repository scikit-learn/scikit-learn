"""
A comparison of different methods in linear_model methods.

Data comes from a random square matrix.

"""
from datetime import datetime
import numpy as np
from scikits.learn import linear_model
from scikits.learn.utils.bench import total_seconds


if __name__ == '__main__':

    import pylab as pl

    n_iter = 20

    time_ridge = np.empty(n_iter)
    time_ols = np.empty(n_iter)
    time_lasso = np.empty(n_iter)

    dimensions = 10 * np.arange(n_iter)

    n_samples, n_features = 100, 100

    X = np.random.randn(n_samples, n_features)
    y = np.random.randn(n_samples)

    start = datetime.now()
    ridge = linear_model.BayesianRidge()
    ridge.fit(X, y)
