"""
A comparison of different methods in GLM

Data comes from a random square matrix.

"""
from datetime import datetime
import numpy as np
from sklearn import linear_model
from sklearn.utils.bench import total_seconds


if __name__ == '__main__':

    import pylab as pl

    n_iter = 40

    time_ridge = np.empty(n_iter)
    time_ols = np.empty(n_iter)
    time_lasso = np.empty(n_iter)

    dimensions = 500 * np.arange(1, n_iter + 1)

    for i in range(n_iter):

        print('Iteration %s of %s' % (i, n_iter))

        n_samples, n_features = 10 * i + 3, 10 * i + 3

        X = np.random.randn(n_samples, n_features)
        Y = np.random.randn(n_samples)

        start = datetime.now()
        ridge = linear_model.Ridge(alpha=1.)
        ridge.fit(X, Y)
        time_ridge[i] = total_seconds(datetime.now() - start)

        start = datetime.now()
        ols = linear_model.LinearRegression()
        ols.fit(X, Y)
        time_ols[i] = total_seconds(datetime.now() - start)

        start = datetime.now()
        lasso = linear_model.LassoLars()
        lasso.fit(X, Y)
        time_lasso[i] = total_seconds(datetime.now() - start)

    pl.figure('scikit-learn GLM benchmark results')
    pl.xlabel('Dimensions')
    pl.ylabel('Time (s)')
    pl.plot(dimensions, time_ridge, color='r')
    pl.plot(dimensions, time_ols, color='g')
    pl.plot(dimensions, time_lasso, color='b')

    pl.legend(['Ridge', 'OLS', 'LassoLars'], loc='upper left')
    pl.axis('tight')
    pl.show()
