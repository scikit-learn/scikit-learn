"""
A comparison of different methods in GLM

Data comes from a random square matrix.

"""
from datetime import datetime
import numpy as np
from scikits.learn import glm
from scikits.learn.utils.bench import total_seconds


if __name__ == '__main__':

    import pylab as pl

    n_iter = 20

    time_ridge   = np.empty (n_iter)
    time_ols     = np.empty (n_iter)
    time_lasso   = np.empty (n_iter)

    dimensions = 10 * np.arange(n_iter)

    for i in range(n_iter):

        print 'Iteration %s of %s' % (i, n_iter)

        n, m = 10*i + 3, 10*i + 3

        X = np.random.randn (n, m)
        Y = np.random.randn (n)

        start = datetime.now()
        ridge = glm.Ridge(alpha=0.)
        ridge.fit (X, Y)
        time_ridge[i] = total_seconds(datetime.now() - start)

        start = datetime.now()
        ols = glm.LinearRegression()
        ols.fit (X, Y)
        time_ols[i] = total_seconds(datetime.now() - start)


        start = datetime.now()
        lasso = glm.LassoLARS()
        lasso.fit (X, Y)
        time_lasso[i] = total_seconds(datetime.now() - start)



    pl.xlabel ('Dimesions')
    pl.ylabel ('Time (in seconds)')
    pl.plot (dimensions, time_ridge, color='r')
    pl.plot (dimensions, time_ols, color='g')
    pl.plot (dimensions, time_lasso, color='b')

    pl.legend (['Ridge', 'OLS', 'LassoLARS'])
    pl.axis ('tight')
    pl.show()
