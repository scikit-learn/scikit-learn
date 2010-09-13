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

    n, m = 100, 100

    X = np.random.randn (n, m) 
    Y = np.random.randn (n)

    start = datetime.now()
    ridge = glm.BayesianRidge()
    ridge.fit (X, Y)
