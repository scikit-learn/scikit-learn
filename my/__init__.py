from math import sqrt
from joblib import Parallel, delayed
import scipy.sparse.linalg as sp_linalg
import numpy as np
import scipy.sparse as sp
from sklearn.datasets.samples_generator import make_sparse_uncorrelated
from sklearn.linear_model.base import LinearRegression
from timeit import Timer
from scipy import sparse
import pylab as pl


ols = LinearRegression()

X, Y = make_sparse_uncorrelated(n_samples=10000, n_features=10,random_state=0)
X = sparse.coo_matrix(X)
ols.fit(X, Y)
print sp.issparse(X)


def multi_fit(X, Y, n_jobs, parallel):
    return(ols.fit(X, Y, n_jobs, parallel))


para1 = []
para2 = []
nopara = []
featuresList = range(100, 2000, 200)
number = 10
n_features = 1000

for n_samples in featuresList:
    X, y = make_sparse_uncorrelated(n_samples=n_samples, n_features=n_features, random_state=0)
    X = sparse.coo_matrix(X)
    Y = np.vstack((y, y)).T
    t = Timer("multi_fit(X,Y,n_jobs=1, parallel=False)", "from __main__ import multi_fit, X, Y")
    nopara.append(t.timeit(number=number))

for n_samples in featuresList:
    X, y = make_sparse_uncorrelated(n_samples=n_samples, n_features=n_features, random_state=0)
    X = sparse.coo_matrix(X)
    Y = np.vstack((y, y)).T
    t = Timer("multi_fit(X,Y,n_jobs=1, parallel=True)", "from __main__ import multi_fit, X, Y")
    para1.append(t.timeit(number=number))

for n_samples in featuresList:
    X, y = make_sparse_uncorrelated(n_samples=n_samples, n_features=n_features, random_state=0)
    X = sparse.coo_matrix(X)
    Y = np.vstack((y, y)).T
    t = Timer("multi_fit(X,Y,n_jobs=2, parallel=True)", "from __main__ import multi_fit, X, Y")
    para2.append(t.timeit(number=number))


pl.plot(featuresList, para1, label='Parallel, cores = 1')
pl.plot(featuresList, para2, color='r', label='Parallel, cores = 2')
pl.plot(featuresList, nopara, color='g', label='without Parallel')
pl.xlabel('number of samples \n number of features = 1000')
pl.ylabel('runtime in sec')
pl.legend()
pl.show()
print para1, para2, nopara





