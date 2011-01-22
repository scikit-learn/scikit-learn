import pylab as pl
import numpy as np
import scipy
from scipy import linalg

from scikits.learn.datasets.samples_generator import multivariate_normal_from_latent_variables

# Check
# -----
# QC routines
def check(m):
    n_components = m.x_loadings_.shape[1]
    print "score correlations:"
    print "X scores\n",np.corrcoef(m.x_scores_,rowvar=0)
    print "Y scores\n",np.corrcoef(m.y_scores_,rowvar=0)
    print "correlations between pairs of scores"
    print [np.corrcoef(m.x_scores_[:,k], m.y_scores_[:,k])[0,1] for k in xrange(n_components)]
    

# Compare loadings and score of the two algo.
def compare(m1, m2):
    n_components = m1.x_loadings_.shape[1]
    print "correlations X loadings"
    print [scipy.corrcoef(m1.x_loadings_[:,k], m2.x_loadings_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations Y loadings"
    print [scipy.corrcoef(m1.y_loadings_[:,k], m2.y_loadings_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations X scores"
    print [scipy.corrcoef(m1.x_scores_[:,k], m2.x_scores_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations Y scores"
    print [scipy.corrcoef(m1.y_scores_[:,k], m2.y_scores_[:,k])[0,1] for k in xrange(n_components)]
 
from scikits.learn.datasets import load_linnerud
d=load_linnerud()
run /home/duchesnay/git/scikit-learn/scikits/learn/pls.py
X = d['data_exercise']
Y = d['data_physiological']

d['header_exercise']
d['header_physiological']

## Canonical
pls_npl = PLS()
pls_npl.fit(X,Y, n_components=2)
check(pls_npl)


pls_svd = PLS(algorithm="svd")
pls_svd.fit(X,Y, n_components=2)
check(pls_svd)

compare(pls_npl,pls_svd)

pls_npl.x_loadings_
array([[-0.58989155,  0.46874883],
       [-0.77134037, -0.56798901],
       [ 0.23887653, -0.67650796]])
pls_npl.y_loadings_
array([[ 0.61330741,  0.74851609],
       [ 0.7469717 ,  0.64704492],
       [ 0.25668522,  0.14510869]])

## Regression
pls_npl = PLS(deflation_mode="regression")
pls_npl.fit(X,Y, n_components=2)
check(pls_npl)
self = pls_npl


compare(pls_npl, pls_svd)


# Play
# ----

# 
latent_coefs = np.array([[1]*5 + [0]*15,
                                 [0]*5  + [1]*5 + [0]*10,
                                 [0]*10 + [1]*5 + [0]*5,
                                 [0]*15 + [1]*5,
                                 [0]*5  + [1]*10 + [0]*5])
print latent_coefs
M, latent = multivariate_normal_from_latent_variables(latent_coefs)
pl.matshow(np.corrcoef(M, rowvar=0)); pl.colorbar()

import numpy.random as rnd
n=100;p=2
#X = rnd.normal(size=n*p).reshape((n,p))
x=rnd.normal(size=n)
y=2*x + rnd.normal(size=n)
X = np.asarray([x,y]).T
X.shape
U, s, Vh = linalg.svd(X, full_matrices = False)

pl.plot(X[:,0], X[:,1], 'g.')
U, s, Vh = linalg.svd(X, full_matrices = False)
pl.plot(U[:,0], U[:,1], 'g.')

## Param
X = M[:,0:10]
Y = M[:,11:19]
n_components = 3
# scale data
X = (X - X.mean(axis=0)) / X.std(axis=0)
Y = (Y - Y.mean(axis=0)) / Y.std(axis=0)

run /home/duchesnay/git/scikit-learn/scikits/learn/pls.py

# Fit PLS using NIPALS
pls_bynipals = PLS(scale_X=False, scale_Y=False, center_X=False, center_Y=False)
pls_bynipals.fit(X,Y, n_components=n_components)

# Fit PLS using SVD
pls_bysvd = PLS(scale_X=False, scale_Y=False, center_X=False, center_Y=False, algorithm="svd")
pls_bysvd.fit(X,Y, n_components=n_components)

    
pls_svd = PLS_SVD(scale_X=False, scale_Y=False, center_X=False, center_Y=False)
pls_svd.fit(X,Y, n_components=n_components)



np.max(np.abs(pls_bynipals.x_loadings - pls_bysvd.x_loadings), axis=0)
np.max(np.abs(pls_bynipals.y_loadings - pls_bysvd.y_loadings), axis=0)
np.max(np.abs(pls_bynipals.x_scores - pls_bysvd.x_scores), axis=0)
np.max(np.abs(pls_bynipals.x_scores - pls_bysvd.x_scores), axis=0)


C = np.dot(X.T,Y)
U, s, Vh = linalg.svd(C, full_matrices = False)

cor_u = [scipy.corrcoef(U[:,k], pls_bynipals.x_loadings[:,k])[0,1] for k in xrange(n_components)]
cor_v = [scipy.corrcoef(Vh.T[:,k],pls_bynipals.y_loadings[:,k])[0,1] for k in xrange(n_components)]
cor_u = [scipy.corrcoef(U[:,k], pls_bysvd.x_loadings[:,k])[0,1] for k in xrange(n_components)]
cor_v = [scipy.corrcoef(Vh.T[:,k],pls_bysvd.y_loadings[:,k])[0,1] for k in xrange(n_components)]

print cor_u, cor_v

np.abs(cor_u)[0]>.99
np.abs(cor_v)[0]>.99

