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
    print [np.corrcoef(m1.x_loadings_[:,k], m2.x_loadings_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations Y loadings"
    print [np.corrcoef(m1.y_loadings_[:,k], m2.y_loadings_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations X scores"
    print [np.corrcoef(m1.x_scores_[:,k], m2.x_scores_[:,k])[0,1] for k in xrange(n_components)]
    print "correlations Y scores"
    print [np.corrcoef(m1.y_scores_[:,k], m2.y_scores_[:,k])[0,1] for k in xrange(n_components)]

import numpy as np
from scikits.learn.datasets import load_linnerud
from scikits.learn.pls import PLS

d=load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# algo nipals vs svd
# ~~~~~~~~~~~~~~~~~~

pls_bynpl = PLS()
pls_bynpl.fit(X,Y, n_components=2)
check(pls_bynpl)


pls_bysvd = PLS(algorithm="svd")
pls_bysvd.fit(X,Y, n_components=2)
check(pls_bysvd)

compare(pls_bynpl, pls_bysvd)

# Non regression
# ~~~~~~~~~~~~~~

pls = PLS(deflation_mode="canonical")
pls.fit(X,Y, n_components=2)

print pls.x_loadings_
[[-0.58989155 -0.78900503]
 [-0.77134037  0.61351764]
 [ 0.23887653  0.03266757]]

print pls.y_loadings_
[[ 0.61330741 -0.25616063]
 [ 0.7469717  -0.11930623]
 [ 0.25668522  0.95924333]]


# check orthogonality of latent scores
print np.corrcoef(pls.x_scores_,rowvar=0)
[[  1.00000000e+00   2.51221165e-17]
 [  2.51221165e-17   1.00000000e+00]]

print np.corrcoef(pls.y_scores_,rowvar=0)
[[  1.00000000e+00  -8.57631722e-17]
 [ -8.57631722e-17   1.00000000e+00]]

## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
## -------------------------------------------------------------

pls2 = PLS(deflation_mode="regression")
pls2.fit(X,Y, n_components=2)

print pls2.x_loadings_
[[-0.58989155  0.46874883]
 [-0.77134037 -0.56798901]
 [ 0.23887653 -0.67650796]]

print pls2.y_loadings_
[[ 0.61330741  0.74851609]
 [ 0.7469717   0.64704492]
 [ 0.25668522  0.14510869]]


