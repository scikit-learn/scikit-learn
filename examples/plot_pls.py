import numpy as np
import numpy.random as rnd
from scikits.learn.datasets import load_linnerud
from scikits.learn.pls import PLS
from scikits.learn.pls import  center_scale_xy

d=load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# retrieve PLS properties
# ~~~~~~~~~~~~~~~~~~~~~~~
plsca = PLS(deflation_mode="canonical")
plsca.fit(X,Y, n_components=3)

T  = plsca.x_scores_
P  = plsca.x_loadings_
Wx = plsca.x_weights_
U  = plsca.y_scores_
Q  = plsca.y_loadings_
Wy = plsca.y_weights_

# On centered data, ()with all components):
# X = TP' + E
# Y = UQ' + F
Xc, Yc, x_mean, y_mean, x_std, y_std =\
         center_scale_xy(X.copy(), Y.copy(), scale=True)
         
np.round(Xc - np.dot(T, P.T), decimals=2)
np.round(Yc - np.dot(U, Q.T), decimals=2)

# scores orthogonality
np.round(np.dot(T.T, T), decimals=2)
np.round(np.dot(U.T, U), decimals=2)

# weights orthogonality
np.round(np.dot(Wx.T, Wx), decimals=2)
np.round(np.dot(Wy.T, Wy), decimals=2)

# T = XW(P'W)^−1 = XW∗
Wxs = np.dot(Wx, linalg.inv(np.dot(P.T,Wx)))
np.round(T - np.dot(Xc, Wxs), decimals=2)

# Transform data
# ~~~~~~~~~~~~~~~~~~~~~
plsca.fit(X,Y, n_components=2)

# with training data check that we get the scores
Xr, Yr = plsca.transform(X, Y)

print Xr - plsca.x_scores_
print Yr - plsca.y_scores_

Xnew = rnd.normal(size=5*X.shape[1]).reshape((5,X.shape[1]))
Ynew = rnd.normal(size=5*Y.shape[1]).reshape((5,Y.shape[1]))

# 
Xnew_r, Ynew_r = plsca.transform(Xnew, Ynew)

np.dot(Xnew, Wxs)

np.dot(Xc, Wxs)

## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
## -------------------------------------------------------------

pls2 = PLS(deflation_mode="regression")
pls2.fit(X,Y, n_components=2)


