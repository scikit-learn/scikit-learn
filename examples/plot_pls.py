import numpy.random as rnd
from scikits.learn.datasets import load_linnerud
from scikits.learn.pls import PLS

d=load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# Transform data
# ~~~~~~~~~~~~~~
plsca = PLS(deflation_mode="canonical")
plsca.fit(X,Y, n_components=2)

Xnew = rnd.normal(size=5*X.shape[1]).reshape((5,X.shape[1]))
Ynew = rnd.normal(size=5*Y.shape[1]).reshape((5,Y.shape[1]))

# 
Xnew_r, Ynew_r = plsca.transform(Xnew, Ynew)



## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
## -------------------------------------------------------------

pls2 = PLS(deflation_mode="regression")
pls2.fit(X,Y, n_components=3)
pls2.predict(X)

