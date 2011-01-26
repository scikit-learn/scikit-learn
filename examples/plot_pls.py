import numpy as np
from scikits.learn.datasets import load_linnerud
from scikits.learn.pls import PLS
from scikits.learn.pls import  _center_scale_xy

d=load_linnerud()
X = d['data_exercise']
Y = d['data_physiological']

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# X = TP' + E
# Y = UQ' + F

## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
## -------------------------------------------------------------

pls2 = PLS(deflation_mode="regression")
pls2.fit(X,Y, n_components=2)


