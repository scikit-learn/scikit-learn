import numpy.random as rnd
from scikits.learn.pls import PLS
import numpy as np
import pylab as pl

## Dataset based latent variables model
## ------------------------------------

n=500
# 2 latents vars:
l1 = rnd.normal(size = n)
l2 = rnd.normal(size = n)
# X block
x0 = l1 + rnd.normal(size = n)
x1 = l1 + rnd.normal(size = n)
x2 = l2 + rnd.normal(size = n)
x3 = l2 + rnd.normal(size = n)
# Y block
y0 = l1 + rnd.normal(size = n)
y1 = l1 + rnd.normal(size = n)
y2 = l2 + rnd.normal(size = n)
y3 = l2 + rnd.normal(size = n)

X = np.array([x0, x1, x2, x3]).T
Y = np.array([y0, y1, y2, y3]).T

Xtrain = X[:n/2,:]
Ytrain = Y[:n/2,:]
Xtest  = X[n/2:,:]
Ytest  = Y[n/2:,:]

print "Corr(X)"
print np.round(np.corrcoef(X.T),2)
print "Corr(Y)"
print np.round(np.corrcoef(Y.T),2)

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# Transform data
# ~~~~~~~~~~~~~~
plsca = PLS(deflation_mode="canonical")
plsca.fit(X,Y, n_components=2)
Xtrain_r, Ytrain_r = plsca.transform(Xtrain, Ytrain)
Xtest_r, Ytest_r = plsca.transform(Xtest, Ytest)

# Scatter plot of scores
# ~~~~~~~~~~~~~~~~~~~~~~
# 1) on diagonal plot X vs Y scores on each components 
pl.subplot(221)
pl.plot(Xtrain_r[:,0], Ytrain_r[:,0], "ob",label="train")
pl.plot(Xtest_r[:,0], Ytest_r[:,0], "or",label="test")
pl.xlabel("y")
pl.ylabel("x")
pl.title('Comp. 1, corr='+\
    str(np.round(np.corrcoef(Xtest_r[:,0], Xtest_r[:,0])[0,1], decimals=2)))
pl.legend()

pl.subplot(224)
pl.plot(Xtrain_r[:,1], Ytrain_r[:,1], "ob",label="train")
pl.plot(Xtest_r[:,1], Ytest_r[:,1], "or",label="test")
pl.xlabel("y")
pl.ylabel("x")
pl.title('Comp. 2, corr='+\
    str(np.round(np.corrcoef(Xtest_r[:,1], Xtest_r[:,1])[0,1], decimals=2)))
pl.legend()


# 2) Off diagonal plot compements 1 vs 2 for X and Y 
pl.subplot(222)  
pl.plot(Xtrain_r[:,0], Xtrain_r[:,1], "*b",label="train")
pl.plot(Xtest_r[:,0], Xtest_r[:,1], "*r",label="test")
pl.xlabel("X comp. 1")
pl.ylabel("X comp. 2")
pl.title('X, corr='+\
    str(np.round(np.corrcoef(Xtest_r[:,0], Xtest_r[:,1])[0,1], decimals=2)))
pl.legend()

pl.subplot(223)
pl.plot(Ytrain_r[:,0], Ytrain_r[:,1], "*b",label="train")
pl.plot(Ytest_r[:,0], Ytest_r[:,1], "*r",label="test")
pl.xlabel("Y comp. 1")
pl.ylabel("Y comp. 2")
pl.title('Y, corr='+\
    str(np.round(np.corrcoef(Ytest_r[:,0], Ytest_r[:,1])[0,1], decimals=2)))
pl.legend()

pl.show()


## Regression PLS (PLS 2 blocks regression mode A known as PLS2)
## -------------------------------------------------------------
pls2 = PLS(deflation_mode="regression")
pls2.fit(X,Y, n_components=2)
pls2.predict(X)


