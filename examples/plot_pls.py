from scikits.learn.pls import PLS
import numpy as np
import pylab as pl

## Dataset based latent variables model
## ------------------------------------

n = 500
# 2 latents vars:
l1 = np.random.normal(size=n)
l2 = np.random.normal(size=n)

latents = np.array([l1, l1, l2, l2]).T
X = latents + np.random.normal(size=4*n).reshape((n, 4))
Y = latents + np.random.normal(size=4*n).reshape((n, 4))

X_train = X[:n/2, :]
Y_train = Y[:n/2, :]
X_test = X[n/2:, :]
Y_test = Y[n/2:, :]

print "Corr(X)"
print np.round(np.corrcoef(X.T), 2)
print "Corr(Y)"
print np.round(np.corrcoef(Y.T), 2)

## Canonical (symetric) PLS (PLS 2 blocks canonical mode A)
## --------------------------------------------------------

# Transform data
# ~~~~~~~~~~~~~~
plsca = PLS(deflation_mode="canonical")
plsca.fit(X_train, Y_train, n_components=2)
X_train_r, Y_train_r = plsca.transform(X_train, Y_train)
X_test_r, Y_test_r = plsca.transform(X_test, Y_test)

# Scatter plot of scores
# ~~~~~~~~~~~~~~~~~~~~~~
# 1) on diagonal plot X vs Y scores on each components
pl.subplot(221)
pl.plot(X_train_r[:, 0], Y_train_r[:, 0], "ob", label="train")
pl.plot(X_test_r[:, 0], Y_test_r[:, 0], "or", label="test")
pl.xlabel("y")
pl.ylabel("x")
pl.title('Comp. 1, corr = %.2f' % 
         np.corrcoef(X_test_r[:, 0], X_test_r[:, 0])[0, 1])
pl.legend()

pl.subplot(224)
pl.plot(X_train_r[:, 1], Y_train_r[:, 1], "ob", label="train")
pl.plot(X_test_r[:, 1], Y_test_r[:, 1], "or", label="test")
pl.xlabel("y")
pl.ylabel("x")
pl.title('Comp. 2, corr = %.2f' % 
         np.corrcoef(X_test_r[:, 1], X_test_r[:, 1])[0, 1])
pl.legend()


# 2) Off diagonal plot components 1 vs 2 for X and Y
pl.subplot(222)
pl.plot(X_train_r[:, 0], X_train_r[:, 1], "*b", label="train")
pl.plot(X_test_r[:, 0], X_test_r[:, 1], "*r", label="test")
pl.xlabel("X comp. 1")
pl.ylabel("X comp. 2")
pl.title('X, corr = %.2f' % np.corrcoef(X_test_r[:, 0], X_test_r[:, 1])[0, 1])
pl.legend()

pl.subplot(223)
pl.plot(Y_train_r[:, 0], Y_train_r[:, 1], "*b", label="train")
pl.plot(Y_test_r[:, 0], Y_test_r[:, 1], "*r", label="test")
pl.xlabel("Y comp. 1")
pl.ylabel("Y comp. 2")
pl.title('Y, corr = %.2f' % np.corrcoef(Y_test_r[:, 0], Y_test_r[:, 1])[0, 1])
pl.legend()
pl.show()


## Regression PLS with 2 multidimentional blocks (PLS2)
## ----------------------------------------------------
n = 1000
q = 3
p = 10
X = np.random.normal(size=n * p).reshape((n, p))
B = np.array([[1, 2] + [0] * (p - 2)] * q).T
# each Yj = 1*X1 + 2*X2 + noize
Y = np.dot(X, B) + np.random.normal(size=n * q).reshape((n, q)) + 5

pls2 = PLS(deflation_mode="regression")
pls2.fit(X, Y, n_components=3)
print "True B (such that: Y = XB + Err)"
print B
# compare pls2.coefs with B
print "Estimated B"
print np.round(pls2.coefs, 1)
pls2.predict(X)

## Regression PLS with 1 dimensional response (PLS1)
## -------------------------------------------------
n = 1000
p = 10
X = np.random.normal(size=n*p).reshape((n, p))
y = X[:, 0] + 2 * X[:, 1] + np.random.normal(size=n * 1) + 5
pls1 = PLS(deflation_mode="regression")
pls1.fit(X, y, n_components=3)
# note that the number of compements exceeds 1 (the dimension of y)
print "Estimated betas"
print np.round(pls1.coefs, 1)

## Canonical (symetric) mode B PLS (CCA)
## --------------------------------------------------------

cca = PLS(deflation_mode="canonical", mode="B")
cca.fit(X_train, Y_train, n_components=2)
X_train_r, Y_train_r = plsca.transform(X_train, Y_train)
X_test_r, Y_test_r = plsca.transform(X_test, Y_test)
