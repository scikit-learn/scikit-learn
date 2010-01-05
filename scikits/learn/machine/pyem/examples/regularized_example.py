#! /usr/bin/env python
# Last Change: Mon Jul 02 04:00 PM 2007 J

__doc__ = """Example of using RegularizedEM with pendigits data.

If you want to do discriminant analysis with pendigits, you quickly have
problems with EM if you use directly the coordinates, because some points are
likely to be on the border, hence the corresponding component can have a
covariance matrix which easily becomes singular. Regularized EM avoids this
problem by using simple regularization on the mixture. You can play with pcount
and pval to see the effect on pdf estimation.

For now, regularized EM is pretty crude, but is enough for simple cases where
you need to avoid singular covariance matrices."""

import numpy as N
import pylab as P

#from scipy.sandbox import pyem
from gauss_mix import GM
from gmm_em import GMM, EM, RegularizedEM
import utils

x, y = utils.get_pendigits()

# Take only the first point of pendigits for pdf estimation
dt1 = N.concatenate([x[:, :1], y[:, :1]], 1)
dt1 = utils.scale(dt1.astype(N.float))

# pcnt is the poportion of samples to use as prior count. Eg if you have 1000
# samples, and pcn is 0.1, then the prior count would be 100, and 1100 samples
# will be considered as overall when regularizing the parameters.
pcnt = 0.05
# You should try different values of pval. If pval is 0.1, then the
# regularization will be strong. If you use something like 0.01, really sharp
# components will appear. If the values are too small, the regularizer may not
# work (singular covariance matrices).
pval = 0.05

# This function train a mixture model with k components, returns the trained
# model and the BIC
def cluster(data, k, mode = 'full'):
    d = data.shape[1]
    gm = GM(d, k, mode)
    gmm = GMM(gm, 'random')
    em = RegularizedEM(pcnt = pcnt, pval = pval)
    em.train(data, gmm, maxiter = 20)
    return gm, gmm.bic(data)

# bc will contain a list of BIC values for each model trained
N.seterr(all = 'warn')
bc = []
mode = 'full'

P.figure()
for k in range(1, 5):
    # Train a model of k component, and plots isodensity curve
    P.subplot(2, 2, k)
    gm, b = cluster(dt1, k = k, mode = mode)
    bc.append(b)

    X, Y, Z, V = gm.density_on_grid(nl = 20)
    P.contour(X, Y, Z, V)
    P.plot(dt1[:, 0], dt1[:, 1], '.')
    P.xlabel('x coordinate (scaled)')
    P.ylabel('y coordinate (scaled)')

print "According to the BIC, model with %d components is better" % (N.argmax(bc) + 1)
P.show()
