#! /usr/bin/env python
# Last Change: Sat Jun 09 07:00 PM 2007 J

# Example of doing pdf estimation with EM algorithm. Requires matplotlib.
import numpy as N
import pylab as P

from scipy.sandbox import pyem
import utils

oldfaithful = utils.get_faithful()

# We want the relationship between d(t) and w(t+1), but get_faithful gives
# d(t), w(t), so we have to shift to get the "usual" faithful data
waiting = oldfaithful[1:, 1:]
duration = oldfaithful[:len(waiting), :1]
dt = N.concatenate((duration, waiting), 1)

# Scale the data so that each component is in [0..1]
dt = utils.scale(dt)

# This function train a mixture model with k components, returns the trained
# model and the BIC
def cluster(data, k, mode = 'full'):
    d = data.shape[1]
    gm = pyem.GM(d, k, mode)
    gmm = pyem.GMM(gm)
    em = pyem.EM()
    em.train(data, gmm, maxiter = 20)
    return gm, gmm.bic(data)

# bc will contain a list of BIC values for each model trained
bc = []
mode = 'full'
for k in range(1, 5):
    # Train a model of k component, and plots isodensity curve
    P.subplot(2, 2, k)
    gm, b = cluster(dt, k = k, mode = mode)
    bc.append(b)

    X, Y, Z, V = gm.density_on_grid()
    P.contour(X, Y, Z, V)
    P.plot(dt[:, 0], dt[:, 1], '.')
    P.xlabel('duration time (scaled)')
    P.ylabel('waiting time (scaled)')

print "According to the BIC, model with %d components is better" % (N.argmax(bc) + 1)
P.show()
