import numpy as N
from numpy.random import seed

from scipy.sandbox.pyem import GM, GMM, EM
import copy

seed(2)

k       = 4
d       = 2
mode    = 'diag'
nframes = 1e3

#+++++++++++++++++++++++++++++++++++++++++++
# Create an artificial GMM model, samples it
#+++++++++++++++++++++++++++++++++++++++++++
w, mu, va   = GM.gen_param(d, k, mode, spread = 1.0)
gm          = GM.fromvalues(w, mu, va)

# Sample nframes frames  from the model
data    = gm.sample(nframes)

#++++++++++++++++++++++++
# Learn the model with EM
#++++++++++++++++++++++++

# List of learned mixtures lgm[i] is a mixture with i+1 components
lgm     = []
kmax    = 6
bics    = N.zeros(kmax)
em      = EM()
for i in range(kmax):
    lgm.append(GM(d, i+1, mode))

    gmm = GMM(lgm[i], 'kmean')
    em.train(data, gmm, maxiter = 30, thresh = 1e-10)
    bics[i] = gmm.bic(data)

print "Original model has %d clusters, bics says %d" % (k, N.argmax(bics)+1)

#+++++++++++++++
# Draw the model
#+++++++++++++++
import pylab as P
P.subplot(3, 2, 1)

for k in range(kmax):
    P.subplot(3, 2, k+1)
    level   = 0.9
    P.plot(data[:, 0], data[:, 1], '.', label = '_nolegend_')

    # h keeps the handles of the plot, so that you can modify 
    # its parameters like label or color
    h   = lgm[k].plot(level = level)
    [i.set_color('r') for i in h]
    h[0].set_label('EM confidence ellipsoides')

    h   = gm.plot(level = level)
    [i.set_color('g') for i in h]
    h[0].set_label('Real confidence ellipsoides')

P.legend(loc = 0)
# depending on your configuration, you may have to call P.show() 
# to actually display the figure
