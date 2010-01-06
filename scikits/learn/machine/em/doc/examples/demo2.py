#! /usr/bin/env python

# Example of use of pyem toolbox. Feel free to change parameters
# such as dimension, number of components, mode of covariance.
#
# You can also try less trivial things such as adding outliers, sampling
# a mixture with full covariance and estimating it with a mixture with diagonal
# gaussians (replace the mode of the learned model lgm)
#
# Later, I hope to add functions for number of component estimation using eg BIC

import numpy as N
from numpy.random import seed

from scipy.sandbox.pyem import GM, GMM, EM
import copy

seed(2)
#+++++++++++++++++++++++++++++
# Meta parameters of the model
#   - k: Number of components
#   - d: dimension of each Gaussian
#   - mode: Mode of covariance matrix: full or diag (string)
#   - nframes: number of frames (frame = one data point = one
#   row of d elements)
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

lgm     = []
kmax    = 6
bics    = N.zeros(kmax)
for i in range(kmax):
    # Init the model with an empty Gaussian Mixture, and create a Gaussian 
    # Mixture Model from it
    lgm.append(GM(d, i+1, mode))
    gmm = GMM(lgm[i], 'kmean')

    # The actual EM, with likelihood computation. The threshold
    # is compared to the (linearly appromixated) derivative of the likelihood
    em      = EM()
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
    # Level is the confidence level for confidence ellipsoids: 1.0 means that
    # all points will be (almost surely) inside the ellipsoid
    level   = 0.8
    if not d == 1:
        P.plot(data[:, 0], data[:, 1], '.', label = '_nolegend_')

        # h keeps the handles of the plot, so that you can modify 
        # its parameters like label or color
        h   = lgm[k].plot(level = level)
        [i.set_color('r') for i in h]
        h[0].set_label('EM confidence ellipsoides')

        h   = gm.plot(level = level)
        [i.set_color('g') for i in h]
        h[0].set_label('Real confidence ellipsoides')
    else:
        # The 1d plotting function is quite elaborate: the confidence
        # interval are represented by filled areas, the pdf of the mixture and
        # the pdf of each component is drawn (optional)
        h   = gm.plot1d(level = level)
        [i.set_color('g') for i in h['pdf']]
        h['pdf'][0].set_label('true pdf')

        h0  = gm0.plot1d(level = level)
        [i.set_color('k') for i in h0['pdf']]
        h0['pdf'][0].set_label('initial pdf')

        hl  = lgm.plot1d(fill = 1, level = level)
        [i.set_color('r') for i in hl['pdf']]
        hl['pdf'][0].set_label('pdf found by EM')

        P.legend(loc = 0)

P.legend(loc = 0)
P.show()
# P.save('2d diag.png')
