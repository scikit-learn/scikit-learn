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

seed(1)
#+++++++++++++++++++++++++++++
# Meta parameters of the model
#   - k: Number of components
#   - d: dimension of each Gaussian
#   - mode: Mode of covariance matrix: full or diag (string)
#   - nframes: number of frames (frame = one data point = one
#   row of d elements)
k       = 2 
d       = 2
mode    = 'diag'
nframes = 1e3

#+++++++++++++++++++++++++++++++++++++++++++
# Create an artificial GMM model, samples it
#+++++++++++++++++++++++++++++++++++++++++++
w, mu, va   = GM.gen_param(d, k, mode, spread = 1.5)
gm          = GM.fromvalues(w, mu, va)

# Sample nframes frames  from the model
data    = gm.sample(nframes)

#++++++++++++++++++++++++
# Learn the model with EM
#++++++++++++++++++++++++

# Init the model
lgm = GM(d, k, mode)
gmm = GMM(lgm, 'kmean')
gmm.init(data)

# Keep a copy for drawing later
gm0 = copy.copy(lgm)

# The actual EM, with likelihood computation. The treshold
# is compared to the (linearly appromixated) derivative of the likelihood
em      = EM()
like    = em.train(data, gmm, maxiter = 30, thresh = 1e-8)

#+++++++++++++++
# Draw the model
#+++++++++++++++
import pylab as P
P.subplot(2, 1, 1)

level   = 0.8
if not d == 1:
    P.plot(data[:, 0], data[:, 1], '.', label = '_nolegend_')

    # h keeps the handles of the plot, so that you can modify 
    # its parameters like label or color
    h   = gm.plot(level = level)
    [i.set_color('g') for i in h]
    h[0].set_label('true confidence ellipsoides')

    # Initial confidence ellipses as found by kmean
    h   = gm0.plot(level = level)
    [i.set_color('k') for i in h]
    h[0].set_label('kmean confidence ellipsoides')

    # Values found by EM
    h   = lgm.plot(level = level)
    [i.set_color('r') for i in h]
    h[0].set_label('kmean confidence ellipsoides')

    P.legend(loc = 0)
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

P.subplot(2, 1, 2)
P.plot(like)
P.title('log likelihood')

P.show()
P.save('2d diag.png')
