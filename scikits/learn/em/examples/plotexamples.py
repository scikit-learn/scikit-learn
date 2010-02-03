#! /usr/bin/env python
# Last Change: Sun Jul 22 12:00 PM 2007 J

# This is a simple test to check whether plotting ellipsoides of confidence and
# isodensity contours match
import numpy as N

import pylab as P

from scikits.learn.machine.em import EM, GM, GMM

# Generate a simple mixture model, plot its confidence ellipses + isodensity
# curves for both diagonal and full covariance matrices
d = 3
k = 3
dim = [0, 2]
# diag model
w, mu, va = GM.gen_param(d, k)
dgm = GM.fromvalues(w, mu, va)
# full model
w, mu, va = GM.gen_param(d, k, 'full', spread = 1)
fgm = GM.fromvalues(w, mu, va)

def plot_model(gm, dim):
    X, Y, Z, V = gm.density_on_grid(dim = dim)
    h = gm.plot(dim = dim)
    [i.set_linestyle('-.') for i in h]
    P.contour(X, Y, Z, V)
    data = gm.sample(200)
    P.plot(data[:, dim[0]], data[:,dim[1]], '.')

# Plot the contours and the ellipsoids of confidence
P.subplot(2, 1, 1)
plot_model(dgm, dim)

P.subplot(2, 1, 2)
plot_model(fgm, dim)

P.show()
