#! /usr/bin/env python
# Last Change: Mon Jul 02 06:00 PM 2007 J

__doc__ = """This example shows how to do pdfestimation for one dimensional
data. It estimates a Gaussian mixture model for several number of components,
and it determines the "best" one according to the Bayesian Information
Criterion.

It uses old faitfhul waiting time as the one dimension data, and plots the best
model as well as the BIC as a function of the number of component."""
# Example of doing pdf estimation with EM algorithm. Requires matplotlib.
import numpy as N
from numpy.testing import set_package_path, restore_path

import pylab as P
import matplotlib as MPL

set_package_path()
import pyem
restore_path()
import utils

oldfaithful = utils.get_faithful()

duration = oldfaithful[:, :1]
waiting = oldfaithful[:, 1:]

#dt = utils.scale(duration)
#dt = duration / 60.
dt = waiting / 60.

# This function train a mixture model with k components, returns the trained
# model and the BIC
def cluster(data, k):
    d = data.shape[1]
    gm = pyem.GM(d, k)
    gmm = pyem.GMM(gm)
    em = pyem.EM()
    em.train(data, gmm, maxiter = 20)
    return gm, gmm.bic(data)

# bc will contain a list of BIC values for each model trained, gml the
# corresponding mixture model
bc = []
gml = []

for k in range(1, 8):
    gm, b = cluster(dt, k = k)
    bc.append(b)
    gml.append(gm)

mbic = N.argmax(bc)

# Below is code to display a figure with histogram and best model (in the BIC
# sense) pdf, with the BIC as a function of the number of components on the
# right.
P.figure(figsize = [12, 7])
#---------------------------
# histogram + pdf estimation
#---------------------------
P.subplot(1, 2, 1)
h = gml[mbic].plot1d(gpdf=True)
# You can manipulate the differents parts of the plot through the returned
# handles
h['gpdf'][0].set_linestyle('-')
h['gpdf'][0].set_label('pdf of the mixture')
h['pdf'][0].set_label('pdf of individual component')
[l.set_linestyle('-') for l in h['pdf']]
[l.set_color('g') for l in h['pdf']]

prop = MPL.font_manager.FontProperties(size='smaller')
P.legend(loc = 'best', prop = prop)

P.hist(dt, 25, normed = 1, fill = False)
P.xlabel('waiting time between consecutive eruption (in min)')

#------------------------------------------
# BIC as a function of number of components
#------------------------------------------
P.subplot(1, 2, 2)
P.plot(N.arange(1, 8), bc, 'o:')
P.xlabel("number of components")
P.ylabel("BIC")
print "According to the BIC, model with %d components is better" % (N.argmax(bc) + 1)
P.show()
