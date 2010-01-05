#! /usr/bin/env python
# Last Change: Mon Jul 09 07:00 PM 2007 J

__doc__ = """Example of doing classification with mixture of Gaussian. Note
that this is really a toy example: we do not use testing testset nor cross
validation.

We use the famous iris database used by Sir R.A. Fisher. You can try to change
the attributes used for classification, number of components used for the
mixtures, etc..."""

import numpy as N
import pylab as P
import matplotlib as MPL

from scipy.sandbox import svm
import utils

from scikits.learn.datasets import german
data = german.load()

features = N.vstack([data['feat']['feat' + str(i)].astype(N.float) for i in range(1, 25)]).T
label = data['label']

t, s = utils.scale(features)

training = svm.LibSvmClassificationDataSet(label, features)

def train_svm(cost, gamma, fold = 5):
    """Train a SVM for given cost and gamma."""
    kernel = svm.RBFKernel(gamma = gamma)
    model = svm.LibSvmCClassificationModel(kernel, cost = cost)
    cv = model.cross_validate(training, fold)
    return cv

c_range = N.exp(N.log(2.) * N.arange(-5, 15))
g_range = N.exp(N.log(2.) * N.arange(-15, 3))

# Train the svm on a log distributed grid
gr = N.meshgrid(c_range, g_range)
c = gr[0].flatten()
g = gr[1].flatten()
cf = N.hstack((c, g))
cv = N.empty(c.size)
for i in range(cv.size):
    print "=============== iteration %d / %d ============" % (i, cv.size)
    cv[i] = train_svm(c[i], g[i])

v = P.contourf(N.log2(gr[0]), N.log2(gr[1]), cv.reshape(g_range.size, c_range.size), 10)
v = P.contour(N.log2(gr[0]), N.log2(gr[1]), cv.reshape(g_range.size, c_range.size), 10)
P.clabel(v, inline = 1, fontsize = 10)
P.show()
