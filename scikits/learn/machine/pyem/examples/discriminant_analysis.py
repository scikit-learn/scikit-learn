#! /usr/bin/env python
# Last Change: Mon Jul 02 09:00 PM 2007 J

__doc__ = """Example of doing classification with mixture of Gaussian. Note
that this is really a toy example: we do not use testing testset nor cross
validation.

We use the famous iris database used by Sir R.A. Fisher. You can try to change
the attributes used for classification, number of components used for the
mixtures, etc..."""

import numpy as N
import pylab as P
import matplotlib as MPL

from scipy.sandbox import pyem
import utils

data = utils.iris.load()
# cnames are the class names
cnames = data.keys()

#--------------------
# Data pre processing
#--------------------
# we use 25 samples of each class (eg half of iris), for
# learning, and the other half for testing. We use sepal width and petal width
# only
ln = 25
tn = 25
xdata = {}
ydata = {}
# learning data
ldata = {}

# you can change here the used attributes (sepal vs petal, width vs height)
for i in cnames:
    xdata[i] = data[i]['sepal width']
    ydata[i] = data[i]['petal width']
    ldata[i] = N.concatenate((xdata[i][:ln, N.newaxis], 
                              ydata[i][:ln, N.newaxis]), 
                             axis = 1)

tx = N.concatenate([xdata[i][ln:] for i in cnames])
ty = N.concatenate([ydata[i][ln:] for i in cnames])
tdata = N.concatenate((tx[:, N.newaxis], ty[:, N.newaxis]), axis = 1)

# tclass makes the correspondance class <-> index in the testing data tdata
tclass = {}
for i in range(3):
    tclass[cnames[i]] = N.arange(tn * i, tn * (i+1))

#----------------------------
# Learning and classification
#----------------------------
# This function train a mixture model with k components
def cluster(data, k, mode = 'full'):
    d = data.shape[1]
    gm = pyem.GM(d, k, mode)
    gmm = pyem.GMM(gm)
    em = pyem.EM()
    em.train(data, gmm, maxiter = 20)
    return gm

# Estimate each class with a mixture of nc components
nc = 2
mode = 'diag'
lmod = {}
for i in cnames:
    lmod[i] = cluster(ldata[i], nc, mode)

# Classifiy the testing data. Of course, the data are not really IID, because
# we did not shuffle the testing data, but in this case, this does not change
# anything.
p = N.empty((len(tdata), 3))
for i in range(3):
    # For each class, computes the likelihood for the testing data
    p[:, i] = lmod[cnames[i]].pdf(tdata)

# We then take the Maximum A Posteriori class (same than most likely model in
# this case, since each class is equiprobable)
cid = N.argmax(p, 1)
classification = {}
for i in range(3):
    classification[cnames[i]] = N.where(cid == i)[0]

correct = {}
incorrect = {}
for i in cnames:
    correct[i] = N.intersect1d(classification[i], tclass[i])
    incorrect[i] = N.setdiff1d(classification[i], tclass[i])

#-----------------
# Plot the results
#-----------------
csym = {'setosa' : 's', 'versicolor' : 'x', 'virginica' : 'o'}
r = 50.
P.figure(figsize = [600/r, 400/r])

prop = MPL.font_manager.FontProperties(size='smaller')

# Plot the learning data with the mixtures
P.subplot(2, 1, 1)
for i in lmod.values():
    #i.plot()
    X, Y, Z, V = i.density_on_grid()
    P.contourf(X, Y, Z, V)

for i in cnames:
    P.plot(ldata[i][:, 0], ldata[i][:, 1], csym[i], label = i + ' (learning)')
P.xlabel('sepal width')
P.ylabel('petal width')
P.legend(loc = 'best')

# Plot the results on test dataset (green for correctly classified, red for
# incorrectly classified)
P.subplot(2, 1, 2)
for i in cnames:
    P.plot(tx[correct[i]], ty[correct[i]], 'g' + csym[i], 
           label = '%s (correctly classified)' % i)
    if len(incorrect[i]) > 0:
        P.plot(tx[incorrect[i]], ty[incorrect[i]], 'r' + csym[i], 
               label = '%s (incorrectly classified)' % i)
P.legend(loc = 'best', prop = prop)
P.xlabel('sepal width')
P.ylabel('petal width')
P.savefig('dclass.png', dpi = 60)
