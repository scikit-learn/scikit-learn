#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
====================================================
Pipelining: chaining a RBM and a logistic regression
====================================================

The BernoulliRBM does unsupervised feature extraction, while the logistic
regression does the prediction.

We use a GridSearchCV to set the number of hidden units and the learning rate
of the Bernoulli Restricted Boltzmann Machine.

We also train a simple logistic regression for comparison. The example shows
that the features extracted by the BernoulliRBM help improve the classification
accuracy.

"""
print __doc__


# Code source: Yann N. Dauphin
# License: BSD

import numpy as np
import pylab as pl

from sklearn import linear_model, datasets, metrics
from sklearn.cross_validation import train_test_split
from sklearn.neural_network import BernoulliRBM
from sklearn.pipeline import Pipeline

###############################################################################
# Setting up

np.random.seed(0xfeeb)

# Load Data
digits = datasets.load_digits()
X = np.asarray(digits.data, 'float32') / digits.data.max()
Y = digits.target
X_train, X_test, Y_train, Y_test = train_test_split(X, Y,
                                                    test_size=0.2,
                                                    random_state=0xfeeb)

# Models we will use
logistic = linear_model.LogisticRegression()
rbm = BernoulliRBM()

classifier = Pipeline(steps=[('rbm', rbm), ('logistic', logistic)])

###############################################################################
# Training

# Hyper-parameters. These were set by cross-validation,
# using a GridSearchCV. Here we are not performing cross-validation to
# save time.
rbm.learning_rate = 0.2
rbm.n_iter = 30
# More components tend to give better prediction performance, but larger
# fitting time
rbm.n_components = 500
logistic.C = 1e4

# Training RBM-Logistic Pipeline
classifier.fit(X_train, Y_train)

# Training Logistic regression
logistic_classifier = linear_model.LogisticRegression(C=1e4)
logistic_classifier.fit(X_train, Y_train)

###############################################################################
# Evaluation

print "Classification report for classifier %s:\n%s\n" % (
    classifier, metrics.classification_report(
        Y_test,
        classifier.predict(X_test)))

print "Classification report for classifier %s:\n%s\n" % (
    logistic_classifier, metrics.classification_report(
        Y_test,
        logistic_classifier.predict(X_test)))

###############################################################################
# Plotting

pl.figure(figsize=(4.2, 4))
for i, comp in enumerate(rbm.components_[:100]):
    pl.subplot(10, 10, i + 1)
    pl.imshow(comp.reshape((8, 8)), cmap=pl.cm.gray_r,
              interpolation='nearest')
    pl.xticks(())
    pl.yticks(())
pl.suptitle('100 components extracted by RBM', fontsize=16)
pl.subplots_adjust(0.08, 0.02, 0.92, 0.85, 0.08, 0.23)

pl.show()
