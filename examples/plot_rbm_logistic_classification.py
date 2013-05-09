#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
===============================================================
Restricted Boltzmann Machine features for digits classification
===============================================================

For greyscale image data where pixel values can be interpreted as degrees of
blackness on a white background, like handwritten digit recognition, the
Bernoulli Restricted Boltzmann machine model (:class:`BernoulliRBM
<sklearn.neural_network.BernoulliRBM>`) can perform effective non-linear
feature extraction.

This example shows how to build a classification pipeline with a BernoulliRBM
feature extractor and a :class:`LogisticRegression
<sklearn.linear_model.LogisticRegression>` classifier.  The hyperparameters
of the entire model (learning rate, hidden layer size, regularization)
were optimized by grid search, but the search is not reproduced here because
of runtime constraints.

Logistic regression on raw pixel values is presented for comparison. The
example shows that the features extracted by the BernoulliRBM help improve the
classification accuracy.
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

# Load Data
digits = datasets.load_digits()
X = np.asarray(digits.data, 'float32') / digits.data.max()
Y = digits.target
X_train, X_test, Y_train, Y_test = train_test_split(X, Y,
                                                    test_size=0.2,
                                                    random_state=0)

# Models we will use
logistic = linear_model.LogisticRegression()
rbm = BernoulliRBM(random_state=0)

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
