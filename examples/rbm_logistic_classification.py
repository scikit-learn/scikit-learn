#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Pipelining: chaining a RBM and a logistic regression
=========================================================

The RestrictedBolzmannMachine does unsupervised feature extraction,
while the logistic regression does the prediction.

We use a GridSearchCV to set the number of hidden units and the learning rate
of the RestrictedBolzmannMachine.

We also train a simple logistic regression for comparison. The example shows
that the features extracted by the RestrictedBolzmannMachine help improve
the classification accuracy.

"""
print __doc__


# Code source: Yann N. Dauphin
# License: BSD


import numpy as np

from sklearn import linear_model, datasets, metrics, preprocessing
from sklearn.rbm import RestrictedBolzmannMachine
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV

###############################################################################
# Setting up

np.random.seed(0xfeeb)

# Load Data
digits = datasets.load_digits()
X = np.asarray(digits.data, 'float32') / digits.data.max()
Y = digits.target
train_size = int(0.8 * X.shape[0])

# Shuffle train-test split
inds = range(X.shape[0])
np.random.shuffle(inds)
X = X[inds]
Y = Y[inds]

# Models we will use
logistic = linear_model.LogisticRegression()
rbm = RestrictedBolzmannMachine()

pipe = Pipeline(steps=[('rbm', rbm), ('logistic', logistic)])

###############################################################################
# Training

# Hyper-parameters
n_components = [350, 400, 450]
learning_rate = [0.01, 0.01]
Cs = np.logspace(-4, 4, 2)
n_iter = [20, 30]

# Training RBM-Logistic Pipeline
estimator = GridSearchCV(pipe,
                         dict(rbm__n_components=n_components,
                              rbm__learning_rate=learning_rate,
                              rbm__n_iter=n_iter,
                              logistic__C=Cs))
estimator.fit(X[:train_size], Y[:train_size])
classifier = estimator.best_estimator_

# Training Logistic regression
logistic_estimator = GridSearchCV(logistic,
                                  dict(C=Cs))
logistic_estimator.fit(X[:train_size], Y[:train_size])
logistic_classifier = logistic_estimator.best_estimator_


###############################################################################
# Evaluation

print "Classification report for classifier %s:\n%s\n" % (
    classifier, metrics.classification_report(Y[train_size:],
        classifier.predict(X[train_size:])))

print "Classification report for classifier %s:\n%s\n" % (
    logistic_classifier, metrics.classification_report(Y[train_size:],
        logistic_classifier.predict(X[train_size:])))
