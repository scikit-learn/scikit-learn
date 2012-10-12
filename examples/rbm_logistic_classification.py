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

Note
----

Much better performance can be achieved by using larger n_components and n_iter
for the RestrictedBolzmannMachine.

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

# Load Data
digits = datasets.fetch_mldata('MNIST original')
X = np.asarray(digits.data, 'float32') / digits.data.max()
Y = digits.target == 8 # Classification of class 8 vs all

# Models we will use
logistic = linear_model.LogisticRegression()
rbm = RestrictedBolzmannMachine()

pipe = Pipeline(steps=[('rbm', rbm), ('logistic', logistic)])
logistic_pipe = Pipeline(steps=[('logistic', logistic)])

###############################################################################
# Training

# Hyper-parameters
n_components = [256]
learning_rate = [0.01]
Cs = [10000]
n_iter = [5]

# Training RBM-Logistic Pipeline
estimator = GridSearchCV(pipe,
                         dict(rbm__n_components=n_components,
                              rbm__learning_rate=learning_rate,
                              rbm__n_iter=n_iter,
                              logistic__C=Cs))
estimator.fit(X[:60000], Y[:60000])
classifier = estimator.best_estimator_

# Training Logistic regression
logistic_estimator = GridSearchCV(logistic_pipe,
                                  dict(logistic__C=np.logspace(-4, 4, 2)))
logistic_estimator.fit(X[:60000], Y[:60000])
logistic_classifier = logistic_estimator.best_estimator_


###############################################################################
# Evaluation

print "Classification report for classifier %s:\n%s\n" % (
    classifier, metrics.classification_report(Y[60000:],
        classifier.predict(X[60000:])))

print "Classification report for classifier %s:\n%s\n" % (
    logistic_classifier, metrics.classification_report(Y[60000:],
        logistic_classifier.predict(X[60000:])))
