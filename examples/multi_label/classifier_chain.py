"""
============================
Classifier Chain
============================

Demonstrate that a single classifier chain
out performs 10 independent classifiers
and that an ensemble of 10 randomly ordered
classifier chains out performs a single classifier chain

"""

import numpy as np
from sklearn.datasets import make_multilabel_classification
from sklearn.multi_label import ClassifierChain
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score

n_classes = 10
X, Y = make_multilabel_classification(n_samples=10000,
                                      n_features=100,
                                      n_classes=n_classes)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y)

# 1. fit an independent logistic regression model for each class
independent_models = [LogisticRegression() for _ in range(n_classes)]
Y_pred_ind = np.zeros(Y_test.shape)
for idx, independent_model in enumerate(independent_models):
    independent_model.fit(X_train, Y_train[:, idx])
    Y_pred_ind[:, idx] = independent_model.predict(X_test)
print("independent models average precision",
      average_precision_score(Y_test, Y_pred_ind, average='weighted'))

# 2. fit a single chain of logistic regression models
classifier_chain = ClassifierChain(LogisticRegression())
classifier_chain.fit(X_train, Y_train)
Y_pred = classifier_chain.predict(X_test)
print("single chain average precision",
      average_precision_score(Y_test, Y_pred, average='weighted'))

# 3. fit an ensemble of classifier chains
chains = [ClassifierChain(LogisticRegression()) for i in range(10)]
for chain in chains:
    chain.fit(X_train, Y_train)
Y_pred_mean = np.array([chain.predict(X_test)
                        for chain in chains]).mean(axis=0)
print("chain ensemble average precision",
      average_precision_score(Y_test, Y_pred_mean, average='weighted'))

# The three above cases should have increasing average precision
