"""
============================
Classifier Chain
============================
An ensemble of 10 logistic regression classifier chains trained on a
multi-label dataset achieves a higher Jaccard similarity score than a set
of independently trained logistic regression models.

"""

import numpy as np
from sklearn.multi_label import ClassifierChain
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import jaccard_similarity_score
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import fetch_mldata

# Load a multi-label dataset
yeast = fetch_mldata('yeast')
X = yeast['data']
Y = yeast['target'].transpose().toarray()
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=.5)

# Fit an independent logistic regression model for each class using the
# OneVsRestClassifier wrapper
ovr = OneVsRestClassifier(LogisticRegression())
ovr.fit(X_train, Y_train)
Y_pred_ovr = ovr.predict(X_test)
print("Independent models Jaccard similarity score:", jaccard_similarity_score(
    Y_test, Y_pred_ovr))

# Fit an ensemble of logistic regression classifier chains and take the
# take the average prediction of all the chains
chains = [ClassifierChain(LogisticRegression()) for i in range(10)]
for chain in chains:
    chain.fit(X_train, Y_train)
Y_pred_ensemble = np.array([chain.predict(X_test) for chain in
                            chains]).mean(axis=0)
print("Classifier chain ensemble Jaccard similarity score:",
      jaccard_similarity_score(Y_test, Y_pred_ensemble >= .5))

