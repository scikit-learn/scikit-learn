"""
============================================================
Classifier Chain algorithm for multi-label classification
============================================================

This examples shows how to use Classifier Chain algorithm to deal with
multi-label classification, which is done using the
:class:`sklearn.multi_label.ClassifierChain` object.

"""

from __future__ import print_function

from sklearn import datasets
from sklearn import metrics
from sklearn.multi_label import ClassifierChain
from sklearn.cross_validation import train_test_split
from sklearn.svm import LinearSVC

print(__doc__)

# Loading the sythetic dataset
X, Y = datasets.make_multilabel_classification(return_indicator=True)

# Split the dataset in two equal parts
X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.5, random_state=0)

clf = ClassifierChain(base_estimator=LinearSVC())
clf.fit(X_train, Y_train)

Y_true, Y_pred = Y_train, clf.predict(X_train)

print("Hamming Loss (Training):")
print(metrics.hamming_loss(Y_true, Y_pred))

Y_true, Y_pred = Y_test, clf.predict(X_test)

print()

print("Hamming Loss (Training):")
print(metrics.hamming_loss(Y_true, Y_pred))
