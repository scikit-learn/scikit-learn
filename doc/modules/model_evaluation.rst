.. _model_evaluation:

===================
Model evaluation
===================

.. TODO

  Metrics
  =======


Random classification accuracy
=================================

.. currentmodule:: sklearn.random_classifier

When doing classification, a simple sanity check consists in comparing one's
classifier against the accuracy of a purely random classifier.
:class:`RandomClassifier` implements three different strategies for random
classification. `stratified` generates predictions by respecting the training
set's class distribution. `most_frequent` always predicts the most frequent
label in the training set.  `uniform` generates predictions uniformly at
random. To illustrate :class:`RandomClassifier`, first let's create an imbalanced
dataset::

  >>> from sklearn.datasets import load_iris
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> y[y != 1] = -1


Next, let's compare the accuracy of `SVC` and `most_frequent`::

  >>> from sklearn.random_classifier import RandomClassifier
  >>> from sklearn.svm import SVC
  >>> clf = SVC(kernel='linear', C=1).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.73...
  >>> clf = RandomClassifier(sampling='most_frequent', random_state=0).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.66...

We see that `SVC` doesn't do much better than a dummy classifier. Now, let's change
the kernel::

  >>> clf = SVC(kernel='rbf', C=1).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.99...

More generally, when the accuracy of a classifier is too close to random classification, it
probably means that something went wrong: features are not helpful, a
hyparameter is not correctly tuned, the classifier is suffering from class
imbalance, etc...
